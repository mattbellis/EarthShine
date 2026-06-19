"""
generate_signal_events_CL.py  (v2)
==================================
Same CLI as before, with the storage layer modernized:

  * shards are still produced by dgt.generate_many_events() exactly as before
  * the merged (COMBINED) output is written into the partitioned tree under
    --data-root with full generation params embedded in the parquet footer
    (see earthshine_io.py); the catalog is rebuilt at the end
  * fixes from review:
      - additional_tag tuple bug (trailing comma) removed
      - shard glob moved out of the trial loop (was re-globbing every trial)
      - total_org_nevents kept as a column for backward compatibility, but
        the canonical copy now lives in footer metadata as n_generated

Legacy shard files in --output-directory are unchanged; delete them after
the merge if you don't need them (or keep them and run migrate on them too).
"""

import argparse
import os
import sys

############################################################################
# Redirect OS-level fd 2 to suppress C++ noise (CUDA, XLA, absl)
# but keep a Python stderr writer so exceptions still show
os.environ['CUDA_VISIBLE_DEVICES'] = '-1'
_real_stderr = os.dup(2)
os.dup2(os.open(os.devnull, os.O_WRONLY), 2)
sys.stderr = os.fdopen(_real_stderr, "w", closefd=False)
############################################################################

import datetime
import glob
import re
import logging
import uuid

import numpy as np
import pandas as pd

import dm_generation_tools as dgt
import earthshine_io as eio


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate dark matter scattering events."
    )
    parser.add_argument("--masses-a", type=float, nargs='+', default=[1])
    parser.add_argument("--masses-dm", type=float, nargs='+', default=[100])
    parser.add_argument("--nevents", type=int, default=10)
    parser.add_argument("--ntrials", type=int, default=3)
    parser.add_argument("--detector-radius", type=float, default=7.5)
    parser.add_argument("--detector-half-len", type=float, default=15)
    parser.add_argument("--inner-detector-radius", type=float, default=1.0)
    parser.add_argument("--inner-detector-half-len", type=float, default=2.5)
    parser.add_argument("--disk-radius", type=float, default=None)
    parser.add_argument("--depth", type=float, nargs='+', default=[])
    parser.add_argument("--dm-model", type=str, default="floating")
    parser.add_argument("--save-only-detected",
                        action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--eloss-dict-name", type=str, default=None)
    parser.add_argument("--save-to-file", action="store_true", default=False)
    parser.add_argument("--additional-tag", type=str, default="")
    parser.add_argument("--output-directory", type=str, default=None,
                        help="Scratch dir for intermediate shard files")
    parser.add_argument("--data-root", type=str, default="data",
                        help="Root of the partitioned output tree "
                             "(combined file + catalog go here)")
    parser.add_argument("--average-eloss", action="store_true", default=False)
    parser.add_argument("--seed", type=int, default=None,
                        help="Base RNG seed; trial i uses seed+i (numpy "
                             "and tensorflow).  Default: auto-generate a "
                             "unique seed and record it in the metadata")
    parser.add_argument("--merge-run-id", type=str, default=None,
                        help="Skip generation; merge the existing shards "
                             "with this run_id from --output-directory "
                             "(pass the SAME physics args and --nevents "
                             "as the original run)")
    parser.add_argument("--verbose", action=argparse.BooleanOptionalAction,
                        default=False)
    return parser.parse_args()


################################################################################
def main():
    logging.getLogger("tensorflow").setLevel(logging.FATAL)

    args = parse_args()

    if len(args.depth) == 1:
        args.depth = args.depth[0]

    # Unique id for THIS invocation.  It is woven into the shard tag, which
    # (a) prevents a second run at the same parameters from overwriting the
    # first run's shards in scratch, and (b) makes the merge glob below
    # match only this run's shards -- never stale ones from earlier runs.
    # It also goes into the combined file's metadata and filename, so
    # repeated runs at the same physics point coexist as separate datasets
    # (read them together with eio.load_many()).
    # timestamp + PID + random hex: PIDs are unique among simultaneously
    # running processes, so concurrent launches from multiple terminals
    # can never collide even within the same second
    run_id = (datetime.datetime.now().strftime("%Y%m%dT%H%M%S")
              + f"p{os.getpid()}" + uuid.uuid4().hex[:4])
    if args.merge_run_id:
        run_id = args.merge_run_id
        print(f"merge-only mode for run_id: {run_id}")
    else:
        print(f"run_id: {run_id}")

    gen_disk_radius = None   # set by generate_many_events if it returns it

    if not args.merge_run_id:
        # Auto-seed if none given: every run gets a unique, RECORDED seed,
        # so runs never duplicate each other's events and any run can be
        # (best-effort) replayed from its catalog entry.  Trial i uses
        # seed+i.  Note: NumPy seeding alone does not control phasespace,
        # which draws from TensorFlow's RNG -- we seed both.  Exact TF
        # replay can still depend on TF version/hardware (GPU reductions
        # may be nondeterministic), so treat the seed as best-effort
        # replay, not a bit-exact guarantee.
        if args.seed is None:
            args.seed = int(np.random.SeedSequence().entropy % 2**31)
            print(f"auto-generated seed: {args.seed}")

        try:
            import tensorflow as tf
            _tf_set_seed = tf.random.set_seed
        except Exception:
            _tf_set_seed = None
            print("WARNING: could not import tensorflow to seed it; "
                  "phasespace draws will not be reproducible")

        for i in range(args.ntrials):
            print(f'\nTrial {i} out of {args.ntrials} '
                  f'--------------------------------------------------------')
            print(f'Generate {args.nevents} events...')

            np.random.seed(args.seed + i)
            if _tf_set_seed is not None:
                _tf_set_seed(args.seed + i)

            additional_tag = f'_{run_id}_{i:06d}'
            if args.additional_tag:
                additional_tag = f'_{args.additional_tag}_{run_id}_{i:06d}'

            result = dgt.generate_many_events(
                MASSES_A=args.masses_a,
                MASSES_DM=args.masses_dm,
                nevents_to_generate=args.nevents,
                detector_radius=args.detector_radius,
                detector_half_len=args.detector_half_len,
                inner_detector_radius=args.inner_detector_radius,
                inner_detector_half_len=args.inner_detector_half_len,
                disk_radius=args.disk_radius,
                depth=args.depth,
                dm_model=args.dm_model,
                SAVE_ONLY_DETECTED=args.save_only_detected,
                eloss_dict_name=args.eloss_dict_name,
                save_to_file=args.save_to_file,
                additional_tag=additional_tag,
                output_directory=args.output_directory,
                AVERAGE_ELOSS=args.average_eloss,
                verbose=args.verbose,
            )
            # newer dm_generation_tools returns (df, tag, disk_radius);
            # older versions return (df, tag) -- support both
            if len(result) == 3:
                df_decays, tag, gen_disk_radius = result
            else:
                df_decays, tag = result
                gen_disk_radius = None

    ##############################################################
    # Merge the shard files  (once, after all trials)
    # The run_id is unique per invocation, so matching on it picks up
    # exactly this run's shards -- never stale ones from earlier runs.
    ##############################################################
    files = sorted(glob.glob(
        f'{args.output_directory}/generated_data_*_{run_id}_[0-9]*.parquet'))
    print(f'{len(files)} shard files matched for merging')
    if len(files) == 0:
        raise SystemExit(f"no shards found for run_id {run_id} in "
                         f"{args.output_directory}")
    if not args.merge_run_id and len(files) != args.ntrials:
        print(f'WARNING: expected {args.ntrials} shards but matched '
              f'{len(files)}.  Stale or missing shard files?  '
              f'Inspect before trusting n_generated.')

    # Empty shards (zero detected events in that trial) are legitimate at
    # low-acceptance parameter points.  They carry no rows, so org_nevents
    # cannot be read from them -- but they still represent args.nevents
    # THROWN events and must be counted, or n_generated is undercounted
    # and the acceptance biased upward.
    dfs, tot_nevents, n_empty = [], 0, 0
    for filename in files:
        df_tmp = pd.read_parquet(filename)
        if len(df_tmp):
            n_this = int(df_tmp['org_nevents'].iloc[0])
            if n_this != args.nevents:
                print(f'WARNING: {filename} has org_nevents={n_this} but '
                      f'--nevents={args.nevents}.  If merging an old run, '
                      f'pass the SAME --nevents it was generated with.')
        else:
            n_this = args.nevents
            n_empty += 1
        tot_nevents += n_this
        dfs.append(df_tmp)   # empty frames keep the schema; concat is fine

    if n_empty:
        print(f'{n_empty} shard(s) had zero detected events '
              f'(counted as {args.nevents} thrown each)')

    print(f'\n Merging {len(dfs)} intermediate files...\n')
    df_decays = pd.concat(dfs, ignore_index=True)
    if len(df_decays) == 0:
        print('NOTE: zero detected events in the entire run -- writing an '
              'empty combined file (n_generated metadata still records the '
              'exposure, which is what a limit needs).')
    # kept as a column for backward compatibility with existing notebooks;
    # the canonical value is n_generated in the footer metadata
    df_decays['total_org_nevents'] = int(tot_nevents)

    ##############################################################
    # Recover the disk radius when it was auto-computed.
    # If --disk-radius is not given, dgt computes one from angular
    # considerations.  Preferred source: the value returned directly by
    # generate_many_events (newer dm_generation_tools).  Fallback: parse
    # the 'diskR_<value>' token from the tag (older dgt) or from a shard
    # filename (merge-only mode, where dgt is never called).
    ##############################################################
    disk_radius = args.disk_radius
    disk_radius_auto = False
    if disk_radius is None:
        if gen_disk_radius is not None:
            disk_radius = float(gen_disk_radius)
            disk_radius_auto = True
            print(f'disk radius auto-computed by generator: {disk_radius} m '
                  f'(returned directly; recorded in metadata)')
        else:
            source = (tag if not args.merge_run_id
                      else os.path.basename(files[0]))
            m = re.search(r'diskR_(-?\d+(?:\.\d+)?)', source)
            if m:
                disk_radius = float(m.group(1))
                disk_radius_auto = True
                print(f'disk radius auto-computed by generator: '
                      f'{disk_radius} m (parsed from tag; recorded in '
                      f'metadata)')
            else:
                print(f'WARNING: --disk-radius not given and no diskR '
                      f'token found in {source!r}; disk_radius will be '
                      f'None in the metadata')

    ##############################################################
    # Write combined file into the partitioned tree, with metadata
    ##############################################################
    params = eio.make_params(
        dm_model=args.dm_model,
        masses_dm=args.masses_dm,
        masses_a=args.masses_a,
        depth=args.depth,
        disk_radius=disk_radius,
        selection="hit_detector" if args.save_only_detected else "all",
        eloss=("ave" if args.average_eloss
               else ("geant_spline" if args.eloss_dict_name else "none")),
        stage="combined",
        detector_radius=args.detector_radius,
        detector_half_len=args.detector_half_len,
        inner_detector_radius=args.inner_detector_radius,
        inner_detector_half_len=args.inner_detector_half_len,
        n_generated=tot_nevents,
        n_events_per_trial=args.nevents,
        n_trials=(len(files) if args.merge_run_id else args.ntrials),
        seed=args.seed,
        extra={"run_id": run_id,
               "disk_radius_auto": disk_radius_auto,
               "eloss_dict_name": args.eloss_dict_name,
               "additional_tag": args.additional_tag or None},
    )

    outpath = eio.write_parquet(
        df_decays, f'{args.data_root}/{eio.rel_path(params)}', params)
    print(f'Wrote {outpath}  ({len(df_decays)} detected events '
          f'from {tot_nevents} generated)')

    eio.build_catalog(args.data_root)

    del df_decays
    del dfs[:]


if __name__ == "__main__":
    main()
