#!/usr/bin/env python
import argparse
import logging
import os
import sys
import subprocess
from datetime import datetime


def main():
    parser = argparse.ArgumentParser(
        description='Submit generate_signal_events_CL.py jobs to SLURM'
    )
    parser.add_argument('-t', '--test', action='store_true',
                        help='Print the submission command but do not submit')

    # Physics parameters (forwarded to generate_signal_events_CL.py)
    parser.add_argument('--masses-a', type=float, nargs='+', required=True)
    parser.add_argument('--masses-dm', type=float, nargs='+', required=True)
    parser.add_argument('--nevents', type=int, default=1000)
    parser.add_argument('--ntrials', type=int, default=50)
    parser.add_argument('--depth', type=float, nargs='+', default=[])
    parser.add_argument('--dm-model', type=str, default='floating')
    parser.add_argument('--save-to-file', action='store_true', default=False)
    parser.add_argument('--output-directory', type=str, default='OUTPUT_FILES')
    parser.add_argument('--data-root', type=str, default='data',
                        help='Root of the partitioned output tree that the '
                             'combined file + catalog are written to (v2 '
                             'generator). Required for the diagnostics/'
                             'acceptance code, which reads via '
                             'eio.read_catalog(data_root).')
    parser.add_argument('--additional-tag', type=str, default='')
    parser.add_argument('--average-eloss', action='store_true', default=False)
    parser.add_argument('--save-only-detected', action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument('--eloss-dict-name', type=str, default=None)
    parser.add_argument('--detector-radius', type=float, default=7.5)
    parser.add_argument('--detector-half-len', type=float, default=15)
    parser.add_argument('--inner-detector-radius', type=float, default=1.0)
    parser.add_argument('--inner-detector-half-len', type=float, default=2.5)
    parser.add_argument('--disk-radius', type=float, default=None)
    parser.add_argument('--verbose', action=argparse.BooleanOptionalAction, default=False)

    # SLURM parameters
    parser.add_argument('--partition', type=str, default='batch',
                        help='SLURM partition (default: batch)')
    parser.add_argument('--mem', type=str, default='4000M',
                        help='Memory per job (default: 4000M)')
    parser.add_argument('--time', type=str, default='24:00:00',
                        help='Wall-clock time limit (default: 24:00:00)')
    parser.add_argument('--mail-user', type=str, default='tamasvami@ucsb.edu')

    # Parallelization parameters
    parser.add_argument('--parallel', action='store_true', default=False,
                        help='Run trials as a SLURM array job (one task per trial)')
    parser.add_argument('--base-seed', type=int, default=42,
                        help='Base seed for reproducibility; trial i gets seed = base_seed + i (default: 42)')
    parser.add_argument('--max-concurrent', type=int, default=None,
                        help='Max simultaneous array tasks (SLURM %%N syntax)')

    args = parser.parse_args()

    logging.basicConfig(
        format='[ submitSignal ][ %(levelname)s ]: %(message)s',
        level=logging.DEBUG
    )

    # Directories
    jobdir = '/home/vamitamas/slurm/jobs/'
    logdir = '/home/vamitamas/slurm/logs/'
    os.makedirs(jobdir, exist_ok=True)
    os.makedirs(logdir, exist_ok=True)

    # Where the script lives
    script_dir = os.path.dirname(os.path.abspath(__file__))
    venv_activate = '/home/vamitamas/EarthShine/venv/bin/activate'

    # Job name tag
    current_time = datetime.now().strftime("%Y%m%d_%H%M%S")
    tag = f'mA{"_".join(str(m) for m in args.masses_a)}_mDM{"_".join(str(m) for m in args.masses_dm)}'
    job_name = f'signal_{tag}_{current_time}'

    # Helper: build the list of optional physics flags shared by both modes
    def _optional_physics_flags():
        parts = []
        if args.depth:
            parts.append(f'--depth {" ".join(str(d) for d in args.depth)}')
        if args.disk_radius is not None:
            parts.append(f'--disk-radius {args.disk_radius}')
        if args.save_to_file:
            parts.append('--save-to-file')
        if args.output_directory:
            parts.append(f'--output-directory {args.output_directory}')
        if args.data_root:
            parts.append(f'--data-root {args.data_root}')
        if args.additional_tag:
            parts.append(f'--additional-tag {args.additional_tag}')
        if args.average_eloss:
            parts.append('--average-eloss')
        if not args.save_only_detected:
            parts.append('--no-save-only-detected')
        if args.eloss_dict_name is not None:
            parts.append(f'--eloss-dict-name {args.eloss_dict_name}')
        if args.verbose:
            parts.append('--verbose')
        return parts

    if not args.parallel:
        ################################################################
        # Sequential mode (original behavior)
        ################################################################
        cmd_parts = [
            f'python {script_dir}/generate_signal_events_CL_v2.py',
            f'--masses-a {" ".join(str(m) for m in args.masses_a)}',
            f'--masses-dm {" ".join(str(m) for m in args.masses_dm)}',
            f'--nevents {args.nevents}',
            f'--ntrials {args.ntrials}',
            f'--dm-model {args.dm_model}',
            f'--detector-radius {args.detector_radius}',
            f'--detector-half-len {args.detector_half_len}',
            f'--inner-detector-radius {args.inner_detector_radius}',
            f'--inner-detector-half-len {args.inner_detector_half_len}',
        ] + _optional_physics_flags()

        python_command = ' \\\n    '.join(cmd_parts)

        job_file = f'{jobdir}/slurm_signal_{job_name}.job'
        with open(job_file, 'w') as f:
            f.write('#!/bin/bash\n\n')
            f.write(f'#SBATCH --job-name={job_name}\n')
            f.write('#SBATCH --nodes=1 --ntasks-per-node=1\n')
            f.write('#SBATCH --ntasks=1\n')
            f.write('#SBATCH --cpus-per-task=1\n')
            f.write(f'#SBATCH --mem={args.mem}\n')
            f.write(f'#SBATCH --error={logdir}/slurm-%A_%a.err\n')
            f.write(f'#SBATCH --output={logdir}/slurm-%A_%a.out\n')
            f.write(f'#SBATCH --time={args.time}\n')
            f.write(f'#SBATCH --partition={args.partition}\n')
            f.write('#SBATCH --mail-type=FAIL\n')
            f.write(f'#SBATCH --mail-user={args.mail_user}\n\n')
            f.write('cd $SLURM_SUBMIT_DIR\n')
            f.write('/bin/hostname\n')
            f.write(f'echo "Job started: $(date)"\n\n')
            f.write(f'source {venv_activate}\n\n')
            f.write(f'{python_command}\n\n')
            f.write('echo "Job finished: $(date)"\n')

        submit_command = f'sbatch -p {args.partition} {job_file}'
        logging.info(f'Job file: {job_file}')
        logging.info(f'Submit command: {submit_command}')

        if args.test:
            logging.info('TEST MODE -- not submitting.')
            logging.info(f'Command that would run:\n{python_command}')
        else:
            subprocess.Popen(submit_command, shell=True).wait()

    else:
        ################################################################
        # Parallel mode: SLURM array job (one task per trial) + merge
        ################################################################
        # NOTE: this array/merge path targets the OLD generate_signal_events_CL.py,
        # which supported --trial-index and --merge-only.  The v2 generator
        # (generate_signal_events_CL_v2.py) does NOT expose those flags -- it
        # assigns a fresh run_id per invocation and merges via --merge-run-id.
        # Until the array/merge builders below are reworked for that model,
        # --parallel is incompatible with the v2 data-root tree; use the
        # sequential mode (default) instead.
        logging.error('--parallel is not supported with the v2 generator '
                      '(no --trial-index/--merge-only). Use sequential mode.')
        sys.exit(1)

        # Validate requirements
        if not args.save_to_file:
            logging.warning('--parallel requires --save-to-file; enabling it automatically.')
            args.save_to_file = True
        if not args.output_directory:
            logging.error('--parallel requires --output-directory to be set.')
            sys.exit(1)

        # --- Array job: each task runs a single trial ---
        cmd_parts = [
            f'python {script_dir}/generate_signal_events_CL.py',
            f'--masses-a {" ".join(str(m) for m in args.masses_a)}',
            f'--masses-dm {" ".join(str(m) for m in args.masses_dm)}',
            f'--nevents {args.nevents}',
            f'--ntrials 1',
            f'--trial-index $SLURM_ARRAY_TASK_ID',
            f'--seed $(( {args.base_seed} + $SLURM_ARRAY_TASK_ID ))',
            f'--dm-model {args.dm_model}',
            f'--detector-radius {args.detector_radius}',
            f'--detector-half-len {args.detector_half_len}',
            f'--inner-detector-radius {args.inner_detector_radius}',
            f'--inner-detector-half-len {args.inner_detector_half_len}',
        ] + _optional_physics_flags()

        python_command = ' \\\n    '.join(cmd_parts)

        # SLURM array spec
        array_spec = f'0-{args.ntrials - 1}'
        if args.max_concurrent is not None:
            array_spec += f'%{args.max_concurrent}'

        # Write the array job file
        job_file = f'{jobdir}/slurm_signal_array_{job_name}.job'
        with open(job_file, 'w') as f:
            f.write('#!/bin/bash\n\n')
            f.write(f'#SBATCH --job-name={job_name}\n')
            f.write(f'#SBATCH --array={array_spec}\n')
            f.write('#SBATCH --nodes=1 --ntasks-per-node=1\n')
            f.write('#SBATCH --ntasks=1\n')
            f.write('#SBATCH --cpus-per-task=1\n')
            f.write(f'#SBATCH --mem={args.mem}\n')
            f.write(f'#SBATCH --error={logdir}/slurm-%A_%a.err\n')
            f.write(f'#SBATCH --output={logdir}/slurm-%A_%a.out\n')
            f.write(f'#SBATCH --time={args.time}\n')
            f.write(f'#SBATCH --partition={args.partition}\n')
            f.write('#SBATCH --mail-type=FAIL\n')
            f.write(f'#SBATCH --mail-user={args.mail_user}\n\n')
            f.write('cd $SLURM_SUBMIT_DIR\n')
            f.write('/bin/hostname\n')
            f.write(f'echo "Array task $SLURM_ARRAY_TASK_ID started: $(date)"\n\n')
            f.write(f'source {venv_activate}\n\n')
            f.write(f'{python_command}\n\n')
            f.write(f'echo "Array task $SLURM_ARRAY_TASK_ID finished: $(date)"\n')

        # --- Merge job: runs after all array tasks complete ---
        merge_cmd_parts = [
            f'python {script_dir}/generate_signal_events_CL.py',
            '--merge-only',
            f'--masses-a {" ".join(str(m) for m in args.masses_a)}',
            f'--masses-dm {" ".join(str(m) for m in args.masses_dm)}',
            f'--dm-model {args.dm_model}',
        ] + _optional_physics_flags()

        merge_python_command = ' \\\n    '.join(merge_cmd_parts)

        merge_job_name = f'merge_{tag}_{current_time}'
        merge_job_file = f'{jobdir}/slurm_signal_merge_{job_name}.job'
        with open(merge_job_file, 'w') as f:
            f.write('#!/bin/bash\n\n')
            f.write(f'#SBATCH --job-name={merge_job_name}\n')
            f.write('#SBATCH --nodes=1 --ntasks-per-node=1\n')
            f.write('#SBATCH --ntasks=1\n')
            f.write('#SBATCH --cpus-per-task=1\n')
            f.write(f'#SBATCH --mem={args.mem}\n')
            f.write(f'#SBATCH --error={logdir}/slurm-merge-%j.err\n')
            f.write(f'#SBATCH --output={logdir}/slurm-merge-%j.out\n')
            f.write(f'#SBATCH --time=01:00:00\n')
            f.write(f'#SBATCH --partition={args.partition}\n')
            f.write('#SBATCH --mail-type=END,FAIL\n')
            f.write(f'#SBATCH --mail-user={args.mail_user}\n\n')
            f.write('cd $SLURM_SUBMIT_DIR\n')
            f.write('/bin/hostname\n')
            f.write(f'echo "Merge job started: $(date)"\n\n')
            f.write(f'source {venv_activate}\n\n')
            f.write(f'{merge_python_command}\n\n')
            f.write('echo "Merge job finished: $(date)"\n')

        # Submit array job, then dependent merge job
        submit_command = f'sbatch -p {args.partition} {job_file}'
        logging.info(f'Array job file: {job_file}')
        logging.info(f'Merge job file: {merge_job_file}')
        logging.info(f'Submit command: {submit_command}')

        if args.test:
            logging.info('TEST MODE -- not submitting.')
            logging.info(f'Array command per task:\n{python_command}')
            logging.info(f'Merge command:\n{merge_python_command}')
            logging.info(f'Merge would run with: --dependency=afterok:<ARRAY_JOB_ID>')
        else:
            result = subprocess.run(submit_command, shell=True, capture_output=True, text=True)
            print(result.stdout.strip())
            if result.returncode != 0:
                logging.error(f'Array job submission failed: {result.stderr}')
                sys.exit(1)

            # Parse job ID from "Submitted batch job 12345"
            array_job_id = result.stdout.strip().split()[-1]
            logging.info(f'Array job ID: {array_job_id}')

            # Submit merge job with dependency on all array tasks completing
            merge_submit = f'sbatch -p {args.partition} --dependency=afterok:{array_job_id} {merge_job_file}'
            logging.info(f'Merge submit command: {merge_submit}')
            merge_result = subprocess.run(merge_submit, shell=True, capture_output=True, text=True)
            print(merge_result.stdout.strip())
            if merge_result.returncode != 0:
                logging.error(f'Merge job submission failed: {merge_result.stderr}')


if __name__ == "__main__":
    main()
