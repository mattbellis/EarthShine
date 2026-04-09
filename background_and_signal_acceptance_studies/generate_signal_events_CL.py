import argparse
#from dm_generation_tools import generate_many_events  # adjust this import

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

import dm_generation_tools as dgt
import pandas as pd
import glob
import logging



def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate dark matter scattering events."
    )

    parser.add_argument("--masses-a", type=float, nargs='+', default=[1],
                        help="List of A masses (default: [1])")
    parser.add_argument("--masses-dm", type=float, nargs='+', default=[100],
                        help="List of DM masses (default: [100])")
    parser.add_argument("--nevents", type=int, default=10,
                        help="Number of events to generate (default: 10)")
    parser.add_argument("--ntrials", type=int, default=3,
                        help="Number of files to generate before merging (default: 3)")
    parser.add_argument("--detector-radius", type=float, default=7.5,
                        help="Detector radius (default: 7.5)")
    parser.add_argument("--detector-half-len", type=float, default=15,
                        help="Detector half-length (default: 15)")
    parser.add_argument("--inner-detector-radius", type=float, default=1.0,
                        help="Inner detector radius (default: 1.0)")
    parser.add_argument("--inner-detector-half-len", type=float, default=2.5,
                        help="Inner detector half-length (default: 2.5)")
    parser.add_argument("--disk-radius", type=float, default=None,
                        help="Disk radius (default: None)")
    parser.add_argument("--depth", type=float, nargs='+', default=[],
                        help="Depth (default: None)")
    parser.add_argument("--dm-model", type=str, default="floating",
                        help="DM model name (default: 'floating')")
    parser.add_argument("--save-only-detected", action=argparse.BooleanOptionalAction,
                        default=True,
                        help="Save only detected events (default: True)")
    parser.add_argument("--eloss-dict-name", type=str, default=None,
                        help="Energy loss dictionary name (default: None)")
    parser.add_argument("--save-to-file", action="store_true", default=False,
                        help="Save output to file (default: False)")
    parser.add_argument("--additional-tag", type=str, default="",
                        help="Additional tag for output files (default: '')")
    parser.add_argument("--output-directory", type=str, default=None,
                        help="Output directory (default: None)")
    parser.add_argument("--average-eloss", action="store_true", default=False,
                        help="Use average energy loss (default: False)")
    parser.add_argument("--verbose", action=argparse.BooleanOptionalAction,
                        default=False,
                        help="Print out some debugging information (default: False)")

    return parser.parse_args()


################################################################################
def main():
    
    logging.getLogger("tensorflow").setLevel(logging.FATAL)
    ############################################################################

    args = parse_args()

    # If only one depth is passed in, then use that
    if len(args.depth)==1:
        args.depth = args.depth[0]

    for i in range(0,args.ntrials):
        print(f'\nTrial {i} out of {args.ntrials} --------------------------------------------------------')
        print(f'Generate {args.nevents} events...')

        additional_tag=f'_{i:06d}',
        if args.additional_tag is not None:
            additional_tag=f'_{args.additional_tag}_{i:06d}'

        df_decays, tag = dgt.generate_many_events(
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
            verbose=args.verbose
        )

        ##############################################################
        # Merge the "trial" files
        ##############################################################
        partial_tag = tag[0:-6]
        if args.verbose:
            print(f'partial_tag: {partial_tag}')

        files = glob.glob(f'{args.output_directory}/generated_data_{partial_tag}[0-9]*.parquet')

        #print('files')
        #print(files)

    # Merge them
    dfs = []
    for i,filename in enumerate(files):
      #print(filename)
      dfs.append( pd.read_parquet(filename))

    print(f'\n Merging {len(dfs)} intermediate files...\n')
    df_decays = pd.concat(dfs)
    df_decays.to_parquet(f'{args.output_directory}/generated_data_{partial_tag}COMBINED.parquet')

    # Delete the intermediate dataframes to free up memory
    del df_decays
    del dfs[:]


if __name__ == "__main__":
    main()
