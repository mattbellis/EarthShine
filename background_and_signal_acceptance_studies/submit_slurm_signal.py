#!/usr/bin/env python
import argparse
import logging
import os
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
    parser.add_argument('--ntrials', type=int, default=3)
    parser.add_argument('--depth', type=float, nargs='+', default=[])
    parser.add_argument('--dm-model', type=str, default='floating')
    parser.add_argument('--save-to-file', action='store_true', default=False)
    parser.add_argument('--output-directory', type=str, default='OUTPUT_FILES')
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

    # Build the python command
    cmd_parts = [
        f'python {script_dir}/generate_signal_events_CL.py',
        f'--masses-a {" ".join(str(m) for m in args.masses_a)}',
        f'--masses-dm {" ".join(str(m) for m in args.masses_dm)}',
        f'--nevents {args.nevents}',
        f'--ntrials {args.ntrials}',
        f'--dm-model {args.dm_model}',
        f'--detector-radius {args.detector_radius}',
        f'--detector-half-len {args.detector_half_len}',
        f'--inner-detector-radius {args.inner_detector_radius}',
        f'--inner-detector-half-len {args.inner_detector_half_len}',
    ]
    if args.depth:
        cmd_parts.append(f'--depth {" ".join(str(d) for d in args.depth)}')
    if args.disk_radius is not None:
        cmd_parts.append(f'--disk-radius {args.disk_radius}')
    if args.save_to_file:
        cmd_parts.append('--save-to-file')
    if args.output_directory:
        cmd_parts.append(f'--output-directory {args.output_directory}')
    if args.additional_tag:
        cmd_parts.append(f'--additional-tag {args.additional_tag}')
    if args.average_eloss:
        cmd_parts.append('--average-eloss')
    if not args.save_only_detected:
        cmd_parts.append('--no-save-only-detected')
    if args.eloss_dict_name is not None:
        cmd_parts.append(f'--eloss-dict-name {args.eloss_dict_name}')
    if args.verbose:
        cmd_parts.append('--verbose')

    python_command = ' \\\n    '.join(cmd_parts)

    # Job name tag
    current_time = datetime.now().strftime("%Y%m%d_%H%M%S")
    tag = f'mA{"_".join(str(m) for m in args.masses_a)}_mDM{"_".join(str(m) for m in args.masses_dm)}'
    job_name = f'signal_{tag}_{current_time}'

    # Write the SLURM job file
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

    # Submit
    submit_command = f'sbatch -p {args.partition} {job_file}'
    logging.info(f'Job file: {job_file}')
    logging.info(f'Submit command: {submit_command}')

    if args.test:
        logging.info('TEST MODE -- not submitting.')
        logging.info(f'Command that would run:\n{python_command}')
    else:
        subprocess.Popen(submit_command, shell=True).wait()


if __name__ == "__main__":
    main()
