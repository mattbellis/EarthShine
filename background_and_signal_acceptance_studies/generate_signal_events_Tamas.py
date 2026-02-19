import matplotlib
matplotlib.use('Agg')

import os
import sys

# Redirect OS-level fd 2 to suppress C++ noise (CUDA, XLA, absl)
# but keep a Python stderr writer so exceptions still show
os.environ['CUDA_VISIBLE_DEVICES'] = '-1'
_real_stderr = os.dup(2)
os.dup2(os.open(os.devnull, os.O_WRONLY), 2)
sys.stderr = os.fdopen(_real_stderr, "w", closefd=False)

import dm_generation_tools as dgt
import pandas as pd
import glob
import logging

logging.getLogger("tensorflow").setLevel(logging.FATAL)


# print(f'{np.__version__ = }')
# print(f'{phasespace.__version__ = }')
# print(f'{tensorflow.__version__ = }')

######################################################################
# Masses in GeV/c^2
MASSES_A = [.220]
# MASSES_DM = [1000, 7000]
MASSES_DM = [7000]

######################################################################
# Set the depths
######################################################################
# Depths are in meters. 0 depth is *below* the bottom of the detector
#
# Values should be negative
#
# When the depth values are a list, then it generates the muon origins uniformly
# over a volume

# Plane
# depths = [-10, -1000]
depths = [-10]

# Volume
#depths = [[-100,-200]]

# When the radius is None, it calculates the radius based on an angle of 91 degrees
#disk_radii = [0, 40, 100, 500, 1000, None]
disk_radii = [1000]

# Dark matter model
# "floating" - unform direction for dark photon direction, but constrained to go "up"
#dm_model = 'floating'

# "core" - dark photons are coming straight up (all momentum in y)
#dm_model = 'core'

# "momentum_constrained" unphysical model with dark photon at rest, but allows
# us to fix the magnitude of the momentum of the muon
# Uniform in direction
dm_model = 'momentum_constrained'

# Output directory for intermediate files and summary file
output_directory = 'OUTPUT_FILES'

# Do you want to save all the events or not?
# If you set this to be true, then try it with a small number of events and only 1 "trial" file
#SAVE_ONLY_DETECTED = False
SAVE_ONLY_DETECTED = True

##################################################
# Size of the detector
#
# Note that we have a text string (tracker_tag) that
# gets added to file names.
##################################################
# Hits the tracker
#detector_radius = 1.5
#half_len = 1
#tracker_tag = 'TRACKER_VOL_'

# Full Detector
detector_radius = 8 #m
half_len = 10
tracker_tag = ''

######################################################################
# Generate!
######################################################################

for disk_radius in disk_radii:
  for depth in depths:
    # Depth is from below the bottom of the detector, so we add 8 meters to get the depth from the CMS origin
    depth = depth - 8 # meters
    print(f'Disk radius: {disk_radius}, depth from CMS origin: {depth}')
    save_to_file = True

    #nevents = 1_000_000
    # For testing
    # nevents = 10
    nevents = 1000000

    # Trials is not the right word, but these are the number of intermediate files that get produced
    # and then combined into one file.
    #
    # This is useful when you only want to save events that hit the detector
    #ntrials = 1000
    # ntrials = 1
    ntrials = 50

    # You might want to have more trials if you are very deep and have a low acceptance
    '''
    if depth < -200 and depth > -999:
      ntrials = 2000
    elif depth <= -999:
      ntrials = 50000
    '''

    ##############################################################
    # Generate all the files!!!!!!
    for i in range(0,ntrials):
      print(f'Trial {i} --------------------------------------------------------')
      print(f'Generate {nevents} events...')
      df_decays, tag = dgt.generate_many_events(
          MASSES_A = MASSES_A,
          MASSES_DM = MASSES_DM,
          nevents_to_generate = nevents,
          disk_radius = disk_radius,
          depth = depth,
          dm_model = dm_model,
          detector_radius = detector_radius,
          half_len = half_len,
          SAVE_ONLY_DETECTED = SAVE_ONLY_DETECTED,
          eloss_dict_name = 'eloss_dictionary_11032025_v2.pkl',
          save_to_file = save_to_file,
          additional_tag = f'_{tracker_tag}{i:06d}',
          output_directory = output_directory
        )
    ##############################################################


    ##############################################################
    # Merge the "trial" files
    ##############################################################
    partial_tag = tag[0:-6]
    print(f'partial_tag: {partial_tag}')

    files = glob.glob(f'{output_directory}/generated_data_{partial_tag}[0-9]*.parquet')

    #print('files')
    #print(files)

    # Merge them
    dfs = []
    for i,filename in enumerate(files):
      #print(filename)
      dfs.append( pd.read_parquet(filename))

    print(f'\n Merging {len(dfs)} intermediate files...\n')
    df_decays = pd.concat(dfs)
    df_decays.to_parquet(f'{output_directory}/generated_data_{partial_tag}COMBINED.parquet')

    # Delete the intermediate dataframes to free up memory
    del df_decays
    del dfs[:]
