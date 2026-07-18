import matplotlib
matplotlib.use('Agg')
import diagnostics as dg
import warnings

# Suppress all warnings
warnings.filterwarnings("ignore")


# Read in a file and make some diagnostic plots
DETECTOR_Y_CONSTRAINED = -0.4
# This tells us what file to read in
'''
generated_data_d_-17.5_r_1003_mDM_1000-7000_mA_0.22_dm_model_momentum_constrained_HIT_DETECTOR_COMBINED.parquet
generated_data_d_-17.5_r_40_mDM_1000-7000_mA_0.22_dm_model_momentum_constrained_HIT_DETECTOR_COMBINED.parquet

generated_data_d_-1000_r_40_mDM_1000-7000_mA_0.22_dm_model_momentum_constrained_HIT_DETECTOR_COMBINED.parquet
generated_data_d_-1000_r_6000_mDM_1000-7000_mA_0.22_dm_model_momentum_constrained_HIT_DETECTOR_COMBINED.parquet
'''
# base_directory = '/home/danyi/earthAsDM/EarthShine/background_and_signal_acceptance_studies'
base_directory = './'


input_directory = f'{base_directory}/OUTPUT_FILES'
# Available masses in OUTPUT_FILES: 100.0, 1000.0, 3000.0, 7000.0, 10000.0, 100000.0
masses = [100.0, 1000.0, 3000.0, 7000.0, 10000.0, 100000.0]
depth = -8.0
diskR = 459

# Process each mass file individually since they are separate files
for mass in masses:
    tag = f'mDM_{mass}_mA_0.22_dmModel_momentum_constrained_HIT_DETECTOR_ave_eloss_TESTING_combined'
    print(f"\nProcessing mass: {mass} GeV")
    df_decays = dg.kinematic_diagnostic(masses=[mass], depth=depth, diskR=diskR, input_directory=input_directory, tag=tag,
                                        DETECTOR_Y_CONSTRAINED=DETECTOR_Y_CONSTRAINED)
    # The function returns the dataframe
    df_decays.info()
