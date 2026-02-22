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
tag = 'mDM_7000_mA_0.22_dmModel_momentum_constrained_HIT_DETECTOR_COMBINED'


'''
depths = [-10]
disk_radii = [40, 100, 500, 1000]
'''

df_decays = dg.kinematic_diagnostic(masses=[7000], depth=-18, diskR=40, input_directory='OUTPUT_FILES', tag=tag,
                                    DETECTOR_Y_CONSTRAINED=DETECTOR_Y_CONSTRAINED)
df_decays = dg.kinematic_diagnostic(masses=[7000], depth=-18, diskR=100, input_directory='OUTPUT_FILES', tag=tag,
                                    DETECTOR_Y_CONSTRAINED=DETECTOR_Y_CONSTRAINED)
df_decays = dg.kinematic_diagnostic(masses=[7000], depth=-18, diskR=500, input_directory='OUTPUT_FILES', tag=tag,
                                    DETECTOR_Y_CONSTRAINED=DETECTOR_Y_CONSTRAINED)
df_decays = dg.kinematic_diagnostic(masses=[7000], depth=-18, diskR=1000, input_directory='OUTPUT_FILES', tag=tag,
                                    DETECTOR_Y_CONSTRAINED=DETECTOR_Y_CONSTRAINED)

# The function returns the dataframe
df_decays.info()
