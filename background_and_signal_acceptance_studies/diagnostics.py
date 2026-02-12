import pandas as pd
import numpy as np
import matplotlib.pylab as plt


def kinematic_diagnostic(d=-7.5, r=500, tag='mDM_200-10000_mA_0.22_dm_model_floating', masses=None, position_only=False, input_directory=None, DETECTOR_Y_CONSTRAINED=None):

    if input_directory is None:
        input_directory = './'

    tag = f'd_{d}_r_{r}_{tag}'
    if DETECTOR_Y_CONSTRAINED is not None:
        out_tag = f'd_{d}_r_{r}_{tag}_DETECTOR_Y_CONSTRAINED_{DETECTOR_Y_CONSTRAINED}'
    else:
        out_tag = f'd_{d}_r_{r}_{tag}'

        
    infile = f'{input_directory}/generated_data_{tag}.parquet'
    df_decays = pd.read_parquet(infile)

    print(infile)
    print(df_decays.shape)

    if masses is None:
        masses = df_decays['M_DM'].unique()

    df_decays['rho0_origin'] = np.sqrt(df_decays['x0']**2 + df_decays['z0']**2) 
    df_decays['phi0_origin'] = np.arctan2(df_decays['z0'],df_decays['x0']) 

    df_decays['p_pt_CMS'] = np.sqrt(df_decays['px1']**2 + df_decays['pz1']**2) 
    df_decays['p_phi_CMS'] = np.arctan(df_decays['pz1'],df_decays['px1']) 
    df_decays['p_eta_CMS'] = 0.5*np.log((np.sqrt(df_decays['pmag1'])+df_decays['pz1'])/(np.sqrt(df_decays['pmag1'])-df_decays['pz1'])) 


    for mass in masses:
        #mass = 1000
        
        filter = (df_decays['efinal_mu1']>1)
        filter = filter & (df_decays['M_DM']==mass)
        #filter = (df_decays['M_DM']==mass)

        # For comparison, constrain energy of muons to be
        # close to 1/2 the mass of the DM
        filter = filter & (np.abs(df_decays['e1']-mass/2)/mass < 0.01)

        if DETECTOR_Y_CONSTRAINED is not None:
            print(f"DETECTOR_Y_CONSTRAINED: {DETECTOR_Y_CONSTRAINED}")
            filter = filter & (df_decays['ip_y0']<DETECTOR_Y_CONSTRAINED)

        print(f'After selections: {len(df_decays[filter])}')

        ########################################################################
        plt.figure(figsize=(12,12))

        plt.subplot(3,3,1)
        df_decays[filter].plot.scatter(x='x0', y='y0', s=0.1, ax=plt.gca())
        plt.xlabel(r'Origin x (m)', fontsize=14)
        plt.ylabel(r'Origin y (m)', fontsize=14)

        plt.subplot(3,3,2)
        df_decays[filter].plot.scatter(x='z0', y='y0', s=0.1, ax=plt.gca())
        plt.xlabel(r'Origin z (m)', fontsize=14)
        plt.ylabel(r'Origin y (m)', fontsize=14)

        plt.subplot(3,3,3)
        df_decays[filter].plot.scatter(x='z0', y='x0', s=0.1, ax=plt.gca())
        plt.xlabel(r'Origin z (m)', fontsize=14)
        plt.ylabel(r'Origin x (m)', fontsize=14)
        plt.xlim(-900,900)
        plt.ylim(-900,900)

        #########################################################################
        max_xval = 150
        min_xval = -150
        nbins = 50
        xranges = (min_xval, max_xval)

        #plt.figure(figsize=(7,9))
        plt.subplot(3,3,4)
        df_decays[filter]['x0'].hist(bins=nbins, range=(-50,50), histtype="step", density=False, linewidth=2.5, label='Origin')
        df_decays[filter]['ip_x0'].hist(bins=nbins, range=(-50,50), histtype="step", density=False, linewidth=2.5, label='Detector entry')
        plt.xlabel(r'Origin x (m)', fontsize=14)
        plt.yscale('log')
        plt.legend()


        plt.subplot(3,3,5)
        df_decays[filter]['y0'].hist(bins=nbins, range=(-50,50), histtype="step", density=False, linewidth=2.5, label='Origin')
        df_decays[filter]['ip_y0'].hist(bins=nbins, range=(-50,50), histtype="step", density=False, linewidth=2.5, label='Detector entry')
        plt.xlabel(r'Origin y (m)', fontsize=14)
        plt.yscale('log')
        plt.legend()

        plt.subplot(3,3,6)
        df_decays[filter]['z0'].hist(bins=200, range=(-50,50), histtype="step", density=False, linewidth=2.5, label='Origin')
        df_decays[filter]['ip_z0'].hist(bins=200, range=(-50,50), histtype="step", density=False, linewidth=2.5, label='Detector entry')
        plt.xlabel(r'Origin z (m)', fontsize=14)
        plt.yscale('log')
        plt.legend()
        
        #plt.tight_layout()


        ########################################################################

        #plt.figure(figsize=(8,3))
        
        plt.subplot(3,3,7)
        df_decays[filter].plot.scatter(y='ip_y0', x='ip_x0', s=0.1, ax=plt.gca())
        plt.xlabel(r'Detector entry x (m)', fontsize=14)
        plt.ylabel(r'Detector entry y (m)', fontsize=14)
        plt.ylim(-2,2)
        plt.xlim(-2,2)

        plt.subplot(3,3,8)
        df_decays[filter].plot.scatter(y='ip_y0', x='ip_z0', s=0.1, ax=plt.gca())
        plt.xlabel(r'Detector entry z (m)', fontsize=14)
        plt.ylabel(r'Detector entry y (m)', fontsize=14)

        plt.subplot(3,3,9)
        df_decays[filter].plot.scatter(y='ip_x0', x='ip_z0', s=0.1, ax=plt.gca())
        plt.xlabel(r'Detector entry z (m)', fontsize=14)
        plt.ylabel(r'Detector entry x (m)', fontsize=14)

        plt.tight_layout()

        outfile = f'origin_and_detector_entry_d_{r}_r_{r}_{out_tag}.png'
        plt.savefig(outfile)
        ########################################################################

        ########################################################################
        layout = [
            ['A', 'A', 'B','B', 'C','C'],  # 
            ['D','D','D', 'E','E','E'],
            ['F','F','F', 'G','G','G']
        ]

        fig, axes = plt.subplot_mosaic(layout, figsize=(12, 12))

        #plt.subplot(3,3,1)
        plt.sca(axes['A'])
        df_decays[filter]['theta1'].hist(bins=64, range=(0,3.2), histtype="step", density=True, linewidth=2.5)
        plt.xlabel(r'$p_{\theta}$', fontsize=14)

        #plt.subplot(3,3,2)
        plt.sca(axes['B'])
        df_decays[filter]['eta1'].hist(bins=64, range=(0,3), density=True,  histtype="step",linewidth=2.5)
        plt.xlabel(r'$p_{\eta}$', fontsize=14)
        plt.yscale('log')

        #plt.subplot(3,3,3)
        plt.sca(axes['C'])
        df_decays[filter]['phi1'].hist(bins=64, range=(0,3.2), density=True,  histtype="step",linewidth=2.5)
        plt.xlabel(r'$p_{\phi}$', fontsize=14)
        plt.yscale('log')

        #plt.tight_layout()

        
        ########################################################################
        #plt.figure(figsize=(12,4))

        max_xval = 5000
        min_xval = 0
        nbins = 50
        xranges = (min_xval, max_xval)

        #plt.subplot(3,2,4)
        plt.sca(axes['D'])
        df_decays[filter]['e1'].hist(bins=nbins, range=xranges, histtype="step", density=True, linewidth=2.5, label='Orig. energy')
        df_decays[filter]['efinal_mu1'].hist(bins=nbins, range=xranges, histtype="step",density=True, linewidth=2.5,  label='Energy at detector')
        plt.xlabel(r'$E_{\mu}$ (GeV)', fontsize=14)
        plt.legend()

        plt.sca(axes['E'])
        df_decays[filter]['e1'].hist(bins=nbins, range=xranges, density=True,  histtype="step",linewidth=2.5, label='Orig. energy')
        df_decays[filter]['efinal_mu1'].hist(bins=nbins, range=xranges, density=True, histtype="step",linewidth=2.5,  label='Energy at detector')
        plt.xlabel(r'$E_{\mu}$ (GeV)', fontsize=14)
        plt.legend()
        plt.yscale('log')

        #plt.tight_layout()

        
        ########################################################################
        #plt.figure(figsize=(12,4))

        max_xval = 5000
        min_xval = 0
        nbins = 50
        xranges = (min_xval, max_xval)
        
        #plt.subplot(3,2,6)
        plt.sca(axes['F'])
        df_decays[filter]['pt1_detector_acceptance'].hist(bins=nbins, range=xranges,  density=True, histtype="step",linewidth=2.5, label='Ignoring eloss')
        df_decays[filter]['pt1_detector_acceptance_eloss'].hist(bins=nbins, range=xranges, density=True, histtype="step",linewidth=2.5,  label='With eloss')
        plt.xlabel(r'$p_{T}$ (GeV)', fontsize=14)
        plt.legend()

        #plt.subplot(3,2,7)
        plt.sca(axes['G'])
        df_decays[filter]['pt1_detector_acceptance'].hist(bins=nbins, range=xranges,  histtype="step", density=True, linewidth=2.5, label='Ignoring eloss')
        df_decays[filter]['pt1_detector_acceptance_eloss'].hist(bins=nbins, range=xranges, histtype="step",density=True, linewidth=2.5,  label='With eloss')
        plt.xlabel(r'$p_{T}$ (GeV)', fontsize=14)
        plt.ylim(1e-4)
        plt.legend()
        plt.yscale('log')
        
        plt.tight_layout()

        outfile = f'kinematic_d_{r}_r_{r}_{out_tag}.png'
        plt.savefig(outfile)

    
    return df_decays
##################################
