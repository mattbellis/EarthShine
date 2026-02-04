import pandas as pd
import numpy as np
import matplotlib.pylab as plt


def kinematic_diagnostic(d=-7.5, r=500, tag='mDM_200-10000_mA_0.22_dm_model_floating', masses=None, position_only=False):

    tag = f'd_{d}_r_{r}_{tag}'
        
    infile = f'generated_data_{tag}.parquet'
    df_decays = pd.read_parquet(infile)

    if masses is None:
        masses = df_decays['M_DM'].unique()

    df_decays['rho0_origin'] = np.sqrt(df_decays['x0']**2 + df_decays['z0']**2) 
    df_decays['phi0_origin'] = np.arctan(df_decays['z0'],df_decays['x0']) 
    df_decays['theta1'] = np.arccos(df_decays['costh1']) 

    df_decays['p_pt_CMS'] = np.sqrt(df_decays['px1']**2 + df_decays['pz1']**2) 
    df_decays['p_phi_CMS'] = np.arctan(df_decays['pz1'],df_decays['px1']) 
    df_decays['p_eta_CMS'] = 0.5*np.log((np.sqrt(df_decays['pmag1'])+df_decays['pz1'])/(np.sqrt(df_decays['pmag1'])-df_decays['pz1'])) 

    #df_decays['theta1'] = np.arccos(df_decays['costh1']) 

    for mass in masses:
        #mass = 1000
        
        #filter = (df_decays['efinal_mu1']>1)
        #filter = filter & (df_decays['M_DM']==mass)
        filter = (df_decays['M_DM']==mass)

        # For comparison
        filter = filter & (np.abs(df_decays['e1']-mass/2)/mass < 0.01)

        print(f'After selections: {len(df_decays[filter])}')

        ########################################################################
        plt.figure(figsize=(8,6))

        plt.subplot(2,3,1)
        df_decays[filter].plot.scatter(x='x0', y='y0', s=0.1, ax=plt.gca())
        plt.xlabel(r'Origin x (m)', fontsize=18)
        plt.ylabel(r'Origin y (m)', fontsize=18)

        plt.subplot(2,3,2)
        df_decays[filter].plot.scatter(x='z0', y='y0', s=0.1, ax=plt.gca())
        plt.xlabel(r'Origin z (m)', fontsize=18)
        plt.ylabel(r'Origin y (m)', fontsize=18)

        plt.subplot(2,3,3)
        df_decays[filter].plot.scatter(x='z0', y='x0', s=0.1, ax=plt.gca())
        plt.xlabel(r'Origin z (m)', fontsize=18)
        plt.ylabel(r'Origin x (m)', fontsize=18)

        plt.tight_layout()

        #########################################################################
        max_xval = 150
        min_xval = -150
        nbins = 50
        xranges = (min_xval, max_xval)

        plt.figure(figsize=(7,9))
        plt.subplot(3,1,1)
        df_decays[filter]['x0'].hist(bins=nbins, range=(-30,20), histtype="step", density=False, linewidth=2.5, label='Origin')
        df_decays[filter]['ip_x0'].hist(bins=nbins, range=(-30,20), histtype="step", density=False, linewidth=2.5, label='Detector IP')
        plt.xlabel(r'Origin x (m)', fontsize=18)
        plt.yscale('log')
        plt.legend()


        plt.subplot(3,1,2)
        df_decays[filter]['y0'].hist(bins=nbins, range=(-30,20), histtype="step", density=False, linewidth=2.5, label='Origin')
        df_decays[filter]['ip_y0'].hist(bins=nbins, range=(-30,20), histtype="step", density=False, linewidth=2.5, label='Detector IP')
        plt.xlabel(r'Origin y (m)', fontsize=18)
        plt.yscale('log')
        plt.legend()

        plt.subplot(3,1,3)
        df_decays[filter]['z0'].hist(bins=200, range=(-30,20), histtype="step", density=False, linewidth=5.5, label='Origin')
        df_decays[filter]['ip_z0'].hist(bins=200, range=(-30,20), histtype="step", density=False, linewidth=2.5, label='Detector IP')
        plt.xlabel(r'Origin z (m)', fontsize=18)
        #plt.yscale('log')
        plt.legend()
        
        plt.tight_layout()


        ########################################################################

        plt.figure(figsize=(8,6))
        
        plt.subplot(2,3,4)
        df_decays[filter].plot.scatter(y='ip_y0', x='ip_x0', s=0.1, ax=plt.gca())
        plt.xlabel(r'Detector IP x (m)', fontsize=18)
        plt.ylabel(r'Detector IP y (m)', fontsize=18)

        plt.subplot(2,3,5)
        df_decays[filter].plot.scatter(y='ip_y0', x='ip_z0', s=0.1, ax=plt.gca())
        plt.xlabel(r'Detector IP z (m)', fontsize=18)
        plt.ylabel(r'Detector IP y (m)', fontsize=18)

        plt.subplot(2,3,6)
        df_decays[filter].plot.scatter(y='ip_x0', x='ip_z0', s=0.1, ax=plt.gca())
        plt.xlabel(r'Detector IP z (m)', fontsize=18)
        plt.ylabel(r'Detector IP x (m)', fontsize=18)

        plt.tight_layout()

        ########################################################################

        plt.figure(figsize=(12,8))
        
        plt.subplot(2,2,1)
        plt.hist(df_decays[filter]['z0'],bins=100, range=(-4000,0))
        plt.xlabel('depth (m)', fontsize=18)

        plt.subplot(2,2,2)
        plt.hist(df_decays[filter]['rho0_origin'],bins=100, range=(-10,500))
        plt.xlabel('radial distance (m)', fontsize=18)
        
        plt.subplot(2,2,3)
        df_decays[filter].plot.scatter(y='y0', x='x0', s=0.1, alpha=0.1, ax=plt.gca())        
        
        plt.subplot(2,2,4)
        #df_decays[filter].plot.scatter(y='efinal_mu1', x='y0', s=0.1, ax=plt.gca())
        df_decays[filter].plot.scatter(y='pt1_detector_acceptance_eloss', x='rho0_origin', s=0.1, ax=plt.gca())
        
        DMstr = 'DM'
        lbracket = '{'
        rbracket = '}'
        plt.gcf().suptitle(f'$M_{lbracket}DM{rbracket}$ {int(mass)} GeV/c$^2$')
        
        plt.tight_layout()
        
        #outfile = 'depth_and_pt_d_{r}_r_{r}_{tag}.png'
        #plt.savefig(outfile)

        
        ########################################################################
        plt.figure(figsize=(12,4))

        max_xval = 5000
        min_xval = 0
        nbins = 50
        xranges = (min_xval, max_xval)

        plt.subplot(1,3,1)
        df_decays[filter]['e1'].hist(bins=nbins, range=xranges, histtype="step", density=True, linewidth=2.5, label='Orig. energy')
        df_decays[filter]['efinal_mu1'].hist(bins=nbins, range=xranges, histtype="step",density=True, linewidth=2.5,  label='Energy at detector')
        plt.xlabel(r'$E_{\mu}$ (GeV)', fontsize=18)
        plt.legend()

        plt.subplot(1,3,2)
        df_decays[filter]['e1'].hist(bins=nbins, range=xranges, density=True,  histtype="step",linewidth=2.5, label='Orig. energy')
        df_decays[filter]['efinal_mu1'].hist(bins=nbins, range=xranges, density=True, histtype="step",linewidth=2.5,  label='Energy at detector')
        plt.xlabel(r'$E_{\mu}$ (GeV)', fontsize=18)
        plt.legend()
        plt.yscale('log')
        #plt.ylim(5e-4)

        plt.subplot(1,3,3)
        df_decays[filter].plot.scatter(y='efinal_mu1', x='rho0_origin', s=0.1, ax=plt.gca())
        plt.xlabel(r'Radial distance (m)', fontsize=18)
        plt.ylabel(r'$E_{\mu}$ (GeV)', fontsize=18)

        plt.tight_layout()

        
        ########################################################################
        plt.figure(figsize=(12,4))

        max_xval = 5000
        min_xval = 0
        nbins = 50
        xranges = (min_xval, max_xval)
        
        plt.subplot(1,3,1)
        df_decays[filter]['pt1_detector_acceptance'].hist(bins=nbins, range=xranges,  histtype="step",linewidth=2.5, label='Ignoring eloss')
        df_decays[filter]['pt1_detector_acceptance_eloss'].hist(bins=nbins, range=xranges, histtype="step",linewidth=2.5,  label='With eloss')
        plt.xlabel(r'$p_{T}$ (GeV)', fontsize=18)
        plt.legend()

        plt.subplot(1,3,2)
        df_decays[filter]['pt1_detector_acceptance'].hist(bins=nbins, range=xranges,  histtype="step", density=True, linewidth=2.5, label='Ignoring eloss')
        df_decays[filter]['pt1_detector_acceptance_eloss'].hist(bins=nbins, range=xranges, histtype="step",density=True, linewidth=2.5,  label='With eloss')
        plt.xlabel(r'$p_{T}$ (GeV)', fontsize=18)
        plt.legend()
        plt.yscale('log')
        
        plt.subplot(1,3,3)
        df_decays[filter].plot.scatter(y='pt1_detector_acceptance', x='rho0_origin', s=0.1, ax=plt.gca())
        plt.xlabel(r'Radial distance (m)', fontsize=18)
        plt.ylabel(r'$p_{T}$ (GeV) (no eloss)', fontsize=18)

        plt.tight_layout()

        outfile = f'depth_and_pt_d_{r}_r_{r}_{tag}.png'
        plt.savefig(outfile)

        ########################################################################
        plt.figure(figsize=(12,4))

        plt.subplot(1,3,1)
        df_decays[filter]['phi1'].hist(bins=50, range=(-1.6,1.6),  histtype="step",linewidth=2.5, label='Muons')
        plt.xlabel(r'$p_{\phi}$ (GeV)', fontsize=18)
        plt.legend()

        plt.subplot(1,3,2)
        df_decays[filter]['pt1_detector_acceptance'].hist(bins=50, range=(-100,3500),  histtype="step", density=True, linewidth=2.5, label='Ignoring eloss')
        df_decays[filter]['pt1_detector_acceptance_eloss'].hist(bins=50, range=(-100,3500), histtype="step",density=True, linewidth=2.5,  label='With eloss')
        plt.xlabel(r'$p_{T}$ (GeV)', fontsize=18)
        plt.legend()
        plt.yscale('log')
        
        plt.subplot(1,3,3)
        df_decays[filter].plot.scatter(y='pt1_detector_acceptance', x='phi0_origin', s=0.1, ax=plt.gca())
        plt.xlabel(r'Radial distance (m)', fontsize=18)
        plt.ylabel(r'$p_{T}$ (GeV) (no eloss)', fontsize=18)

        plt.tight_layout()

        outfile = f'SECOND_depth_and_pt_d_{r}_r_{r}_{tag}.png'
        plt.savefig(outfile)

        ######################################################################################################
        plt.figure(figsize=(8,6))

        plt.subplot(2,2,1)
        df_decays[filter]['phi1'].hist(bins=50, range=(-1.6,1.6),  histtype="step",linewidth=2.5, label='Muons')
        plt.xlabel(r'$p_{\phi}$ (GeV)', fontsize=18)
        plt.legend()

        plt.subplot(2,2,2)
        df_decays[filter]['phi1'].hist(bins=50, range=(-1.6,1.6),  histtype="step",linewidth=2.5, label='Muons')
        plt.xlabel(r'$p_{\phi}$ (GeV)', fontsize=18)
        plt.yscale('log')
        plt.legend()


        plt.subplot(2,2,3)
        df_decays[filter]['phi0_origin'].hist(bins=50, range=(-1.6,1.6),  histtype="step",linewidth=2.5, label='Muons')
        plt.xlabel(r'Origin $\phi$', fontsize=18)
        plt.legend()

        plt.subplot(2,2,4)
        df_decays[filter]['phi0_origin'].hist(bins=50, range=(-1.6,1.6),  histtype="step",linewidth=2.5, label='Muons')
        plt.xlabel(r'Origin $\phi$', fontsize=18)
        plt.yscale('log')
        plt.legend()

        plt.tight_layout()

        ######################################################################################################
        plt.figure(figsize=(8,6))

        plt.subplot(2,2,1)
        df_decays[filter]['costh1'].hist(bins=100, range=(-0.1,1.1),  histtype="step",linewidth=2.5, label='Muons')
        plt.xlabel(r'$p_{\cos \theta}$ (GeV)', fontsize=18)
        plt.legend()

        plt.subplot(2,2,2)
        df_decays[filter]['theta1'].hist(bins=50, range=(0,2),  histtype="step",linewidth=2.5, label='Muons')
        plt.xlabel(r'$p_{\theta}$ (GeV)', fontsize=18)
        plt.legend()

        plt.subplot(2,2,3)
        df_decays[filter]['costh1'].hist(bins=100, range=(-.1,1.1),  histtype="step",linewidth=2.5, label='Muons')
        plt.xlabel(r'$p_{\cos \theta}$ (GeV)', fontsize=18)
        plt.yscale('log')
        plt.legend()

        plt.subplot(2,2,4)
        df_decays[filter]['theta1'].hist(bins=50, range=(0.,2),  histtype="step",linewidth=2.5, label='Muons')
        plt.xlabel(r'$p_{\theta}$ (GeV)', fontsize=18)
        plt.yscale('log')
        plt.legend()


        plt.tight_layout()

    
    return df_decays
##################################
