import pandas as pd
import numpy as np
import matplotlib.pylab as plt


##########################################################################
def kinematic_diagnostic(depth=-8, diskR=500, tag='mDM_200-10000_mA_0.22_dm_model_floating', masses=None, position_only=False, \
        INNER_DETECTOR_FILTER=False, input_directory=None, DETECTOR_Y_CONSTRAINED=None, CONSTRAIN_ENERGY=False, \
        plotdir='plots'):

    if input_directory is None:
        input_directory = './'

    tag = f'depth_{depth}_diskR_{diskR}_{tag}'

    if DETECTOR_Y_CONSTRAINED is not None:
        out_tag = f'{tag}_DETECTOR_Y_CONSTRAINED_{DETECTOR_Y_CONSTRAINED}'
    else:
        out_tag = tag

    if INNER_DETECTOR_FILTER is True:
        out_tag = f'{out_tag}_INN_TK_CONST'
        
        
    infile = f'{input_directory}/generated_data_{tag}.parquet'
    df_decays = pd.read_parquet(infile)

    print(infile)
    print(df_decays.shape)

    print('masses in file: ')
    print(df_decays['M_DM'].unique())
    if masses is None:
        masses = df_decays['M_DM'].unique()

    df_decays['rho0_origin'] = np.sqrt(df_decays['x0']**2 + df_decays['y0']**2)
    df_decays['phi0_origin'] = np.arctan2(df_decays['y0'],df_decays['x0'])

    df_decays['p_pt_CMS'] = np.sqrt(df_decays['px1']**2 + df_decays['py1']**2)
    df_decays['p_phi_CMS'] = np.arctan2(df_decays['py1'],df_decays['px1']) 
    df_decays['p_eta_CMS'] = 0.5*np.log((np.sqrt(df_decays['pmag1'])+df_decays['pz1'])/(np.sqrt(df_decays['pmag1'])-df_decays['pz1'])) 


    for mass in masses:
        #mass = 1000
        
        filter = (df_decays['efinal_mu1']>1)
        filter = filter & (df_decays['M_DM']==mass)
        #filter = (df_decays['M_DM']==mass)

        # For comparison purposes, constrain energy of muons to be
        # close to 1/2 the mass of the DM
        if CONSTRAIN_ENERGY:
            print(f"CONSTRAIN_ENERGY: {CONSTRAIN_ENERGY}")
            filter = filter & (np.abs(df_decays['e1']-mass/2)/mass < 0.01)

        if DETECTOR_Y_CONSTRAINED is not None:
            print(f"DETECTOR_Y_CONSTRAINED: {DETECTOR_Y_CONSTRAINED}")
            filter = filter & (df_decays['ip_y0']<DETECTOR_Y_CONSTRAINED)

        if INNER_DETECTOR_FILTER:
            print(f"Requiring that inner detector is hit")
            filter = filter & (df_decays['hit_inner_detector'])

        print(f'After selections: {len(df_decays[filter])}')

        ########################################################################
        plt.figure(figsize=(12,12))

        x = df_decays[filter]['x0']
        y = df_decays[filter]['y0']
        z = df_decays[filter]['z0']

        plt.subplot(3,3,1)
        #df_decays[filter].plot.scatter(x='x0', y='y0', s=0.1, ax=plt.gca())
        plt.hist2d(x=x, y=y, bins=100, range=([-1*float(diskR), float(diskR)], [-2500,0]))
        plt.xlabel(r'Origin x (m)', fontsize=14)
        plt.ylabel(r'Origin y (m)', fontsize=14)

        plt.subplot(3,3,2)
        #df_decays[filter].plot.scatter(x='z0', y='y0', s=0.1, ax=plt.gca())
        plt.hist2d(x=z, y=y, bins=100, range=([-1*float(diskR), float(diskR)], [-2500,0]))
        plt.xlabel(r'Origin z (m)', fontsize=14)
        plt.ylabel(r'Origin y (m)', fontsize=14)

        plt.subplot(3,3,3)
        #df_decays[filter].plot.scatter(x='z0', y='x0', s=0.1, ax=plt.gca())
        plt.hist2d(x=z, y=x, bins=100, range=([-1*float(diskR), float(diskR)], [-1*float(diskR),float(diskR)]))
        plt.xlabel(r'Origin z (m)', fontsize=14)
        plt.ylabel(r'Origin x (m)', fontsize=14)
        #plt.xlim(-900,900)
        #plt.ylim(-900,900)

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
        
        '''
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
        '''
        mask = (df_decays[filter]['hit_inner_detector'])

        #'''
        plt.subplot(3,3,7)
        plt.gca().set_aspect('equal')
        df_decays[filter].plot.scatter(x='ip_x0', y='ip_y0', s=1,  ax=plt.gca(), label='muon system')
        df_decays[filter].plot.scatter(x='ip_inner_x0', y='ip_inner_y0', color='r', marker='s',s=1, ax=plt.gca(), label='inner tk')
        df_decays[filter][mask].plot.scatter(x='ip_x0', y='ip_y0', color='k', s=20,  ax=plt.gca(), label='muon sys. when inner tk')
        plt.legend()
        plt.xlim(-9, 9)
        plt.ylim(-9, 9)
        plt.gca().set_xlabel('x (meters)', fontsize=14)
        plt.gca().set_ylabel('y (meters)', fontsize=14)


        plt.subplot(3,3,8)
        plt.gca().set_aspect('equal')
        df_decays[filter].plot.scatter(y='ip_x0', x='ip_z0', s=1,  ax=plt.gca(), label='muon system')
        df_decays[filter].plot.scatter(y='ip_inner_x0', x='ip_inner_z0', color='r', marker='s',s=1, ax=plt.gca(), label='inner tk')
        df_decays[filter][mask].plot.scatter(y='ip_x0', x='ip_z0', color='k', s=20,  ax=plt.gca(), label='muon sys. when inner tk')
        plt.legend(loc='upper left')
        plt.xlim(-20, 20)
        plt.ylim(-20, 20)
        plt.gca().set_ylabel('x (meters)', fontsize=14)
        plt.gca().set_xlabel('z (meters)', fontsize=14)

        plt.subplot(3,3,9)
        plt.gca().set_aspect('equal')
        plt.gca().hist2d(df_decays[filter][mask]['ip_z0'], df_decays[filter][mask]['ip_x0'], bins=100, range=([-20, 20], [-20, 20]), norm='log')
        plt.xlim(-20, 20)
        plt.ylim(-20, 20)
        plt.xlabel('z (meters)', fontsize=14)
        plt.ylabel('x (meters)', fontsize=14)
        #'''
    

        plt.tight_layout()

        outfile = f'{plotdir}/origin_{depth}_diskR_{diskR}_{out_tag}.png'
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

        #max_xval = 5000
        max_xval = 1.1*mass/2
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

        #max_xval = 5000
        max_xval = 1.1*mass/2
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
        #plt.ylim(1e-5, 1e-2)
        plt.legend()
        plt.yscale('log')
        
        plt.tight_layout()

        outfile = f'{plotdir}/kinematic_depth_{depth}_diskR_{diskR}_{out_tag}.png'
        plt.savefig(outfile)

    
    return df_decays
##################################

##########################################################################
def detector_hits_plots(depth=-8, diskR=500, tag='mDM_200-10000_mA_0.22_dm_model_floating', masses=None, position_only=False, \
        INNER_DETECTOR_FILTER=False, input_directory=None, DETECTOR_Y_CONSTRAINED=None, CONSTRAIN_ENERGY=False, \
        plotdir='plots', xrange=20, yrange=20, zrange=None):

    if input_directory is None:
        input_directory = './'

    tag = f'depth_{depth}_diskR_{diskR}_{tag}'

    if DETECTOR_Y_CONSTRAINED is not None:
        out_tag = f'{tag}_DETYCON_{DETECTOR_Y_CONSTRAINED}'
    else:
        out_tag = tag

    if INNER_DETECTOR_FILTER is True:
        out_tag = f'{out_tag}_INNCON'

    out_tag = out_tag.replace('_HIT_DETECTOR','').replace('_COMBINED','')
    out_tag = out_tag.replace('momentum_constrained','mom_con')
        
    infile = f'{input_directory}/generated_data_{tag}.parquet'
    df_decays = pd.read_parquet(infile)

    print(infile)
    print(df_decays.shape)

    print('masses in file: ')
    print(df_decays['M_DM'].unique())
    if masses is None:
        masses = df_decays['M_DM'].unique()

    df_decays['rho0_origin'] = np.sqrt(df_decays['x0']**2 + df_decays['y0']**2)
    df_decays['phi0_origin'] = np.arctan2(df_decays['y0'],df_decays['x0'])

    df_decays['p_pt_CMS'] = np.sqrt(df_decays['px1']**2 + df_decays['py1']**2)
    df_decays['p_phi_CMS'] = np.arctan2(df_decays['py1'],df_decays['px1']) 
    df_decays['p_eta_CMS'] = 0.5*np.log((np.sqrt(df_decays['pmag1'])+df_decays['pz1'])/(np.sqrt(df_decays['pmag1'])-df_decays['pz1'])) 


    dict_hist = {}

    for mass in masses:
        #mass = 1000
        
        filter = (df_decays['efinal_mu1']>1)
        filter = filter & (df_decays['M_DM']==mass)
        #filter = (df_decays['M_DM']==mass)

        # For comparison purposes, constrain energy of muons to be
        # close to 1/2 the mass of the DM
        if CONSTRAIN_ENERGY:
            print(f"CONSTRAIN_ENERGY: {CONSTRAIN_ENERGY}")
            filter = filter & (np.abs(df_decays['e1']-mass/2)/mass < 0.01)

        if DETECTOR_Y_CONSTRAINED is not None:
            print(f"DETECTOR_Y_CONSTRAINED: {DETECTOR_Y_CONSTRAINED}")
            filter = filter & (df_decays['ip_y0']<DETECTOR_Y_CONSTRAINED)

        if INNER_DETECTOR_FILTER:
            print(f"Requiring that inner detector is hit")
            filter = filter & (df_decays['hit_inner_detector'])

        print(f'After selections: {len(df_decays[filter])}')

        ########################################################################
        plt.figure(figsize=(5,4))

        x = df_decays[filter]['x0']
        y = df_decays[filter]['y0']
        z = df_decays[filter]['z0']

        ip_x0 = df_decays[filter]['ip_x0']
        ip_y0 = df_decays[filter]['ip_y0']
        ip_z0 = df_decays[filter]['ip_z0']
        if xrange is None:
            xrange=float(diskR)
        if yrange is None:
            yrange=1.1*(abs(float(depth)))
        if zrange is None:
            zrange=float(diskR)

        #print('Ranges...')
        #print(xrange, yrange, zrange)
        #print(type(xrange), type(yrange), type(zrange))

        norm = None
        #norm = 'log'

        plt.subplot(1,1,1)
        #df_decays[filter].plot.scatter(x='x0', y='y0', s=0.1, ax=plt.gca())
        plt.hist2d(x=ip_z0, y=ip_x0, bins=50, range=([-1*xrange, xrange], [-zrange,zrange]), norm=norm)
        plt.xlabel(r'Origin x (m)', fontsize=14)
        plt.ylabel(r'Origin y (m)', fontsize=14)
        plt.colorbar(label='Counts')


        dmstr = '{DM}'
        title = f'$M_{dmstr}$={int(mass)}  depth={depth}  diskR={diskR}'
        print(title)
        plt.suptitle(title, fontsize=14)

        plt.tight_layout()

        outfile = f'{plotdir}/detector_hits_mass_{mass}_{out_tag}.png'
        plt.savefig(outfile)


    df_hist = pd.DataFrame.from_dict(dict_hist)
    outfile = f'HISTOGRAM_FILES_FOR_RATES/detector_hits_masses_{masses[0]}_{masses[-1]}_{out_tag}.parquet'
    print(f'Saving file to {outfile}')
    df_hist.to_parquet(outfile)

    return df_decays

######################################################################################################
def origin_plots(depth=-8, diskR=500, tag='mDM_200-10000_mA_0.22_dm_model_floating', masses=None, position_only=False, \
        INNER_DETECTOR_FILTER=False, input_directory=None, DETECTOR_Y_CONSTRAINED=None, CONSTRAIN_ENERGY=False, \
        plotdir='plots', xrange=None, yrange=None, zrange=None):

    if input_directory is None:
        input_directory = './'

    tag = f'depth_{depth}_diskR_{diskR}_{tag}'

    if DETECTOR_Y_CONSTRAINED is not None:
        out_tag = f'{tag}_DETYCON_{DETECTOR_Y_CONSTRAINED}'
    else:
        out_tag = tag

    if INNER_DETECTOR_FILTER is True:
        out_tag = f'{out_tag}_INNCON'

    out_tag = out_tag.replace('_HIT_DETECTOR','').replace('_COMBINED','')
    out_tag = out_tag.replace('momentum_constrained','mom_con')
        
    infile = f'{input_directory}/generated_data_{tag}.parquet'
    df_decays = pd.read_parquet(infile)

    print(infile)
    print(df_decays.shape)

    print('masses in file: ')
    print(df_decays['M_DM'].unique())
    if masses is None:
        masses = df_decays['M_DM'].unique()

    df_decays['rho0_origin'] = np.sqrt(df_decays['x0']**2 + df_decays['y0']**2)
    df_decays['phi0_origin'] = np.arctan2(df_decays['y0'],df_decays['x0'])

    df_decays['p_pt_CMS'] = np.sqrt(df_decays['px1']**2 + df_decays['py1']**2)
    df_decays['p_phi_CMS'] = np.arctan2(df_decays['py1'],df_decays['px1']) 
    df_decays['p_eta_CMS'] = 0.5*np.log((np.sqrt(df_decays['pmag1'])+df_decays['pz1'])/(np.sqrt(df_decays['pmag1'])-df_decays['pz1'])) 


    dict_hist = {}

    for mass in masses:
        #mass = 1000
        
        filter = (df_decays['efinal_mu1']>1)
        filter = filter & (df_decays['M_DM']==mass)
        #filter = (df_decays['M_DM']==mass)

        # For comparison purposes, constrain energy of muons to be
        # close to 1/2 the mass of the DM
        if CONSTRAIN_ENERGY:
            print(f"CONSTRAIN_ENERGY: {CONSTRAIN_ENERGY}")
            filter = filter & (np.abs(df_decays['e1']-mass/2)/mass < 0.01)

        if DETECTOR_Y_CONSTRAINED is not None:
            print(f"DETECTOR_Y_CONSTRAINED: {DETECTOR_Y_CONSTRAINED}")
            filter = filter & (df_decays['ip_y0']<DETECTOR_Y_CONSTRAINED)

        if INNER_DETECTOR_FILTER:
            print(f"Requiring that inner detector is hit")
            filter = filter & (df_decays['hit_inner_detector'])

        print(f'After selections: {len(df_decays[filter])}')

        ########################################################################
        plt.figure(figsize=(12,4))

        x = df_decays[filter]['x0']
        y = df_decays[filter]['y0']
        z = df_decays[filter]['z0']

        if xrange is None:
            xrange=float(diskR)
        if yrange is None:
            yrange=1.1*(abs(float(depth)))
        if zrange is None:
            zrange=float(diskR)

        #print('Ranges...')
        #print(xrange, yrange, zrange)
        #print(type(xrange), type(yrange), type(zrange))

        norm = None
        #norm = 'log'

        plt.subplot(1,3,1)
        #df_decays[filter].plot.scatter(x='x0', y='y0', s=0.1, ax=plt.gca())
        plt.hist2d(x=x, y=y, bins=50, range=([-1*xrange, xrange], [-yrange,0]), norm=norm)
        plt.xlabel(r'Origin x (m)', fontsize=14)
        plt.ylabel(r'Origin y (m)', fontsize=14)
        plt.colorbar(label='Counts')


        plt.subplot(1,3,2)
        #df_decays[filter].plot.scatter(x='z0', y='y0', s=0.1, ax=plt.gca())
        plt.hist2d(x=z, y=y, bins=50, range=([-1*zrange, zrange], [-yrange,0]), norm=norm)
        plt.xlabel(r'Origin z (m)', fontsize=14)
        plt.ylabel(r'Origin y (m)', fontsize=14)
        plt.colorbar(label='Counts')

        plt.subplot(1,3,3)
        #df_decays[filter].plot.scatter(x='z0', y='x0', s=0.1, ax=plt.gca())
        plt.hist2d(x=z, y=x, bins=50, range=([-1*zrange, zrange], [-1*xrange,xrange]), norm=norm)
        plt.xlabel(r'Origin z (m)', fontsize=14)
        plt.ylabel(r'Origin x (m)', fontsize=14)
        plt.colorbar(label='Counts')

        dmstr = '{DM}'
        title = f'$M_{dmstr}$={int(mass)}  depth={depth}  diskR={diskR}'
        print(title)
        plt.suptitle(title, fontsize=14)

        plt.tight_layout()

        outfile = f'{plotdir}/origin_mass_{mass}_{out_tag}.png'
        plt.savefig(outfile)

        # ALSO DO PROJECTION
        plt.figure(figsize=(12,4))

        plt.subplot(1,2,1)
        plt.hist(y,bins=100, range=(-yrange,0), density=True)
        plt.xlabel('Depth (m)', fontsize=18)
        plt.tight_layout()

        plt.subplot(1,2,2)
        #plt.hist(y,bins=100)
        bins = np.logspace(np.log10(1), np.log10(4000), 34)
        counts, bin_edges, _ = plt.hist(-y,bins=bins, density=True)
        plt.xscale('log')

        dict_hist[f'{mass:d} bin_start'] =  bin_edges[:-1]
        dict_hist[f'{mass:d} count'] = counts 

        plt.xlabel('Depth (m) [ZOOMED IN]', fontsize=18)
        plt.tight_layout()

        outfile = f'{plotdir}/origin_depth_1D_mass_{mass}_{out_tag}.png'
        plt.savefig(outfile)

    df_hist = pd.DataFrame.from_dict(dict_hist)
    outfile = f'HISTOGRAM_FILES_FOR_RATES/dist_from_depths_at_detector_masses_{masses[0]}_{masses[-1]}_{out_tag}.parquet'
    print(f'Saving file to {outfile}')
    df_hist.to_parquet(outfile)

    return df_decays

######################################################################################################
##########################################################################
def e_and_pt_plots(depth=-8, diskR=500, tag='mDM_200-10000_mA_0.22_dm_model_floating', masses=None, position_only=False, \
        INNER_DETECTOR_FILTER=False, input_directory=None, DETECTOR_Y_CONSTRAINED=None, CONSTRAIN_ENERGY=False, \
        plotdir='plots', max_xval=None):

    if input_directory is None:
        input_directory = './'

    tag = f'depth_{depth}_diskR_{diskR}_{tag}'

    if DETECTOR_Y_CONSTRAINED is not None:
        out_tag = f'{tag}_DETYCON_{DETECTOR_Y_CONSTRAINED}'
    else:
        out_tag = tag

    if INNER_DETECTOR_FILTER is True:
        out_tag = f'{out_tag}_INN_TK_CONST'
        
    out_tag = out_tag.replace('_HIT_DETECTOR','').replace('_COMBINED','')
    out_tag = out_tag.replace('momentum_constrained','mom_con')
        
    infile = f'{input_directory}/generated_data_{tag}.parquet'
    df_decays = pd.read_parquet(infile)

    print(infile)
    print(df_decays.shape)

    print('masses in file: ')
    print(df_decays['M_DM'].unique())
    if masses is None:
        masses = df_decays['M_DM'].unique()

    df_decays['rho0_origin'] = np.sqrt(df_decays['x0']**2 + df_decays['y0']**2)
    df_decays['phi0_origin'] = np.arctan2(df_decays['y0'],df_decays['x0'])

    df_decays['p_pt_CMS'] = np.sqrt(df_decays['px1']**2 + df_decays['py1']**2)
    df_decays['p_phi_CMS'] = np.arctan2(df_decays['py1'],df_decays['px1']) 
    df_decays['p_eta_CMS'] = 0.5*np.log((np.sqrt(df_decays['pmag1'])+df_decays['pz1'])/(np.sqrt(df_decays['pmag1'])-df_decays['pz1'])) 


    for mass in masses:
        #mass = 1000
        print(f'mass: {mass}')
        
        filter = (df_decays['efinal_mu1']>1)
        filter = filter & (df_decays['M_DM']==mass)
        #filter = (df_decays['M_DM']==mass)

        # For comparison purposes, constrain energy of muons to be
        # close to 1/2 the mass of the DM
        if CONSTRAIN_ENERGY:
            print(f"CONSTRAIN_ENERGY: {CONSTRAIN_ENERGY}")
            filter = filter & (np.abs(df_decays['e1']-mass/2)/mass < 0.01)

        if DETECTOR_Y_CONSTRAINED is not None:
            print(f"DETECTOR_Y_CONSTRAINED: {DETECTOR_Y_CONSTRAINED}")
            filter = filter & (df_decays['ip_y0']<DETECTOR_Y_CONSTRAINED)

        if INNER_DETECTOR_FILTER:
            print(f"Requiring that inner detector is hit")
            filter = filter & (df_decays['hit_inner_detector'])

        print(f'After selections: {len(df_decays[filter])}')

        ########################################################################
        plt.figure(figsize=(12,4))

        min_xval = 0
        nbins = 50
        if max_xval is None:
            xranges = (min_xval, 1.5*mass/2)

        print('Ranges...')
        print(max_xval)
        print(type(max_xval))

        dmstr = '{DM}'

        plt.subplot(1,2,1)
        df_decays[filter]['e1'].hist(bins=nbins, range=xranges, density=True,  histtype="step",linewidth=2.5, label='Orig. energy')
        df_decays[filter]['efinal_mu1'].hist(bins=nbins, range=xranges, density=True, histtype="step",linewidth=2.5,  label='Energy at detector')
        plt.xlabel(r'$E_{\mu}$ (GeV)', fontsize=14)
        plt.legend()
        plt.yscale('log')

        plt.subplot(1,2,2)
        df_decays[filter]['pt1_detector_acceptance'].hist(bins=nbins, range=xranges,  histtype="step", density=True, linewidth=2.5, label='Ignoring eloss')
        df_decays[filter]['pt1_detector_acceptance_eloss'].hist(bins=nbins, range=xranges, histtype="step",density=True, linewidth=2.5,  label='With eloss')
        plt.xlabel(r'$p_{T}$ (GeV)', fontsize=14)
        plt.legend()
        plt.yscale('log')

        text = f'$M_{dmstr}$={int(mass)}\ndepth={depth}\ndiskR={diskR}'
        plt.gca().text(0.7, 0.1, text,transform=plt.gca().transAxes)

        title = f'$M_{dmstr}$={int(mass)}  depth={depth}  diskR={diskR}'
        print(title)
        plt.suptitle(title, fontsize=14)

        plt.tight_layout()

        outfile = f'{plotdir}/e_and_pt_mass_{mass}_{out_tag}.png'
        plt.savefig(outfile)

    return df_decays

######################################################################################################
##########################################################################
def compare_pt_plots(depth=-8, diskR=500, infiles=[], labels=[], masses=None, position_only=False, \
        INNER_DETECTOR_FILTER=False, input_directory=None, DETECTOR_Y_CONSTRAINED=None, CONSTRAIN_ENERGY=False, \
        plotdir='plots', max_xval=None):

    if input_directory is None:
        input_directory = './'

    #tag = f'depth_{depth}_diskR_{diskR}_{tag}'

    #if DETECTOR_Y_CONSTRAINED is not None:
        #out_tag = f'{tag}_DETYCON_{DETECTOR_Y_CONSTRAINED}'
    #else:
        ##out_tag = tag

    #if INNER_DETECTOR_FILTER is True:
        #out_tag = f'{out_tag}_INN_TK_CONST'
        
    #out_tag = out_tag.replace('_HIT_DETECTOR','').replace('_COMBINED','')
    #out_tag = out_tag.replace('momentum_constrained','mom_con')
        
    infile1 = f'{input_directory}/{infiles[0]}'
    infile2 = f'{input_directory}/{infiles[1]}'
    df_decays1 = pd.read_parquet(infile1)
    df_decays2 = pd.read_parquet(infile2)

    print(infile1)
    print(df_decays1.shape)

    print('masses in file: ')
    print(df_decays1['M_DM'].unique())
    if masses is None:
        masses = df_decays['M_DM'].unique()

    plt.figure(figsize=(6,4))

    for idx, df_decays in enumerate([df_decays1, df_decays2]):
        df_decays['rho0_origin'] = np.sqrt(df_decays['x0']**2 + df_decays['y0']**2)
        df_decays['phi0_origin'] = np.arctan2(df_decays['y0'],df_decays['x0'])

        df_decays['p_pt_CMS'] = np.sqrt(df_decays['px1']**2 + df_decays['py1']**2)
        df_decays['p_phi_CMS'] = np.arctan2(df_decays['py1'],df_decays['px1']) 
        df_decays['p_eta_CMS'] = 0.5*np.log((np.sqrt(df_decays['pmag1'])+df_decays['pz1'])/(np.sqrt(df_decays['pmag1'])-df_decays['pz1'])) 


        for mass in masses:
            #mass = 1000
            print(f'mass: {mass}')
            
            filter = (df_decays['efinal_mu1']>1)
            filter = filter & (df_decays['M_DM']==mass)
            #filter = (df_decays['M_DM']==mass)

            # For comparison purposes, constrain energy of muons to be
            # close to 1/2 the mass of the DM
            if CONSTRAIN_ENERGY:
                print(f"CONSTRAIN_ENERGY: {CONSTRAIN_ENERGY}")
                filter = filter & (np.abs(df_decays['e1']-mass/2)/mass < 0.01)

            if DETECTOR_Y_CONSTRAINED is not None:
                print(f"DETECTOR_Y_CONSTRAINED: {DETECTOR_Y_CONSTRAINED}")
                filter = filter & (df_decays['ip_y0']<DETECTOR_Y_CONSTRAINED)

            if INNER_DETECTOR_FILTER:
                print(f"Requiring that inner detector is hit")
                filter = filter & (df_decays['hit_inner_detector'])

            print(f'After selections: {len(df_decays[filter])}')

            ########################################################################

            min_xval = 0
            nbins = 50
            if max_xval is None:
                xranges = (min_xval, 1.5*mass/2)

            print('Ranges...')
            print(max_xval)
            print(type(max_xval))

            dmstr = '{DM}'

            plt.subplot(1,1,1)
            #df_decays[filter]['pt1_detector_acceptance'].hist(bins=nbins, range=xranges,  histtype="step", density=True, linewidth=2.5, label='Ignoring eloss')
            df_decays[filter]['pt1_detector_acceptance_eloss'].hist(bins=nbins, range=xranges, histtype="step",density=True, linewidth=2.5,  label=labels[idx])

    plt.xlabel(r'$p_{T}$ (GeV)', fontsize=14)
    plt.legend()
    plt.yscale('log')

    text = f'$M_{dmstr}$={int(mass)}'#\ndepth={depth}\ndiskR={diskR}'
    plt.gca().text(0.7, 0.1, text,transform=plt.gca().transAxes)

    title = f'$M_{dmstr}$={int(mass)}'#  depth={depth}  diskR={diskR}'
    print(title)
    plt.suptitle(title, fontsize=14)

    plt.tight_layout()

    #outfile = f'{plotdir}/compare_pt_mass_{mass}_{out_tag}.png'
    outfile = f'{plotdir}/compare_pt_mass_{mass}.png'
    plt.savefig(outfile)

    return df_decays1,df_decays2


######################################################################################################
