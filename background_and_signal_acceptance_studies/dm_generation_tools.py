import numpy as np
import pandas as pd
import detector_simulation_tools as dst
import eloss_tools

import eloss_average

import phasespace

import time
import pickle
import os

################################################################################
def make_directory(directory_name = None):

    if directory_name is None:
        print("No directory name passed in!\n")
        return -1


    # Create the directory
    try:
        os.mkdir(directory_name)
        print(f"Directory '{directory_name}' created successfully.")
    except FileExistsError:
        print(f"Directory '{directory_name}' already exists.")
    except PermissionError:
        print(f"Permission denied: Unable to create '{directory_name}'.")
    except Exception as e:
        print(f"An error occurred: {e}")

    return 0

################################################################################

def opening_angle(p4s):
  p0mag = np.sqrt(p4s[0][0]**2 + p4s[0][1]**2 + p4s[0][2]**2)
  p1mag = np.sqrt(p4s[1][0]**2 + p4s[1][1]**2 + p4s[1][2]**2)

  dot_product = p4s[0][0]*p4s[1][0] + p4s[0][1]*p4s[1][1] + p4s[0][2]*p4s[1][2]

  theta = np.arccos(dot_product/(p0mag*p1mag))

  return theta


def distance(v1, v2):
  dx = v1[0] - v2[0]
  dy = v1[1] - v2[1]
  dz = v1[2] - v2[2]

  d = np.sqrt(dx**2 + dy**2 + dz**2)
  return d

################################################################################

################################################################################
import numpy as np

def energy_loss(E0, particle, material, distance_cm=None):
    """
    Vectorized version. E0 can be a float or a numpy array of energies in GeV.

    Parameters:
    - E0: float or np.ndarray, initial energy in GeV
    - particle: 'muon', 'pion', or 'kaon'
    - material: one of 'rock', 'copper', 'lead', 'water', 'aluminum', 'iron', 'carbon', 'ice'
    - distance_cm: optional float. If None, return range in cm; else return final energy in GeV

    Returns:
    - If distance_cm is None: range(s) in cm
    - If distance_cm is provided: final energy/energies in GeV
    """
    particle = particle.lower()
    material = material.lower()
    E0 = np.asarray(E0)  # support scalar or array

    material_db = {
        "rock":     {"rho": 2.65,  "a_mu": 2.0e-3, "b_mu": 4.0e-6,  "lambda_pi": 120, "lambda_K": 160},
        "copper":   {"rho": 8.96,  "a_mu": 1.75e-3, "b_mu": 1.4e-5, "lambda_pi": 106, "lambda_K": 130},
        "lead":     {"rho": 11.34, "a_mu": 1.5e-3, "b_mu": 5.4e-5, "lambda_pi": 106, "lambda_K": 130},
        "water":    {"rho": 1.0,   "a_mu": 2.0e-3, "b_mu": 3.0e-6,  "lambda_pi": 120, "lambda_K": 160},
        "aluminum": {"rho": 2.7,   "a_mu": 1.9e-3, "b_mu": 6.5e-6,  "lambda_pi": 118, "lambda_K": 155},
        "iron":     {"rho": 7.87,  "a_mu": 1.75e-3, "b_mu": 1.7e-5, "lambda_pi": 109, "lambda_K": 135},
        "carbon":   {"rho": 2.267, "a_mu": 1.9e-3, "b_mu": 0.9e-6,  "lambda_pi": 120, "lambda_K": 160},
        "ice":      {"rho": 0.917, "a_mu": 2.0e-3, "b_mu": 3.0e-6,  "lambda_pi": 120, "lambda_K": 160},
    }

    if material not in material_db:
        raise ValueError(f"Unsupported material '{material}'.")

    mat = material_db[material]
    rho = mat["rho"]

    if particle == "muon":
        a_cm = mat["a_mu"] * rho
        b_cm = mat["b_mu"] * rho

        if distance_cm is None:
            if b_cm > 0:
                return (1 / b_cm) * np.log(1 + (b_cm * E0 / a_cm)), rho
            else:
                return E0 / a_cm, rho
        else:
            if b_cm > 0:
                epsilon = a_cm / b_cm
                return np.maximum((E0 + epsilon) * np.exp(-b_cm * distance_cm) - epsilon, 0), rho
            else:
                return np.maximum(E0 - a_cm * distance_cm, 0), rho

    elif particle in ["pion", "kaon"]:
        a_ion = 2.0e-3
        dEdx_ion = a_ion * rho
        f = 0.6
        lambda_gcm2 = mat["lambda_pi"] if particle == "pion" else mat["lambda_K"]
        lambda_cm = lambda_gcm2 / rho

        def simulate_range_or_final_energy(E_start):
            E = E_start
            x = 0
            dist_since_int = 0
            step = 1  # cm
            if distance_cm is None:
                while E > 0:
                    E -= dEdx_ion * step
                    dist_since_int += step
                    x += step
                    if dist_since_int >= lambda_cm:
                        E *= f
                        dist_since_int = 0
                #print("HERE", x, rho, type(x), type(rho))
                return x, rho
            else:
                while x < distance_cm and E > 0:
                    E -= dEdx_ion * step
                    dist_since_int += step
                    x += step
                    if dist_since_int >= lambda_cm:
                        E *= f
                        dist_since_int = 0
                #print('THERE', max(E,0), rho)
                return max(E, 0), rho

        # Vectorize the above simulation for numpy arrays
        return np.vectorize(simulate_range_or_final_energy)(E0)#, rho

    else:
        raise ValueError("Unsupported particle. Choose 'muon', 'pion', or 'kaon'.")

################################################################################
##################################################################

def generate_dm_decays(MASSES_A=[.250,1,5], MASSES_DM=[10,100,1000], nevents_to_generate=10, dm_model='core', 
                       verbose=False):
    # Loop over different values

    decays = {}
    decays['M_DM'] = []
    decays['M_A'] = []

    decays['px1'] = []
    decays['py1'] = []
    decays['pz1'] = []
    decays['e1'] = []

    decays['px2'] = []
    decays['py2'] = []
    decays['pz2'] = []
    decays['e2'] = []

    #MASSES_A = [250,1000,5000]
    #MUON_MASS = 105.11

    #MASSES_A = [.250,1,5]
    #MASSES_DM = [10,100,1000]

    MUON_MASS = 0.10511

    thetas = []
    pmags = []

    muon_p = []

    #nevents_to_generate = 1000

    #pmags_GeV = [10,100,1000]

    for MASSES_A in MASSES_A:
      for pmag in MASSES_DM:

        #pmag = pmag_GeV*1000

        #print(f'M_A: {MASSES_A}   pmag: {pmag}')

        # Here is the boost vectors to boost it to the moving frame of the dark photon, pmag
        boost_vector, boost_vectors = None,None

        if dm_model=='core':
            # The core model has all the photons with their momentum in the z-direction
            # Old
            #boost_vector = np.array([0,0, pmag, np.sqrt(pmag**2 + MASSES_A**2)])
            # Shifted CMS defintion so that y is up and z is along beamline
            boost_vector = np.array([0, pmag, 0, np.sqrt(pmag**2 + MASSES_A**2)])
            boost_vectors = np.tile(boost_vector, (nevents_to_generate,1))

        elif dm_model=='floating':
            # In the floating model, the decay of the dark photon is uniform
            # Generate in spherical coordinates
            # Generate theta from 0 to 1. This should constrain them to pointing up
            #theta = np.arccos(np.random.random(nevents_to_generate))
            #phi = np.pi*np.random.random(nevents_to_generate)# + np.pi # From pi to 2*pi
            #r = pmag*np.ones(nevents_to_generate)
            
            # Claude
            cos_alpha = np.random.random(nevents_to_generate)#rng.uniform(0, 1, size=n)  # cos(alpha) in [0, 1] -> alpha in [0, pi/2]
            alpha = np.arccos(cos_alpha)           # angle from +y axis
            #beta = rng.uniform(0, 2 * np.pi, size=n)  # rotation around y-axis
            beta = 2*np.pi*np.random.random(nevents_to_generate)

            # Unit direction vector
            px_unit = np.sin(alpha) * np.cos(beta)
            py_unit = cos_alpha  # = cos(alpha), always positive
            pz_unit = np.sin(alpha) * np.sin(beta)

            # Scale by momentum magnitude
            px = pmag * px_unit
            py = pmag * py_unit
            pz = pmag * pz_unit


            E = np.sqrt(pmag**2 + MASSES_A**2)*np.ones(nevents_to_generate)
            
            # Convert to Cartesian for phasespace
            #px = r*np.sin(theta)*np.cos(phi)
            #py = r*np.sin(theta)*np.sin(phi)
            #pz = r*np.cos(theta)
            #print(r)
            #print(E)
            #px, py, pz, p, pt, theta, eta = dst.cms_to_cartesian(p=r, phi=phi, theta=theta)

            boost_vectors = np.array([px, py, pz, E]).T
            # Switching this so that y is up
            #boost_vectors = np.array([px, py, py, E]).T
            #print(boost_vectors)

            #boost_vector = np.array([0,0, pmag, np.sqrt(pmag**2 + MASSES_A**2)])
            #boost_vectors = np.tile(boost_vector, (nevents_to_generate,1))

        # Generate a bunch of vectors of mass 0 assuming there is no boost from the photon. 
        # Instead, we get muons with fixed momentum of 1/2 DM_MASS
        elif dm_model=='momentum_constrained':

            px = np.zeros(nevents_to_generate)
            py = np.zeros(nevents_to_generate)
            pz = np.zeros(nevents_to_generate)
            E = pmag*np.ones(nevents_to_generate)
            boost_vectors = np.array([px, py, pz, E]).T

        else:
            print(f"Dark matter model '{dm_model}' has not been implemented!")
            print("No events will be generated")
            return None

        #if verbose:
        print(f'M_A: {MASSES_A:5.3f}   pmag: {pmag:8.1f}     generating {nevents_to_generate:8d} decays')

        weights, particles = phasespace.nbody_decay(MASSES_A,
                                                    [MUON_MASS, MUON_MASS]).generate(n_events=nevents_to_generate, boost_to=boost_vectors)

        if verbose:
            print(f"Generated the decays")

        if verbose:
            print(f"Pulling out the values from the decays and filling our dataframe...")

        p0 = particles['p_0'][:].numpy().T
        p1 = particles['p_1'][:].numpy().T

        #data[MASSES_A][pmag] = [p0.T, p1.T]
        decays['M_DM'] += (pmag*np.ones(nevents_to_generate)).tolist()
        decays['M_A'] += (MASSES_A*np.ones(nevents_to_generate)).tolist()

        #print(len(decays['M_A']))

        decays['px1'] += p0[0].tolist()
        decays['py1'] += p0[1].tolist()
        decays['pz1'] += p0[2].tolist()
        decays['e1'] +=  p0[3].tolist()

        decays['px2'] += p1[0].tolist()
        decays['py2'] += p1[1].tolist()
        decays['pz2'] += p1[2].tolist()
        decays['e2'] +=  p1[3].tolist()
        if verbose:
            print(f"Pulled out the values from the decays and filled our dataframe...")

    dfdec1 = pd.DataFrame.from_dict(decays)

    # Generate a few other entries
    if verbose:
        print("Generated all the events and now calculating a few new values!")


    # pmag, theta (degrees), phi
    px1 = dfdec1['px1'].values
    py1 = dfdec1['py1'].values
    pz1 = dfdec1['pz1'].values
    e1 = dfdec1['e1'].values

    # This should still work for CMS coordinates
    #pt1 = np.sqrt(px1*px1 + py1*py1)
    #pmag1 = mag([px1,py1,pz1])
    #theta1 = np.atan2(pt1/pz1)
    #phi1 = np.arctan2(py1,px1)
    p, pt, eta, phi, theta = dst.cartesian_to_cms(px1, py1, pz1)

    if verbose:
        print("Calculated pmag1, theta1, phi1")

    dfdec1['pmag1'] = p
    dfdec1['pt1'] = pt
    dfdec1['eta1'] = eta
    dfdec1['theta1'] = theta
    dfdec1['phi1'] = phi

    #####################################

    px2 = dfdec1['px2'].values
    py2 = dfdec1['py2'].values
    pz2 = dfdec1['pz2'].values
    e2 = dfdec1['e2'].values

    p, pt, eta, phi, theta = dst.cartesian_to_cms(px2, py2, pz2)

    #pt2 = np.sqrt(px1*px1 + py1*py1)
    #pmag2 = mag([px2,py2,pz2])
    #theta2 = np.atan2(pt2/pz2)
    #phi2 = np.arctan2(py2,px2)

    if verbose:
        print("Calculated pmag2, theta2, phi2")

    #dfdec1['pmag2'] = pmag2
    #dfdec1['theta2'] = theta2
    #dfdec1['phi2'] = phi2

    dfdec1['pmag2'] = p
    dfdec1['pt2'] = pt
    dfdec1['eta2'] = eta
    dfdec1['theta2'] = theta
    dfdec1['phi2'] = phi

    dfdec1['costh1'] = np.cos(dfdec1['theta1'])
    dfdec1['costh2'] = np.cos(dfdec1['theta2'])

    # Opening angle
    p4s = [[px1,py1,pz1,e1], [px2,py2,pz2,e2]]
    thetas = np.rad2deg(opening_angle(p4s))
    dfdec1['opening angle'] = thetas
    if verbose:
        print("Calculated opening values")

    #dfdec1.sample(5)

    return dfdec1



################################################################################

def generate_many_events(MASSES_A=[1], MASSES_DM=[100], nevents_to_generate=10, \
                         detector_radius=7.5, detector_half_len=10.5, 
                         inner_detector_radius=1.0, inner_detector_half_len=2.5, 
                         disk_radius=None,  depth=-7.5, 
                         #radius=20, depth=-7.5, 
                         dm_model='floating', 
                         SAVE_ONLY_DETECTED=True, eloss_dict_name=None, 
                             save_to_file=False, additional_tag="", output_directory=None, 
                         AVERAGE_ELOSS=False, verbose=False):

    #######################################################################
    # Create the output directory
    #######################################################################
    if output_directory is not None:
        retval = make_directory(output_directory)
        if retval == -1:
            print("Output directory not created. Stopping DM generation.\n")
            return -1
    else:
        output_directory = './'
    
    #######################################################################
    # Generate the decays of the dark photons
    #######################################################################
    
    start_time = time.time()
    #df_decays = generate_dm_decays(MASSES_A=MASSES_A, MASSES_DM=MASSES_DM, nevents_to_generate=nevents)
    df_decays = generate_dm_decays(MASSES_A=MASSES_A, MASSES_DM=MASSES_DM, nevents_to_generate=nevents_to_generate, dm_model=dm_model)
    
    print(f"\nTime to generate events {time.time() - start_time:.2} seconds  ------")

    #######################################################################
    # Place the decay points
    #######################################################################

    nentries = len(df_decays)

    #origins, directions = dst.generate_origins_and_directions(nevents=nentries, radius=radius, depth=depth, dm_model=dm_model)
    origins, disk_radius = dst.generate_origins(nevents=nentries, disk_radius=disk_radius, depth=depth)

    #print(origins)
    
    #######################################################################
    # What strikes the detector
    #######################################################################

    # These directions have magnitude 1
    directions1, directions2 = dst.directions_from_momentum(df_decays)

    # We multiply the directions so that they are incredibly long vectors to make sure the line
    # is long enough to potentially intersect the detector
    # JUST DOING THIS FOR ONE MUON FOR NOW
    # NEED TO FIX THIS LATER
    #detector_entry_pts,detector_exit_pts = dst.intersect_finite_cylinder_x_np(origins=origins, directions=1e6*directions1)
    # For tracker
    #pts0,pts1 = dst.intersect_finite_cylinder_x_np(origins=origins, directions=1e6*directions1, \
    #                                               radius=detector_radius, half_len=detector_half_len)

    # From Claude with CMS coordinates
    detector_entry_pts,detector_exit_pts = dst.ray_cylinder_intersection(origins=origins, directions=directions1, \
                                                   radius=detector_radius, half_length=detector_half_len)

    # Did it hit the inner detector
    detector_entry_pts_inner,detector_exit_pts_inner = dst.ray_cylinder_intersection(origins=origins, directions=directions1, \
                                                   radius=inner_detector_radius, half_length=inner_detector_half_len)
    
    hit_detector_idx = ~np.isnan(detector_entry_pts.T[0])
    hit_inner_detector_idx = ~np.isnan(detector_entry_pts_inner.T[0])
#=======
    ##detector_entry_pts,detector_exit_pts = dst.intersect_finite_cylinder_x_np(origins=origins, directions=1e6*directions1, \
            ##                                               radius=detector_radius, half_len=half_len)
#
    ## From Claude with CMS coordinates
    #detector_entry_pts,detector_exit_pts = dst.ray_cylinder_intersection(origins=origins, directions=directions1, \
            #radius=detector_radius, half_length=half_len)
    ##print(detector_entry_pts)
#>>>>>>> 6e6c78dba70c28c55320ef121e77eb13ef9c0833
    
    hit_detector_idx = ~np.isnan(detector_entry_pts.T[0])
    n_muons_hit = np.sum(hit_detector_idx)

    if verbose:
        print(f"Number of muons that struck the detector: {n_muons_hit} / {len(hit_detector_idx)} ({100*n_muons_hit/len(hit_detector_idx):.2f}%)")

    df_decays['hit_detector'] = hit_detector_idx
    df_decays['hit_inner_detector'] = hit_inner_detector_idx
    
    df_decays['x0'] = origins.T[0]
    df_decays['y0'] = origins.T[1]
    df_decays['z0'] = origins.T[2]
    
    #distance_to_detector = np.sqrt(df_decays['x0']*df_decays['x0'] + df_decays['y0']*df_decays['y0'] + df_decays['z0']*df_decays['z0'])
    #df_decays['distance_to_detector'] = distance_to_detector
    
    # Save the intersection points
    df_decays['ip_x0'] = detector_entry_pts.T[0]
    df_decays['ip_y0'] = detector_entry_pts.T[1]
    df_decays['ip_z0'] = detector_entry_pts.T[2]

    df_decays['ip_x1'] = detector_exit_pts.T[0]
    df_decays['ip_y1'] = detector_exit_pts.T[1]
    df_decays['ip_z1'] = detector_exit_pts.T[2]

    # Save the intersection points for the inner detector
    df_decays['ip_inner_x0'] = detector_entry_pts_inner.T[0]
    df_decays['ip_inner_y0'] = detector_entry_pts_inner.T[1]
    df_decays['ip_inner_z0'] = detector_entry_pts_inner.T[2]

    df_decays['ip_inner_x1'] = detector_exit_pts_inner.T[0]
    df_decays['ip_inner_y1'] = detector_exit_pts_inner.T[1]
    df_decays['ip_inner_z1'] = detector_exit_pts_inner.T[2]

    dx = df_decays['x0']-df_decays['ip_x0']
    dy = df_decays['y0']-df_decays['ip_y0']
    dz = df_decays['z0']-df_decays['ip_z0']

    distance_to_detector = np.sqrt(dx*dx + dy*dy + dz*dz)
    df_decays['distance_to_detector'] = distance_to_detector
    

    #######################################################################
    # Do we save only the ones that strike the detector?
    #######################################################################
    detected_tag = ""
    if SAVE_ONLY_DETECTED:
        if verbose:
            print("\nSaving only the tracks that strike the detector")
            print(f"Initially dataframe is {df_decays.shape}")
        filter = df_decays['hit_detector'] == True
        df_temp = df_decays[filter]
        del df_decays
        df_decays = df_temp
        if verbose:
            print(f"After only keeping some tracks, dataframe is {df_decays.shape}\n")
        detected_tag = "_HIT_DETECTOR"
        
    #######################################################################
    # Calculate eloss
    #######################################################################

    # New stuff with precomputed splines
    start = time.time()

    e_final_vals = None

    if AVERAGE_ELOSS is False:

        '''
        # Load the spline object from the file
        with open(eloss_dict_name, 'rb') as file_handle:
            loaded_eloss_dict = pickle.load(file_handle)
        
        #print(f'Time to read energy loss spline file is {time.time() - start} seconds')

        lookup = eloss_tools.RaggedSplineIndexLookup(loaded_eloss_dict)

        ########### NEED TO DO THIS FOR BOTH MUONS!!!!!!!!! ######################
        #E_query = 10000 + 2000*np.random.random(nvals)
        #d_query = 1000*np.random.random(nvals)
        E_query = df_decays['e1']
        d_query = df_decays['distance_to_detector']

        #print(d_query)
        
        E_idx, D_idx = lookup.get_indices(E_query, d_query)

        print(f'Processing eloss for {len(E_query)} muons')
        start = time.time()
        
        uniq, counts = eloss_tools.count_pairs(E_idx, D_idx)    
        print(f'{len(uniq)=}, {len(counts)=}, {sum(counts)=}')
        
        eloss_data = {}
        for i,(u,c) in enumerate(zip(uniq, counts)):
            #print(u,c)
            if i%1000==0:
                print(i)
            ei,di = u
            #di = u[1]
            E_index = lookup.energy_keys[ei]
            d_index = lookup.dist_keys[ei][di]
            #print(u,c,E_index, d_index)
            spl = loaded_eloss_dict[E_index][d_index]
            de_vals = eloss_tools.generate_data_from_spline(spl, c)
            #print(de_vals)
            #eloss_data[(E_index,d_index)] = de_vals
            eloss_data[(ei,di)] = de_vals
        
        #print(f"Calculated all the uniq and counts")
        
        e_final_vals = []
        de = -1
        for i in range(len(E_query)):
            ei = E_idx[i]
            di = D_idx[i]
            # CHECK THIS BUT I THINK di is -1 when the distance is greater than the 
            # the max distance for that energy
            de_vals = eloss_data[(ei,di)]#.pop()
            #print(ei,di,de_vals)
            if len(de_vals)>0 and di>-1:
                de = de_vals.pop()
            else:
                # Make the delta E the same as the initial energy
                de = E_query.iloc[i]
            #print(i, len(E_query))
            efin = E_query.iloc[i] - de
            #print(E_query.iloc[i], di, de, efin)
            if efin<0:
                efin = 0
            e_final_vals.append(efin)
        '''

        E_query = df_decays['e1']
        d_query = df_decays['distance_to_detector']
        e_final_vals = eloss_tools.calculate_eloss_from_geant_splines(E_query, d_query, eloss_dict_name='eloss_dictionary_11032025_v2.pkl', verbose=verbose)


    else:
        E_query = df_decays['e1']
        d_query = df_decays['distance_to_detector']
        #e_final_vals = eloss_average.propagate_muon(E_query, d_query)
        e_final_vals = eloss_average.propagate_muon_CMSSW(E_query, d_query)

        
    if verbose:
        print(f'\nTime to calculate eloss for {len(E_query)} muons is {time.time() - start} seconds')

    df_decays['efinal_mu1'] = e_final_vals

    #######################################################################
    # Calculate pt as seen by the detector
    #######################################################################
    #pt = np.cos(theta) * df_decays['pmag1']
    # Get the intersection points back from the dataframe after reductions
    detector_entry_pts = np.array([df_decays['ip_x0'], df_decays['ip_y0'], df_decays['ip_z0']]).T 
    detector_exit_pts = np.array([df_decays['ip_x1'], df_decays['ip_y1'], df_decays['ip_z1']]).T 

    projection,vmag = dst.projection_length_on_plane_from_points(detector_entry_pts, detector_exit_pts, normal=(0.0, 0.0, 1.0), verbose=verbose)    
    #projection,vmag = dst.projection_length_on_plane_from_points(detector_entry_pts[hit_detector_idx], detector_exit_pts[hit_detector_idx])

    #print(projection)
    #print(vmag)
    ######## SHOULD BE CONSISTENT ABOUT USING PMAG OR ENERGY
    transverse_scaling = projection/vmag
    pt = transverse_scaling * df_decays['pmag1']
    pt_eloss = transverse_scaling * df_decays['efinal_mu1']
    
    #print(len(pt))
    #print(len(pt_eloss))
    
    df_decays['pt1_detector_acceptance'] = pt
    df_decays['pt1_detector_acceptance_eloss'] = pt_eloss
    
    #df_decays['pt1_scaling_detector_acceptance'] = transverse_scaling
    
    #df_decays['costh1_detector_acceptance'] = np.cos(theta)
    #df_decays['phi1_detector_acceptance'] = phi
    #df_decays['pmag1_detector_acceptance'] = pmag
    #'''

    ##########################################################################
    # Generate a tag that represents all this information
    depth_tag = dst.return_tag(depth)
    m_dm_tag = dst.return_tag(MASSES_DM)
    m_a_tag = dst.return_tag(MASSES_A)

    radius_tag = 'AD' # Angle dependent
    if disk_radius is not None:
        radius_tag = f'{disk_radius}'

    average_eloss_tag = '' 

    if AVERAGE_ELOSS is True:
        average_eloss_tag = f'_ave_eloss'

    tag = f'depth_{depth_tag}_diskR_{radius_tag}_mDM_{m_dm_tag}_mA_{m_a_tag}_dmModel_{dm_model}{detected_tag}{average_eloss_tag}{additional_tag}'

    # Do we write it all out to a file?
    if save_to_file:
        outfile = f'{output_directory}/generated_data_{tag}.parquet'
        df_decays.to_parquet(outfile)#, key='df')
        print(f'Saved to file {outfile}')
    ##########################################################################


    return df_decays, tag



