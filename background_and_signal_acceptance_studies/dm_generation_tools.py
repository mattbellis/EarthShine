import numpy as np
import pandas as pd
import matplotlib.pylab as plt

import seaborn as sns

from skspatial.objects import Line, Plane
from skspatial.plotting import plot_3d

from skspatial.objects import Line, Cylinder, Point
from skspatial.plotting import plot_3d

import detector_simulation_tools as dst

import phasespace

import tensorflow

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
#$def intersect_CMS(df):
#    #''

'''
def throw_muons_at_CMS(df_input, ndecays=None, MAKE_PLOTS=False):

    # Define CMS
    origin_CMS = [0, 0, 0]

    nmuons = 0

    # CMS, units are meters. x is direction of beam and z is up
    cylinder = Cylinder.from_points([-10.5, 0, 8], [10.5, 0, 8], 7.5)

    if MAKE_PLOTS:
    fig1 = plt.figure(figsize=(6,6))
    ax1 = fig1.add_subplot(1,1,1,projection='3d')

    fig2 = plt.figure(figsize=(12,12))
    ax2 = fig2.add_subplot(1,1,1,projection='3d')

    fig3 = plt.figure(figsize=(4,4))
    ax3 = fig3.add_subplot(1,1,1)

    # Draw CMS
    cylinder.plot_3d(ax2, alpha=0.2)

    # Which to use?
    #mask = (dfdec2['DM_MASS']==100)
    #mask = mask & (dfdec2['M_A']==5)
    #dftmp = dfdec2[mask]
    
    #mask = (dfdec1['pmag']==1000000)
    #mask = mask & (dfdec1['M_A']==5000)
    #dftmp = dfdec1[mask]

    #dftmp = df_input
    
    if ndecays is None:
        ndecays = len(df_input)
        print(f"Will run over {ndecays} decays")

    distances = []
    pmag_origin = []
    pmag_cms = []

    # Range over which to generate the interaction where the muon pairs are created
    limits = 100
    xlo, xhi = -limits,limits
    ylo, yhi = -limits,limits
    zlo, zhi = -limits, -1

    xwidth = xhi-xlo
    ywidth = yhi-ylo
    zwidth = zhi-zlo

    # Generate the points
    #npts = 3000
    #for i in dftmp.index[0:npts]:
    for i in range(0,ndecays):

        #print(i)

        if i%10000==0:
            print(i)

        # Generate the interaction point of the muon pairs
        origin = [xhi-xwidth*np.random.random(), yhi-ywidth*np.random.random(), -zhi-zwidth*np.random.random()]

        # For debugging
        #origin = [0, 0, -10]

        # Print all the origins
        #Point(origin).plot_3d(ax1, c='k')

        df = dftmp.iloc[i]
        pmag = None

        # Loop over the pairs
        for j in [0,1]:
            nmuons += 1
            # Need to mix up z and y
            if j==0:
                px1,py1,pz1 = df['px1'],df['py1'],df['pz1']
                pmag = df['pmag1']
                dir = np.array([px1, py1, pz1])
            else:
                px2,py2,pz2 = df['px_mu2'],df['py_mu2'],df['pz_mu2']
                pmag = df['pmag2']
                dir = np.array([px2, py2, pz2])

      # Skip downward traveling muons
      # Skip downward traveling muons
      if dir[2]<0:
        #print("skipping ",j,dir)
        continue

      # Need to do this to draw the lines correction
      m = mag(dir)
      dir /= m
      dir *= xwidth

      #print("origin")
      #print(origin)
      #print('dir')
      #print(dir)

      line = Line(point=origin, direction=dir) # Point and direction vector

      #print("here")
      point_a,point_b = None, None
      try:
        point_a, point_b = cylinder.intersect_line(line, infinite=False)
      except ValueError:
        1
        #print('does not')

      #print('which:  ',j)
      #print('dir:    ',dir)
      #print('origin: ', origin)
      #print('points: ',point_a, point_b)

      if point_a is not None:
        1
        #point_a.plot_3d(ax2, c='r',s=10)
        #line.plot_3d(ax2, c='k'),#, t_2=100),
        #Point(origin).plot_3d(ax2, c='b')
        #print(origin)
        #print(dir)

      if point_b is not None:
        1
        #point_b.plot_3d(ax2, c='g',s=10)

      if point_b is None or point_a is None:
        1
        #line.plot_3d(ax2, c='y', linestyle='--'),#, t_2=100),
        #Point(origin).plot_3d(ax2, c='y')
      else:
        d = distance(point_a, origin)
        distances.append(d)
        pmag_origin.append(pmag)
        pfinal = eloss_interp(pmag,d)[0][0]
        pmag_cms.append(pfinal)

  if MAKE_PLOTS:
    ax2.set_xlim(-20,20)
    ax2.set_ylim(-20,20)
    ax2.set_zlim(-20,20)
    ax2.set_xlabel('x')
    ax2.set_ylabel('y')
    ax2.set_zlabel('z')

    plt.sca(ax3)
    plt.hist(distances);

  nhits = len(distances)
  print(f"{nhits}  {ndecays}   {nhits/ndecays:.6f}")

  return pmag_origin, pmag_cms, distances, nmuons

'''

