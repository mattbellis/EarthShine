import numpy as np
import pandas as pd
import matplotlib.pylab as plt

import seaborn as sns

from skspatial.objects import Line, Plane
from skspatial.plotting import plot_3d

from skspatial.objects import Line, Cylinder, Point
from skspatial.plotting import plot_3d

import phasespace

import tensorflow

import pickle

################################################################################
def mag(p3):
  #print(p3)
  p = np.sqrt(p3[0]*p3[0] + p3[1]*p3[1] + p3[2]*p3[2])
  #if p<5000:
  #  print(p3,p)
  return p

def invmass(p4s):
  E,px,py,pz = 0,0,0,0

  for p4 in p4s:
    E += p4[3]
    px += p4[0]
    py += p4[1]
    pz += p4[2]

  m2 = E**2 - (px**2 + py**2 + pz**2)
  if m2>=0:
    return np.sqrt(m2)
  else:
    return -np.sqrt(-m2)

def invmass_cols(p4):

  E = p4[3]
  px = p4[0]
  py = p4[1]
  pz = p4[2]

  m2 = E**2 - (px**2 + py**2 + pz**2)
  m = -999*np.ones(len(E))
  mask = m2>=0
  #print(mask[mask])
  #print(mask[~mask])

  m[mask] = np.sqrt(m2[mask])
  m[~mask] = -np.sqrt(-m2[~mask])

  return m


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

def generate_dm_decays(MASS_A=[.250,1,5], DM_MASSES=[10,100,1000], nevents_to_generate=10):
    # Loop over different values

    decays = {}
    decays['M_DM'] = []
    decays['M_A'] = []

    decays['px_mu1'] = []
    decays['py_mu1'] = []
    decays['pz_mu1'] = []
    decays['e_mu1'] = []

    decays['px_mu2'] = []
    decays['py_mu2'] = []
    decays['pz_mu2'] = []
    decays['e_mu2'] = []

    #MASS_A = [250,1000,5000]
    #MUON_MASS = 105.11

    #MASS_A = [.250,1,5]
    #DM_MASSES = [10,100,1000]

    MUON_MASS = 0.10511

    thetas = []
    pmags = []

    muon_p = []

    #nevents_to_generate = 1000

    #pmags_GeV = [10,100,1000]

    for MASS_A in MASS_A:
      for pmag in DM_MASSES:

        #pmag = pmag_GeV*1000

        print(MASS_A,pmag)

        # Here is the boost vectors to boost it to the moving frame of the dark photon, pmag
        boost_vector = np.array([0,0, pmag, np.sqrt(pmag**2 + MASS_A**2)])
        boost_vectors = np.tile(boost_vector, (nevents_to_generate,1))

        print(f"Generating {nevents_to_generate} decays")
        weights, particles = phasespace.nbody_decay(MASS_A,
                                                    [MUON_MASS, MUON_MASS]).generate(n_events=nevents_to_generate, boost_to=boost_vectors)

        print(f"Generated the decays")

        print(f"Pulling out the values from the decays and filling our dataframe...")

        p0 = particles['p_0'][:].numpy().T
        p1 = particles['p_1'][:].numpy().T

        #data[MASS_A][pmag] = [p0.T, p1.T]
        decays['M_DM'] += (pmag*np.ones(nevents_to_generate)).tolist()
        decays['M_A'] += (MASS_A*np.ones(nevents_to_generate)).tolist()

        #print(len(decays['M_A']))

        decays['px_mu1'] += p0[0].tolist()
        decays['py_mu1'] += p0[1].tolist()
        decays['pz_mu1'] += p0[2].tolist()
        decays['e_mu1'] +=  p0[3].tolist()

        decays['px_mu2'] += p1[0].tolist()
        decays['py_mu2'] += p1[1].tolist()
        decays['pz_mu2'] += p1[2].tolist()
        decays['e_mu2'] +=  p1[3].tolist()
        print(f"Pulled out the values from the decays and filled our dataframe...")

    dfdec1 = pd.DataFrame.from_dict(decays)

    # Generate a few other entries
    print("Generated all the events and now calculating a few new values!")

    # pmag, theta (degrees), phi
    px1 = dfdec1['px_mu1'].values
    py1 = dfdec1['py_mu1'].values
    pz1 = dfdec1['pz_mu1'].values
    e1 = dfdec1['e_mu1'].values

    pmag1 = mag([px1,py1,pz1])
    theta1 = np.rad2deg(np.arccos(pz1/pmag1))
    phi1 = np.arctan2(py1,pz1)
    print("Calculated pmag1, theta1, phi1")

    dfdec1['pmag1'] = pmag1
    dfdec1['theta1'] = theta1
    dfdec1['phi1'] = phi1

    px2 = dfdec1['px_mu2'].values
    py2 = dfdec1['py_mu2'].values
    pz2 = dfdec1['pz_mu2'].values
    e2 = dfdec1['e_mu2'].values

    pmag2 = mag([px2,py2,pz2])
    theta2 = np.rad2deg(np.arccos(pz2/pmag2))
    phi2 = np.arctan2(py2,pz2)
    print("Calculated pmag2, theta2, phi2")

    dfdec1['pmag2'] = pmag2
    dfdec1['theta2'] = theta2
    dfdec1['phi2'] = phi2

    # Opening angle
    p4s = [[px1,py1,pz1,e1], [px2,py2,pz2,e2]]
    thetas = np.rad2deg(opening_angle(p4s))
    dfdec1['opening angle'] = thetas
    print("Calculated opening values")

    #dfdec1.sample(5)

    return dfdec1



################################################################################
#$def intersect_CMS(df):
#    #'''

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
                px1,py1,pz1 = df['px_mu1'],df['py_mu1'],df['pz_mu1']
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


################################################################################

# Got this code from ChatGPT and verified the output with the skspatial tools.
# 
# This code runs much faster though because I believe it is specialized for this
# purpose while the skspatial tools seem to be more generalized. 
def intersect_finite_cylinder_x_np(origins, directions,
                                   radius=7.5, half_len=10.5, eps=1e-12):
    """
    Vectorized intersection of N rays with a finite cylinder along the x-axis.

    Parameters
    ----------
    origins : (N,3) array_like
        Ray start points.
    directions : (N,3) array_like
        Ray direction vectors.
    radius : float
        Cylinder radius.
    half_len : float
        Half the length of the cylinder along x (so x ∈ [-half_len, +half_len]).
    eps : float
        Threshold for treating a coefficient as zero.

    Returns
    -------
    pts0, pts1 : each an (N,3) array
        The first and second intersection points.  If a ray has
        <1 intersection, that row is NaN; if exactly 1, pts1 is NaN.
    """
    O = np.asarray(origins, dtype=float)
    D = np.asarray(directions, dtype=float)
    Ox, Oy, Oz = O[:,0], O[:,1], O[:,2]
    Dx, Dy, Dz = D[:,0], D[:,1], D[:,2]
    N = len(O)

    # --- barrel (side) intersections ---
    a = Dy**2 + Dz**2
    b = 2*(Oy*Dy + Oz*Dz)
    c = Oy**2 + Oz**2 - radius**2

    disc = b*b - 4*a*c
    real = (disc >= 0) & (a > eps)
    sqrt_disc = np.sqrt(np.clip(disc, 0, None))
    inv2a   = 0.5 / np.where(a>eps, a, 1.0)    # avoid div0

    # two roots
    t_barrel0 = (-b - sqrt_disc) * inv2a
    t_barrel1 = (-b + sqrt_disc) * inv2a

    # keep only those within x‐slab
    x0 = Ox + t_barrel0*Dx
    x1 = Ox + t_barrel1*Dx
    ok0 = real & (x0 >= -half_len) & (x0 <= half_len)
    ok1 = real & (x1 >= -half_len) & (x1 <= half_len)

    # --- cap intersections ---
    # avoid division by zero
    nonpara = np.abs(Dx) > eps
    t_cap_pos = np.where(nonpara, ( half_len - Ox)/Dx, np.nan)
    t_cap_neg = np.where(nonpara, (-half_len - Ox)/Dx, np.nan)

    # check disk‐inclusion
    y_pos = Oy + t_cap_pos*Dy
    z_pos = Oz + t_cap_pos*Dz
    y_neg = Oy + t_cap_neg*Dy
    z_neg = Oz + t_cap_neg*Dz

    ok_pos = nonpara & (y_pos*y_pos + z_pos*z_pos <= radius*radius)
    ok_neg = nonpara & (y_neg*y_neg + z_neg*z_neg <= radius*radius)

    # --- stack all candidates ---
    # shape (N,4)
    t_cand = np.stack([t_barrel0, t_barrel1, t_cap_pos, t_cap_neg], axis=1)
    valid  = np.stack([ok0,        ok1,        ok_pos,    ok_neg],   axis=1)

    # mask out invalids to NaN (so they sort to the end)
    t_cand = np.where(valid, t_cand, np.nan)

    # --- pick the two smallest t values per ray ---
    order = np.argsort(t_cand, axis=1)           # NaNs go last
    idx0  = order[:, 0]
    idx1  = order[:, 1]

    # gather t0, t1
    t0 = t_cand[np.arange(N), idx0]
    t1 = t_cand[np.arange(N), idx1]

    # --- recover 3D points; NaNs propagate automatically ---
    P0 = O + t0[:,None] * D
    P1 = O + t1[:,None] * D

    return P0, P1

################################################################################
################################################################################
def draw_points(point_a, point_b, ax=plt.gca(), color='r'):

    pa,pb,line = None, None, None
    if point_a[0]==point_a[0]:
        pa = Point(point_a)
        pa.plot_3d(ax, s=75, c=color)

    if point_b[0]==point_b[0]:
        pb = Point(point_b)
        pb.plot_3d(ax, s=75, c=color)

    if point_b[0]==point_b[0] and point_a[0]==point_a[0]:
        #line = Line(point_a, point_b)
        line = Line(pa, pb-pa)

        line.plot_3d(ax, c=color)

    return pa,pb,line

################################################################################
