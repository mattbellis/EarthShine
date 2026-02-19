import numpy as np
import pandas as pd
import matplotlib.pylab as plt

from skspatial.objects import Line, Plane
from skspatial.plotting import plot_3d

from skspatial.objects import Line, Cylinder, Point, Points
from skspatial.plotting import plot_3d



##########################################################################
# CMS coordinate transformations
#
# From Claude
##########################################################################
import numpy as np

################################################################################
def cartesian_to_cms(px: np.ndarray, py: np.ndarray, pz: np.ndarray) -> dict:
    """
    Convert Cartesian momentum (px, py, pz) to CMS coordinates (pt, theta, phi, eta).

    CMS coordinate conventions:
    - z-axis: beam direction
    - x-axis: horizontal, pointing to center of LHC
    - y-axis: vertical, pointing up
    - theta: polar angle from +z axis [0, pi]
    - phi: azimuthal angle in x-y plane from +x axis [-pi, +pi]
    - eta: pseudorapidity = -ln(tan(theta/2))
    - pt: transverse momentum = sqrt(px^2 + py^2)

    Parameters
    ----------
    px, py, pz : array-like
        Cartesian momentum components (can be scalars or arrays)

    Returns
    -------
    dict with keys:
        'pt': transverse momentum
        'theta': polar angle in radians [0, pi]
        'phi': azimuthal angle in radians [-pi, pi]
        'eta': pseudorapidity
        'p': total momentum magnitude
    """
    px = np.atleast_1d(np.asarray(px, dtype=float))
    py = np.atleast_1d(np.asarray(py, dtype=float))
    pz = np.atleast_1d(np.asarray(pz, dtype=float))

    # Transverse momentum
    pt = np.sqrt(px**2 + py**2)

    # Total momentum
    p = np.sqrt(px**2 + py**2 + pz**2)

    # Azimuthal angle phi: angle in x-y plane from x-axis
    # atan2 returns values in [-pi, pi]
    phi = np.arctan2(py, px)

    # Polar angle theta: angle from z-axis
    # theta = arctan(pt / pz), but we use arctan2 for proper quadrant handling
    # For a particle along +z: theta = 0
    # For a particle along -z: theta = pi
    # For a particle in the transverse plane: theta = pi/2
    theta = np.arctan2(pt, pz)  # This gives [0, pi] range as desired

    # Pseudorapidity eta = -ln(tan(theta/2))
    # Equivalent formula that's numerically more stable:
    # eta = arctanh(pz/p) = 0.5 * ln((p + pz)/(p - pz))
    # Or: eta = -ln(tan(theta/2))

    # Handle edge cases where pt = 0 (particle along beam axis)
    with np.errstate(divide='ignore', invalid='ignore'):
        # Use the stable formula: eta = arctanh(cos(theta)) = arctanh(pz/p)
        # But arctanh diverges at ±1, so we use the log formula with protection
        eta = np.where(
            pt > 1e-10,
            -np.log(np.tan(theta / 2)),
            np.sign(pz) * np.inf  # eta -> ±inf along beam axis
        )

        # Alternative stable computation using momentum components directly
        # eta = 0.5 * np.log((p + pz) / (p - pz))  # Same result

    #return {
        #'pt': pt,
        #'theta': theta,
        #'phi': phi,
        #'eta': eta,
        #'p': p
    #}
    return p, pt, eta, phi, theta

################################################################################

def cms_to_cartesian(
    pt: np.ndarray = None,
    phi: np.ndarray = None,
    eta: np.ndarray = None,
    theta: np.ndarray = None,
    p: np.ndarray = None
) -> dict:
    """
    Convert CMS coordinates to Cartesian (px, py, pz).

    Two input modes:
    1. (pt, phi, eta) or (pt, phi, theta) - transverse momentum + angles
    2. (p, phi, theta) - total momentum + angles

    Parameters
    ----------
    pt : array-like, optional
        Transverse momentum (use with phi and eta/theta)
    phi : array-like
        Azimuthal angle in radians [-pi, pi] (required)
    eta : array-like, optional
        Pseudorapidity (use with pt)
    theta : array-like, optional
        Polar angle in radians [0, pi]
    p : array-like, optional
        Total momentum magnitude (use with phi and theta)

    Returns
    -------
    dict with keys:
        'px': x-component of momentum
        'py': y-component of momentum
        'pz': z-component of momentum
        'p': total momentum magnitude
        'pt': transverse momentum
        'theta': polar angle
        'eta': pseudorapidity
    """
    # Validate inputs
    if phi is None:
        raise ValueError("'phi' is required")

    phi = np.atleast_1d(np.asarray(phi, dtype=float))

    # Mode 1: (p, phi, theta) - total momentum with angles
    if p is not None:
        if theta is None:
            raise ValueError("When using 'p', must also provide 'theta'")
        if pt is not None:
            raise ValueError("Provide either 'pt' or 'p', not both")
        if eta is not None:
            raise ValueError("When using 'p', use 'theta' not 'eta'")

        p = np.atleast_1d(np.asarray(p, dtype=float))
        theta = np.atleast_1d(np.asarray(theta, dtype=float))

        # Standard spherical to Cartesian
        px = p * np.sin(theta) * np.cos(phi)
        py = p * np.sin(theta) * np.sin(phi)
        pz = p * np.cos(theta)
        pt = p * np.sin(theta)

        # Compute eta from theta
        with np.errstate(divide='ignore', invalid='ignore'):
            eta = np.where(
                np.abs(np.sin(theta)) > 1e-10,
                -np.log(np.tan(theta / 2)),
                np.sign(np.cos(theta)) * np.inf
            )

    # Mode 2: (pt, phi, eta/theta) - transverse momentum with angles
    elif pt is not None:
        if (eta is None) == (theta is None):
            raise ValueError("When using 'pt', must provide exactly one of 'eta' or 'theta'")

        pt = np.atleast_1d(np.asarray(pt, dtype=float))

        if eta is not None:
            eta = np.atleast_1d(np.asarray(eta, dtype=float))
            # Convert eta to theta: theta = 2 * arctan(exp(-eta))
            theta = 2 * np.arctan(np.exp(-eta))
        else:
            theta = np.atleast_1d(np.asarray(theta, dtype=float))
            # Convert theta to eta
            with np.errstate(divide='ignore', invalid='ignore'):
                eta = np.where(
                    np.abs(np.sin(theta)) > 1e-10,
                    -np.log(np.tan(theta / 2)),
                    np.sign(np.cos(theta)) * np.inf
                )

        # Cartesian components
        px = pt * np.cos(phi)
        py = pt * np.sin(phi)

        # z-component: pz = pt / tan(theta) = pt * cot(theta)
        with np.errstate(divide='ignore', invalid='ignore'):
            pz = np.where(
                np.abs(np.sin(theta)) > 1e-10,
                pt / np.tan(theta),
                np.sign(np.cos(theta)) * np.inf
            )

        # Total momentum
        p = np.sqrt(px**2 + py**2 + pz**2)

    else:
        raise ValueError("Must provide either 'pt' or 'p'")

    #return {
        #'px': px,
        #'py': py,
        #'pz': pz,
        #'p': p,
        #'pt': pt,
        #'theta': theta,
        #'eta': eta
    #}

    return px, py, pz, p, pt, theta, eta


################################################################################
# Convenience wrappers for vector inputs
def xyz_to_cms(momentum_xyz: np.ndarray) -> dict:
    """
    Convert momentum vectors from Cartesian to CMS coordinates.

    Parameters
    ----------
    momentum_xyz : ndarray of shape (N, 3) or (3,)
        Momentum vectors as [px, py, pz]

    Returns
    -------
    dict with 'pt', 'theta', 'phi', 'eta', 'p' arrays
    """
    momentum_xyz = np.atleast_2d(momentum_xyz)
    return cartesian_to_cms(
        momentum_xyz[:, 0],
        momentum_xyz[:, 1],
        momentum_xyz[:, 2]
    )


################################################################################

def cms_to_xyz(pt: np.ndarray, phi: np.ndarray,
               eta: np.ndarray = None, theta: np.ndarray = None) -> np.ndarray:
    """
    Convert CMS coordinates to momentum vectors in Cartesian coordinates.

    Parameters
    ----------
    pt : array-like
        Transverse momentum
    phi : array-like
        Azimuthal angle in radians
    eta : array-like, optional
        Pseudorapidity (provide this OR theta)
    theta : array-like, optional
        Polar angle in radians (provide this OR eta)

    Returns
    -------
    ndarray of shape (N, 3)
        Momentum vectors as [px, py, pz]
    """
    result = cms_to_cartesian(pt, phi, eta=eta, theta=theta)
    return np.column_stack([result['px'], result['py'], result['pz']])

##########################################################################
def mag(p3):
  #print(p3)
  p = np.sqrt(p3[0]*p3[0] + p3[1]*p3[1] + p3[2]*p3[2])
  #if p<5000:
  #  print(p3,p)
  return p

##########################################################################
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

##########################################################################
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


##########################################################################
def opening_angle(p4s):
  p0mag = np.sqrt(p4s[0][0]**2 + p4s[0][1]**2 + p4s[0][2]**2)
  p1mag = np.sqrt(p4s[1][0]**2 + p4s[1][1]**2 + p4s[1][2]**2)

  dot_product = p4s[0][0]*p4s[1][0] + p4s[0][1]*p4s[1][1] + p4s[0][2]*p4s[1][2]

  theta = np.arccos(dot_product/(p0mag*p1mag))
  return theta


##########################################################################
def distance(v1, v2):
  dx = v1[0] - v2[0]
  dy = v1[1] - v2[1]
  dz = v1[2] - v2[2]

  d = np.sqrt(dx**2 + dy**2 + dz**2)
  return d

##################################################################

def draw_detector(pt0=[0, 0, -15], pt1=[0, 0, 15], length=7.5):
    # Define detector
    #origin_detector = [0, 0, 0]

    nmuons = 0

    # detector, units are meters. x is direction of beam and z is up
    #cylinder = Cylinder.from_points([-10.5, 0, 0], [10.5, 0, 0], 7.5)
    # Mock detectors
    #cylinder = Cylinder.from_points([0, 0, -15], [0, 0, 15], 7.5)
    cylinder = Cylinder.from_points(pt0, pt1, length)

    #if MAKE_PLOTS:
    #fig1 = plt.figure(figsize=(6,6))
    #ax1 = fig1.add_subplot(1,1,1,projection='3d')

    fig = plt.figure(figsize=(12,6))
    ax = fig.add_subplot(1,1,1,projection='3d')

    #fig3 = plt.figure(figsize=(4,4))
    #ax3 = fig3.add_subplot(1,1,1)

    # Draw detector
    cylinder.plot_3d(ax, alpha=0.2)
    ax.set_xlim(-100,100)
    ax.set_ylim(-100,100)
    ax.set_zlim(-100,100)

    return fig, ax, cylinder

##########################################################################
################################################################################

# Got this code from ChatGPT and verified the output with the skspatial tools.
# 
# This code runs much faster though because I believe it is specialized for this
# purpose while the skspatial tools seem to be more generalized. 
def intersect_finite_cylinder_x_np(origins, directions,
                                   radius=8, half_len=10.5, eps=1e-12):
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
    intersection_point0, intersection_point1 : each an (N,3) array
        The first and second intersection points.  If a ray has
        <1 intersection, that row is NaN; if exactly 1, intersection_point1 is NaN.
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
    intersection_point0 = O + t0[:,None] * D
    intersection_point1 = O + t1[:,None] * D

    return intersection_point0, intersection_point1

################################################################################
################################################################################
def draw_points(point_a, point_b, ax=plt.gca(), color='r'):

    pa,pb,line = None, None, None
    if point_a[0]==point_a[0]:
        pa = Point(point_a)
        pa.plot_3d(ax, s=75, c=color)

    if point_b[0]==point_b[0]:
        pb = Point(point_b)
        pb.plot_3d(ax, s=75, c=color, alpha=0.2)

    if point_b[0]==point_b[0] and point_a[0]==point_a[0]:
        #line = Line(point_a, point_b)
        line = Line(pa, pb-pa)

        line.plot_3d(ax, c=color)

    return pa,pb,line

################################################################################

##########################################################################
def is_iterable(obj):
    try:
        iter(obj)
        return True
    except TypeError:
        return False
##########################################################################

##########################################################################
def return_tag(var):
    var_tag = ""
    if type(var)==int or type(var)==float or type(var)==np.float64:
        var_tag = f'{var}'
    elif is_iterable(var):
        if len(var)>1:
            var_tag = f'{var[0]}-{var[-1]}'
        else:
            var_tag = f'{var[0]}'

    return var_tag

##########################################################################

##########################################################################
def generate_origins(nevents=100, disk_radius=None, depth=0):
    # Generate random origins in a circule
    origin_phi = 2*np.pi*np.random.random(nevents)

    if disk_radius is None:
        angle = np.deg2rad(91)
        disk_radius = int(np.ceil(np.abs(np.abs(depth)*np.tan(angle))))

        if disk_radius>6000:
            disk_radius = 6000
        
    print(f'Origin points will be uniformly scattered on a plane with a disk radius: {disk_radius} m')
    origin_radial_coord = disk_radius*np.sqrt(np.random.random(nevents))

    #x = origin_r*np.cos(origin_phi)
    #y = origin_r*np.sin(origin_phi)
    #z = None
    # Switch to CMS coordinates
    x = origin_radial_coord*np.cos(origin_phi)
    z = origin_radial_coord*np.sin(origin_phi)
    y = None

    if type(depth)==int or type(depth)==float or type(depth)==np.float64:
        y = depth*np.ones(nevents)
    elif is_iterable(depth):
        width = depth[1] - depth[0]
        y = depth[0] + (width*np.random.random(nevents))
    else:
        print(f"depth variable {depth} is not int,float, or iterable")
        print(f'{depth} {type(depth)}')
        return 0

    origins = np.array([x,y,z]).T

    return origins, disk_radius

##########################################################################
# Extract the momenta and calculate unit vectors for the directions
# of travel of the muons:wq
##########################################################################
def directions_from_momentum(df_decays):

    directions = []
    for i in [1,2]:
        px = df_decays[f'px{i}']
        py = df_decays[f'py{i}']
        pz = df_decays[f'pz{i}']

        #pmag = np.sqrt(px*px + py*py + pz*pz)

        #px /= pmag
        #py /= pmag
        #pz /= pmag

        directions.append(np.array([px, py, pz]).T)

    return directions
##########################################################################


##########################################################################
##########################################################################
def generate_origins_and_directions(nevents=100, radius=10, depth=-7.5, dm_model='floating'):
    # Generate random origins in a circule
    origin_phi = 2*np.pi*np.random.random(nevents)
    origin_r = radius*np.sqrt(np.random.random(nevents))
    #x = origin_r*np.cos(origin_phi)
    #y = origin_r*np.sin(origin_phi)
    #z = None

    # CMS coordinates
    x = origin_r*np.cos(origin_phi)
    z = origin_r*np.sin(origin_phi)
    y = None

    if type(depth)==int or type(depth)==float:
        z = depth*np.ones(nevents)
    elif is_iterable(depth):
        width = depth[1] - depth[0]
        z = depth[0] + (width*np.random.random(nevents))
    else:
        print(f"depth variable {depth} is not int,float, or iterable")
        return 0
    origins = np.array([x,y,z]).T

    # Generate momenta 
    # Generate betwen 0 and 1
    theta = None
    if dm_model == 'floating':   
        # Uniform in all directions
        theta = np.arccos(np.random.random(nevents))
    elif dm_model == 'core':
        # Coming up from underground
        theta = 0*np.ones(nevents)
    phi = 2*np.pi*np.random.random(nevents)
    # Make sure the line is long enough to intersect the cylinder
    r = 1000*radius*np.ones(nevents)

    px = r*np.sin(theta)*np.cos(phi)
    py = r*np.sin(theta)*np.sin(phi)
    pz = r*np.cos(theta)

    directions = np.array([px, py, pz]).T

    return origins, directions

##########################################################################
# Developed with ChatGPT to project vector onto endcap plane
def projection_length_on_plane_from_points(detector_entry_pts, detector_exit_pts, normal=(0.0, 0.0, 1.0)):
# def projection_length_on_plane_from_points(detector_entry_pts, detector_exit_pts, normal=(1.0, 0.0, 0.0)):
    """
    Vectorized: magnitude of the projection(s) of segments detector_entry_pts->detector_exit_pts onto the plane with given normal.

    Parameters
    ----------
    detector_entry_pts, detector_exit_pts : array-like, shape (..., 3)
        Each row (or last-dimension triple) is a 3D point. detector_entry_pts and detector_exit_pts must be broadcastable.
    normal : array-like, shape (3,)
        Plane normal (need not be unit). The same normal is used for all rows.

    Returns
    -------
    lengths : ndarray, shape (...,)
        Magnitude(s) of the projected vector(s) onto the plane.
    """
    detector_entry_pts = np.asarray(detector_entry_pts, dtype=float)
    detector_exit_pts = np.asarray(detector_exit_pts, dtype=float)
    normal = np.asarray(normal, dtype=float)
    normal = normal / np.linalg.norm(normal)  # normalize

    particle_direction = detector_exit_pts - detector_entry_pts                                    # (..., 3)
    vmag = np.sqrt((particle_direction*particle_direction).sum(axis=1))
    
    print(f"[projection] number of segments = {particle_direction.shape[0]}")
    print(f"[projection] particle_direction = detector_exit_pts - detector_entry_pts, first 3 rows:\n{particle_direction[:3]}")
    print(f"[projection] |v| (full 3D magnitude), first 3: {vmag[:3]}")

    v_dot_n = np.sum(particle_direction * normal, axis=-1, keepdims=True)  # (..., 1)
    v_plane = particle_direction - v_dot_n * normal                      # (..., 3)
    proj_lengths = np.linalg.norm(v_plane, axis=-1)
    print(f"[projection] v·n̂ (component along normal), first 3: {v_dot_n[:3].flatten()}")
    print(f"[projection] v_plane (in-plane component), first 3 rows:\n{v_plane[:3]}")
    print(f"[projection] |v_plane| (projected length), first 3: {proj_lengths[:3]}")
    print(f"[projection] ratio |v_plane|/|v| (how much is in-plane), first 3: {(proj_lengths[:3] / vmag[:3])}")

    return proj_lengths, vmag       # (...,)

def projection_length_on_cylinder_base(detector_entry_pts, detector_exit_pts):
    """
    Specialized for cylinder bases perpendicular to x-axis.
    Equivalent to projecting onto any plane parallel to the end-caps.
    """
    detector_entry_pts = np.asarray(detector_entry_pts, dtype=float)
    detector_exit_pts = np.asarray(detector_exit_pts, dtype=float)
    dy = detector_exit_pts[..., 1] - detector_entry_pts[..., 1]
    dz = detector_exit_pts[..., 2] - detector_entry_pts[..., 2]
    return np.hypot(dy, dz)


# Example
#detector_entry_pts = (-0.3, 0.2, -0.4)
#detector_exit_pts = ( 0.4, 0.6,  0.8)
#print(projection_length_on_cylinder_base(detector_entry_pts, detector_exit_pts))                 # -> 1.264911064...
#print(projection_length_on_plane_from_points(detector_entry_pts, detector_exit_pts))             # same result

##########################################################################


################################################################################
# Claude
################################################################################

################################################################################
def ray_cylinder_intersection(
    origins: np.ndarray,      # Shape (N, 3) - starting points
    directions: np.ndarray,   # Shape (N, 3) - momentum/direction vectors (not necessarily normalized)
    radius: float,            # Cylinder radius
    half_length: float        # Half the cylinder length (extends from -half_length to +half_length on z-axis)
) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute entry and exit points for rays intersecting a cylinder centered at origin,
    aligned along the z-axis.
    
    Parameters
    ----------
    origins : ndarray of shape (N, 3)
        Starting points of the rays (outside the cylinder)
    directions : ndarray of shape (N, 3)
        Direction vectors (momentum vectors) of the rays
    radius : float
        Radius of the cylinder
    half_length : float
        Half-length of the cylinder (z ranges from -half_length to +half_length)
    
    Returns
    -------
    entry_points : ndarray of shape (N, 3)
        Entry points into the cylinder (NaN for rays that miss)
    exit_points : ndarray of shape (N, 3)
        Exit points from the cylinder (NaN for rays that miss)
    """
    origins = np.atleast_2d(origins)
    directions = np.atleast_2d(directions)
    
    N = origins.shape[0]
    
    # Initialize output arrays with NaN
    entry_points = np.full((N, 3), np.nan)
    exit_points = np.full((N, 3), np.nan)
    
    # Extract components
    ox, oy, oz = origins[:, 0], origins[:, 1], origins[:, 2]
    dx, dy, dz = directions[:, 0], directions[:, 1], directions[:, 2]
    
    # --- Curved surface intersection ---
    # Solve (ox + t*dx)² + (oy + t*dy)² = R²
    # This is a quadratic: a*t² + b*t + c = 0
    a = dx**2 + dy**2
    b = 2 * (ox * dx + oy * dy)
    c = ox**2 + oy**2 - radius**2
    
    discriminant = b**2 - 4 * a * c
    
    # For each ray, we'll collect valid t values
    # We need to handle curved surface and end caps separately
    
    t_candidates = np.full((N, 4), np.inf)  # Up to 4 candidate t values per ray
    
    # Curved surface intersections (where discriminant >= 0 and a != 0)
    valid_curved = (discriminant >= 0) & (np.abs(a) > 1e-12)
    sqrt_disc = np.sqrt(np.maximum(discriminant, 0))
    
    t1 = np.where(valid_curved, (-b - sqrt_disc) / (2 * a), np.inf)
    t2 = np.where(valid_curved, (-b + sqrt_disc) / (2 * a), np.inf)
    
    # Check if curved surface intersections are within z bounds
    z1 = oz + t1 * dz
    z2 = oz + t2 * dz
    
    t1_valid = valid_curved & (t1 >= 0) & (np.abs(z1) <= half_length)
    t2_valid = valid_curved & (t2 >= 0) & (np.abs(z2) <= half_length)
    
    t_candidates[:, 0] = np.where(t1_valid, t1, np.inf)
    t_candidates[:, 1] = np.where(t2_valid, t2, np.inf)
    
    # --- End cap intersections ---
    # Top cap: z = +half_length, solve oz + t*dz = half_length
    # Bottom cap: z = -half_length
   
    # Avoid division by zero
    dz_safe = np.where(np.abs(dz) > 1e-12, dz, np.inf)
    
    t_top = (half_length - oz) / dz_safe
    t_bottom = (-half_length - oz) / dz_safe
    
    # Check if cap intersections are within radius
    x_top = ox + t_top * dx
    y_top = oy + t_top * dy
    x_bottom = ox + t_bottom * dx
    y_bottom = oy + t_bottom * dy
    
    r2_top = x_top**2 + y_top**2
    r2_bottom = x_bottom**2 + y_bottom**2
    
    t_top_valid = (t_top >= 0) & (r2_top <= radius**2) & (np.abs(dz) > 1e-12)
    t_bottom_valid = (t_bottom >= 0) & (r2_bottom <= radius**2) & (np.abs(dz) > 1e-12)
    
    t_candidates[:, 2] = np.where(t_top_valid, t_top, np.inf)
    t_candidates[:, 3] = np.where(t_bottom_valid, t_bottom, np.inf)
    
    # --- Find the two smallest positive t values (entry and exit) ---
    t_sorted = np.sort(t_candidates, axis=1)
    
    t_entry = t_sorted[:, 0]
    t_exit = t_sorted[:, 1]
    
    # Rays that hit have finite entry and exit times
    valid_hit = np.isfinite(t_entry) & np.isfinite(t_exit)
    
    # Compute intersection points
    entry_points[valid_hit] = origins[valid_hit] + t_entry[valid_hit, np.newaxis] * directions[valid_hit]
    exit_points[valid_hit] = origins[valid_hit] + t_exit[valid_hit, np.newaxis] * directions[valid_hit]
    
    return entry_points, exit_points


################################################################################

def visualize_ray_cylinder_intersections(
    origins: np.ndarray,
    directions: np.ndarray,
    entry_points: np.ndarray,
    exit_points: np.ndarray,
    radius: float,
    half_length: float,
    ax: plt.Axes = None,
    cylinder_alpha: float = 0.2,
    cylinder_color: str = 'cyan',
    ray_colors: list = None,
    figsize: tuple = (12, 10),
    extend_factor: float = 1.5,
) -> plt.Axes:
    """
    Visualize rays intersecting a cylinder.

    Parameters
    ----------
    origins : ndarray of shape (N, 3)
        Starting points of the rays
    directions : ndarray of shape (N, 3)
        Direction vectors of the rays
    entry_points : ndarray of shape (N, 3)
        Entry points into the cylinder (can contain NaN for misses)
    exit_points : ndarray of shape (N, 3)
        Exit points from the cylinder (can contain NaN for misses)
    radius : float
        Radius of the cylinder
    half_length : float
        Half-length of the cylinder
    ax : matplotlib 3D axes, optional
        Existing axes to plot on. If None, creates new figure.
    cylinder_alpha : float
        Transparency of the cylinder surface
    cylinder_color : str
        Color of the cylinder
    ray_colors : list, optional
        List of colors for each ray. If None, uses a colormap.
    figsize : tuple
        Figure size if creating new figure
    extend_factor : float
        How far to extend the dashed ray line beyond the cylinder

    Returns
    -----
    ax : matplotlib 3D axes
    """
    origins = np.atleast_2d(origins)
    directions = np.atleast_2d(directions)
    entry_points = np.atleast_2d(entry_points)
    exit_points = np.atleast_2d(exit_points)

    N = origins.shape[0]

    if ax is None:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, projection='3d')

    # Generate ray colors if not provided
    if ray_colors is None:
        cmap = plt.cm.tab10
        ray_colors = [cmap(i % 10) for i in range(N)]

    # --- Draw the cylinder ---
    # Curved surface
    theta = np.linspace(0, 2 * np.pi, 50)
    z_cyl = np.linspace(-half_length, half_length, 30)
    theta_grid, z_grid = np.meshgrid(theta, z_cyl)
    x_cyl = radius * np.cos(theta_grid)
    y_cyl = radius * np.sin(theta_grid)

    ax.plot_surface(x_cyl, y_cyl, z_grid, alpha=cylinder_alpha,
                    color=cylinder_color, edgecolor='none')

    # End caps
    r_cap = np.linspace(0, radius, 15)
    theta_cap = np.linspace(0, 2 * np.pi, 50)
    r_cap_grid, theta_cap_grid = np.meshgrid(r_cap, theta_cap)
    x_cap = r_cap_grid * np.cos(theta_cap_grid)
    y_cap = r_cap_grid * np.sin(theta_cap_grid)

    # Top cap
    z_top = np.full_like(x_cap, half_length)
    ax.plot_surface(x_cap, y_cap, z_top, alpha=cylinder_alpha,
                    color=cylinder_color, edgecolor='none')

    # Bottom cap
    z_bottom = np.full_like(x_cap, -half_length)
    ax.plot_surface(x_cap, y_cap, z_bottom, alpha=cylinder_alpha,
                    color=cylinder_color, edgecolor='none')

    # Draw cylinder wireframe edges for clarity
    theta_wire = np.linspace(0, 2 * np.pi, 50)
    ax.plot(radius * np.cos(theta_wire), radius * np.sin(theta_wire),
            np.full_like(theta_wire, half_length), 'k-', alpha=0.3, linewidth=0.5)
    ax.plot(radius * np.cos(theta_wire), radius * np.sin(theta_wire),
            np.full_like(theta_wire, -half_length), 'k-', alpha=0.3, linewidth=0.5)

    # --- Draw each ray ---
    for i in range(N):
        color = ray_colors[i]
        origin = origins[i]
        direction = directions[i]
        entry = entry_points[i]
        exit_pt = exit_points[i]

        # Check if this ray hit the cylinder
        hit = not np.any(np.isnan(entry))

        # Draw origin point
        ax.scatter(*origin, color=color, s=100, marker='o',
                   edgecolors='black', linewidths=1, label=f'Ray {i+1} origin' if i < 5 else None)

        if hit:
            # Draw entry and exit points
            ax.scatter(*entry, color=color, s=50, marker='o',
                       edgecolors='black', linewidths=1.5, zorder=5)
            ax.scatter(*exit_pt, color=color, s=50, marker='o',
                       edgecolors='black', linewidths=1.5, zorder=5)

            # Draw solid line through cylinder (entry to exit)
            ax.plot([entry[0], exit_pt[0]],
                    [entry[1], exit_pt[1]],
                    [entry[2], exit_pt[2]],
                    color=color, linewidth=2.5, solid_capstyle='round')

            # Calculate how far to extend the dashed line
            dir_norm = direction / np.linalg.norm(direction)

            # Dashed line from origin to entry
            ax.plot([origin[0], entry[0]],
                    [origin[1], entry[1]],
                    [origin[2], entry[2]],
                    color=color, linewidth=1.5, linestyle='--', alpha=0.7)

            # Dashed line extending beyond exit
            extend_dist = extend_factor * np.linalg.norm(exit_pt - entry)
            end_point = exit_pt + dir_norm * extend_dist
            ax.plot([exit_pt[0], end_point[0]],
                    [exit_pt[1], end_point[1]],
                    [exit_pt[2], end_point[2]],
                    color=color, linewidth=1.5, linestyle='--', alpha=0.7)
        else:
            # Ray missed - draw extended dashed line
            dir_norm = direction / np.linalg.norm(direction)
            extend_dist = 4 * half_length  # Extend far enough to show it misses
            end_point = origin + dir_norm * extend_dist
            ax.plot([origin[0], end_point[0]],
                    [origin[1], end_point[1]],
                    [origin[2], end_point[2]],
                    color=color, linewidth=1.5, linestyle=':', alpha=0.5)

    # --- Set labels and adjust view ---
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    # Set equal aspect ratio
    max_range = max(radius, half_length) * 1.5

    # Include origins in the range calculation
    all_points = [origins]
    if not np.all(np.isnan(entry_points)):
        all_points.append(entry_points[~np.isnan(entry_points).any(axis=1)])
    if not np.all(np.isnan(exit_points)):
        all_points.append(exit_points[~np.isnan(exit_points).any(axis=1)])

    all_points = np.vstack(all_points)

    x_range = [min(all_points[:, 0].min(), -radius) - 2,
               max(all_points[:, 0].max(), radius) + 2]
    y_range = [min(all_points[:, 1].min(), -radius) - 2,
               max(all_points[:, 1].max(), radius) + 2]
    z_range = [min(all_points[:, 2].min(), -half_length) - 2,
               max(all_points[:, 2].max(), half_length) + 2]

    ax.set_xlim(x_range)
    ax.set_ylim(y_range)
    ax.set_zlim(z_range)

    # Try to set equal aspect ratio (works in newer matplotlib)
    try:
        ax.set_aspect('equal')
    except NotImplementedError:
        pass

    ax.set_title('Muon Trajectories through CMS Detector')

    return ax

################################################################################

################################################################################

def demo_visualization():
    """Run a demonstration of the ray-cylinder intersection and visualization."""

    # CMS-like dimensions (in meters)
    cms_radius = 7.0
    cms_half_length = 10.5

    # Create test muons coming from underground (negative y)
    np.random.seed(42)

    origins = np.array([
        [0, -20, 0],           # Directly below center
        [3, -18, 5],           # Off-center
        [-2, -22, -3],         # Another off-center
        [8, -15, 0],           # Near edge
        [20, -20, 0],          # Will miss (too far in x)
        [0, -25, 15],          # Will miss (too far in z)
    ])

    # Direction vectors (roughly pointing upward with some spread)
    directions = np.array([
        [0, 1, 0],             # Straight up
        [-0.1, 1, -0.2],       # Slight angle
        [0.15, 1, 0.1],        # Slight angle
        [-0.3, 1, 0],          # Angled toward center
        [0, 1, 0],             # Straight up (but will miss)
        [0, 1, -0.3],          # Angled (but will miss)
    ])

    # Compute intersections
    entry_points, exit_points = ray_cylinder_intersection(
        origins, directions, cms_radius, cms_half_length
    )

    print("Entry points:")
    print(entry_points)
    print("\nExit points:")
    print(exit_points)

    # Visualize
    ax = visualize_ray_cylinder_intersections(
        origins, directions, entry_points, exit_points,
        cms_radius, cms_half_length
    )
    ax.view_init(vertical_axis='y')
    ax.set_xlim(-20,20)
    ax.set_ylim(-20,20)
    ax.set_zlim(-20,20)

    plt.tight_layout()
    plt.show()

    return ax

################################################################################
