import numpy as np
import pandas as pd
import matplotlib.pylab as plt

from skspatial.objects import Line, Plane
from skspatial.plotting import plot_3d

from skspatial.objects import Line, Cylinder, Point, Points
from skspatial.plotting import plot_3d

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

def draw_detector():
    # Define detector
    #origin_detector = [0, 0, 0]

    nmuons = 0

    # detector, units are meters. x is direction of beam and z is up
    #cylinder = Cylinder.from_points([-10.5, 0, 0], [10.5, 0, 0], 7.5)
    # Mock detectors
    cylinder = Cylinder.from_points([-15, 0, 0], [15, 0, 0], 7.5)

    #if MAKE_PLOTS:
    #fig1 = plt.figure(figsize=(6,6))
    #ax1 = fig1.add_subplot(1,1,1,projection='3d')

    fig2 = plt.figure(figsize=(12,6))
    ax2 = fig2.add_subplot(1,1,1,projection='3d')

    #fig3 = plt.figure(figsize=(4,4))
    #ax3 = fig3.add_subplot(1,1,1)

    # Draw detector
    cylinder.plot_3d(ax2, alpha=0.2)
    ax2.set_xlim(-100,100)
    ax2.set_ylim(-100,100)
    ax2.set_zlim(-100,20)

    return fig2, ax2, cylinder

##########################################################################
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
        pb.plot_3d(ax, s=75, c=color, alpha=0.2)

    if point_b[0]==point_b[0] and point_a[0]==point_a[0]:
        #line = Line(point_a, point_b)
        line = Line(pa, pb-pa)

        line.plot_3d(ax, c=color)

    return pa,pb,line

################################################################################

