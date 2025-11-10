import numpy as np
import pandas as pd
import matplotlib.pylab as plt

import seaborn as sns

from skspatial.objects import Line, Plane
from skspatial.plotting import plot_3d


from skspatial.objects import Line, Cylinder, Point, Points
from skspatial.plotting import plot_3d

import phasespace

import tensorflow

import bisect
import numpy as np
import matplotlib.pylab as plt
import pandas as pd

import seaborn as sns

import numpy as np
from sklearn.mixture import GaussianMixture
from scipy.stats import multivariate_normal

import numpy as np
from scipy.interpolate import griddata
from scipy.integrate import quad, trapezoid
from scipy.interpolate import CubicSpline

import matplotlib.pylab as plt
from scipy import stats
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import matplotlib

from scipy.interpolate import LinearNDInterpolator
import scipy

from scipy import integrate 

from scipy.stats import gaussian_kde
from sklearn.neighbors import KernelDensity

import numpy as np
from typing import Dict, Any, List




import dm_generation_tools as dgt
import detector_simulation_tools as dst


import time

####################################
import warnings
# Suppress all warnings
warnings.filterwarnings("ignore")


import pickle

##########################################################################
def generate_data_from_spline(spl, npts):

    #print("Generating random points..........")

    min_val,max_val = spl.x[0],spl.x[-1]

    start = time.time()
    #npts = 5000
    values = []
    values_random_numbers = np.random.random(npts)

    for i in range(npts):

        solutions = spl.solve(values_random_numbers[i])
        filter = (solutions>=min_val) & (solutions<=max_val)
        solutions = solutions[filter]

        if len(solutions)==1:
            #ynew=ynew[0]
            #print(ynew)
            values.append(solutions[0])

    return values
##########################################################################


class RaggedSplineIndexLookup:
    """
    Vectorized nearest (E, d) lookup that returns **indices**:

      - energy_idx: index into self.energy_keys
      - dist_idx:   index into self.dist_keys[energy_idx]

    If the requested distance is GREATER than the largest distance
    for that energy, we return dist_idx = -1 for that entry.
    """

    def __init__(self, data: Dict[float, Dict[float, Any]]):
        energies: List[float] = []
        dist_keys: List[np.ndarray] = []
        dist_maps: List[Dict[float, int]] = []

        # normalize and store
        for E in sorted(data.keys(), key=float):
            E_f = float(E)
            inner = data[E]
            d_sorted = sorted(float(x) for x in inner.keys())
            energies.append(E_f)
            dist_keys.append(np.array(d_sorted, dtype=float))
            dist_maps.append({d: i for i, d in enumerate(d_sorted)})

        self.energy_keys = np.array(energies, dtype=float)  # (nE,)
        self.dist_keys = dist_keys                          # list of np.ndarray
        self.dist_maps = dist_maps                          # not strictly needed here
        self._has_data = self.energy_keys.size > 0

    # --------------------------------------------------
    # 1. snap to ideal piecewise grid (vectorized)
    # --------------------------------------------------
    @staticmethod
    def snap_energies(E: np.ndarray) -> np.ndarray:
        E = np.asarray(E, dtype=float)
        out = np.empty_like(E)

        m0 = E < 100
        m1 = (E >= 100) & (E < 1000)
        m2 = (E >= 1000) & (E < 10_000)
        m3 = (E >= 10_000) & (E < 100_000)
        m4 = E >= 100_000

        out[m0] = 100.0

        if np.any(m1):
            tmp = E[m1]
            base = 100.0
            step = 10.0
            idx = np.rint((tmp - base) / step).astype(int)
            val = base + idx * step
            val = np.minimum(val, 1000.0)
            out[m1] = val

        if np.any(m2):
            tmp = E[m2]
            base = 1000.0
            step = 100.0
            idx = np.rint((tmp - base) / step).astype(int)
            val = base + idx * step
            val = np.minimum(val, 10_000.0)
            out[m2] = val

        if np.any(m3):
            tmp = E[m3]
            base = 10_000.0
            step = 1000.0
            idx = np.rint((tmp - base) / step).astype(int)
            val = base + idx * step
            val = np.minimum(val, 100_000.0)
            out[m3] = val

        out[m4] = 100_000.0

        return out

    # --------------------------------------------------
    # 2. map snapped energy → closest actual energy INDEX
    # --------------------------------------------------
    def map_to_actual_energy_index(self, E_snap: float) -> int:
        if not self._has_data:
            return -1

        ek = self.energy_keys
        pos = np.searchsorted(ek, E_snap)
        if pos == 0:
            return 0
        if pos == ek.size:
            return ek.size - 1

        before = ek[pos - 1]
        after = ek[pos]
        if abs(before - E_snap) <= abs(after - E_snap):
            return pos - 1
        else:
            return pos

    # --------------------------------------------------
    # 3. nearest distance index for ONE energy index
    #    with "cap at max" -> return -1
    # --------------------------------------------------
    def _nearest_distance_index_for_energy_idx(self, e_idx: int, d_query: np.ndarray) -> np.ndarray:
        """
        For this energy (by index), find nearest distance index for each query.
        BUT: if d_query > max_available_distance → return -1 for that element.
        """
        real_ds = self.dist_keys[e_idx]   # shape (n_real,)
        q = np.asarray(d_query, dtype=float)

        max_d = real_ds[-1]  # largest distance we actually have for this energy

        # start with all -1 (meaning: no match)
        out_idx = np.full(q.shape, -1, dtype=int)

        # mask for those that are within range
        m_valid = q <= max_d
        if not np.any(m_valid):
            return out_idx  # all out of range

        qv = q[m_valid]

        idx = np.searchsorted(real_ds, qv, side="left")

        n_real = real_ds.shape[0]
        idx_left = np.clip(idx - 1, 0, n_real - 1)
        idx_right = np.clip(idx, 0, n_real - 1)

        left_vals = real_ds[idx_left]
        right_vals = real_ds[idx_right]

        left_diff = np.abs(qv - left_vals)
        right_diff = np.abs(qv - right_vals)

        choose_right = right_diff < left_diff
        chosen_idx = np.where(choose_right, idx_right, idx_left)

        # fill only valid positions
        out_idx[m_valid] = chosen_idx
        # positions where q > max_d stay -1

        return out_idx

    # --------------------------------------------------
    # 4. public vectorized interface
    # --------------------------------------------------
    def get_indices(self, E_arr, d_arr):
        """
        returns (energy_idx_array, dist_idx_array)
        dist_idx_array entries are -1 if:
           - energy not found (closest)
           - OR distance > max available for that energy
        """
        E_arr = np.asarray(E_arr, dtype=float)
        d_arr = np.asarray(d_arr, dtype=float)
        assert E_arr.shape == d_arr.shape

        N = E_arr.size
        E_snap = self.snap_energies(E_arr)

        energy_idx_out = np.full(N, -1, dtype=int)
        dist_idx_out = np.full(N, -1, dtype=int)

        unique_Es, inverse = np.unique(E_snap, return_inverse=True)

        for i, E_snap_val in enumerate(unique_Es):
            mask = (inverse == i)
            idxs = np.nonzero(mask)[0]
            d_sub = d_arr[mask]

            e_idx = self.map_to_actual_energy_index(E_snap_val)
            if e_idx == -1:
                # no energies at all
                continue

            d_idx_sub = self._nearest_distance_index_for_energy_idx(e_idx, d_sub)

            energy_idx_out[idxs] = e_idx
            dist_idx_out[idxs] = d_idx_sub

        return energy_idx_out, dist_idx_out

##########################################################################
# Does it in chunks
def count_pairs(E_idx: np.ndarray, D_idx: np.ndarray):
    """
    Given two 1D arrays of the same length, return
    - pairs: (k, 2) array of unique (E_idx, D_idx)
    - counts: (k,) array of counts for each pair
    """
    E_idx = np.asarray(E_idx)
    D_idx = np.asarray(D_idx)
    assert E_idx.shape == D_idx.shape

    # stack into (N, 2)
    pairs = np.stack([E_idx, D_idx], axis=1)  # shape (N, 2)

    # view as a single dtype so np.unique can work along axis=0
    # make a structured view
    dtype = np.dtype([('e', E_idx.dtype), ('d', D_idx.dtype)])
    pairs_view = pairs.view(dtype)

    uniq, counts = np.unique(pairs_view, return_counts=True)

    # convert uniq back to a regular (k, 2) array
    uniq_arr = np.column_stack([uniq['e'], uniq['d']])
    return uniq_arr, counts

##########################################################################
