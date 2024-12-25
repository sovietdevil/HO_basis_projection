import numpy as np
from matplotlib import pyplot as plt
import sys
from structure_generator import (
    construction_periodic,
    set_dopants,
    construct_amorphous)
from GS_waves import (
    defocus_stack,
    select_defocus_profile,
    select_freq_range)
from ase import Atoms
from abtem.waves import PlaneWave
import circle_fit as cf
from tqdm import tqdm

def sketch_defocus_circle(ew, start, end, step, energy, sampling, width=1, ecllipse_corr=False):
    def_stack = defocus_stack(ew, energy, sampling, start, end, step, progress=False)
    n_row, n_col = ew.shape
    if width > 0:
        profile = select_defocus_profile(def_stack, n_row//2, n_col//2, width=width)
    else:
        profile = def_stack[:, n_row//2, n_col//2]
    circle_fit = np.vstack((np.real(profile), np.imag(profile))).T
    if ecllipse_corr:
        pass
    else:
        xc, yc, r, sigma = cf.least_squares_circle(circle_fit)
    return profile, xc, yc, r, sigma

def calculate_curvature(coordinates):
    x = coordinates[:,0]; y = coordinates[:,1]
    dx = np.gradient(x)
    dy = np.gradient(y)
    ddx = np.gradient(dx)
    ddy = np.gradient(dy)
    curvature = np.abs(ddx * dy - dx * ddy) / (dx**2 + dy**2)**1.5
    return curvature

def curvature_based_selection(defocus_profile, threshold=2):
    return 0

def plot_circle(xc, yc, r, color, linewidth=1):
    theta = np.linspace(0, 2*np.pi, 100)
    plt.plot(xc+r*np.cos(theta), yc+r*np.sin(theta), color=color, linewidth=linewidth)