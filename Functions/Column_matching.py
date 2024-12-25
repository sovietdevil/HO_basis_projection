import numpy as np
import scipy
import hyperspy.api as hs
import atomap.api as am
from scipy.special import j0, factorial, hermite
from scipy.optimize import minimize
from abtem.core.energy import energy2wavelength

def match_peaks(data, refine=False, pca=False, threshold=None, separation=5, remove_bounds=True, bounds_sep=None, fix_negative_values=True):
    data = (data - np.mean(data))/np.std(data)
    signal = hs.signals.Signal2D(data)
    peak_sites = am.get_atom_positions(signal, pca=pca, separation=separation)
    values = []
    for index in range(peak_sites.shape[0]):
        i_row = int(peak_sites[index,1])
        i_col = int(peak_sites[index,0])
        value = data[i_row, i_col]
        values.append(value)
    values = np.array(values)
    if threshold is not None:
        peak_sites = peak_sites[values>=threshold]
    sublattice = am.Sublattice(peak_sites, image=data, fix_negative_values=fix_negative_values)
    sublattice.construct_zone_axes()
    if refine:
        sublattice.find_nearest_neighbors()
        sublattice.refine_atom_positions_using_center_of_mass()
        sublattice.refine_atom_positions_using_2d_gaussian()
        peak_sites = sublattice.atom_positions
    if remove_bounds:
        if bounds_sep is None:
            bounds_sep = separation
        m, n = data.shape
        condition_x = np.logical_and(peak_sites[:,0]-bounds_sep>1, peak_sites[:,0]+bounds_sep<m)
        condition_y = np.logical_and(peak_sites[:,1]-bounds_sep>1, peak_sites[:,1]+bounds_sep<n)
        condition = np.logical_and(condition_x, condition_y)
        filtered_pos = np.array([[peak_sites[i, 0], peak_sites[i, 1]] for i in range(len(peak_sites[:,0])) if condition[i]])
        peak_sites = filtered_pos
        lattice = np.array(sublattice.zones_axis_average_distances)
    return peak_sites, lattice

def select_column(ew, peak_sites, distance, index):
    n_col, n_row = peak_sites[index,:]
    dist = int(distance)
    lx, ux = n_col - dist, n_col + dist
    ly, uy = n_row - dist, n_row + dist
    ew_sel = ew[ly:uy, lx:ux]
    return ew_sel

def HO_transform(results, krange, param, n_basis=2):
    Basis = np.zeros(n_basis,).astype(np.complex64)
    Proj = np.zeros(results.shape).astype(np.complex64)
    length = len(results)
    assert length == len(krange)
    k_sampling = krange[1] - krange[0]
    for n in range(n_basis):
        poly = np.exp(-(krange/param)**2/2)\
            /np.sqrt(2**n*factorial(n)*np.sqrt(np.pi)*param)\
            *hermite(n)(krange/param)
        y = poly * results
        result = np.sum(y)*k_sampling
        Basis[n] = result
        Proj = Proj + result * poly
    Proj = Proj * (krange[-1]-krange[0])/length
    return Basis, Proj

def Bessel_transform(func, sampling, kmin, kmax, length, x0=0, y0=0):
    results = []
    component = np.ones(func.shape)*(0+0j)
    krange = np.linspace(kmin, kmax, length)*np.pi*2
    m, n = func.shape; Area = m * n
    x = np.linspace(-sampling*m/2, sampling*m/2, m) - x0
    y = np.linspace(-sampling*n/2, sampling*n/2, n) - y0
    X, Y = np.meshgrid(x, y)
    R = np.sqrt(X**2 + Y**2)
    for k in krange:
        y = np.sqrt(k) * j0(k * R) * func
        result = np.sum(y)*sampling**2
        results.append(result)
        component += result * np.sqrt(k) * j0(k * R)
    component = np.array(component) * (kmax-kmin) / length
    results = np.array(results)/np.sqrt(krange)
    return results, component

def inv_Bessel(results, m, n, sampling, kmin, kmax, length, x0=0, y0=0):
    func = np.ones((m, n))*(0+0j)
    krange = np.linspace(kmin, kmax, length)*np.pi*2
    x = np.linspace(-sampling*m/2, sampling*m/2, m) - x0
    y = np.linspace(-sampling*n/2, sampling*n/2, n) - y0
    X, Y = np.meshgrid(x, y)
    R = np.sqrt(X**2 + Y**2)
    for n, k in enumerate(krange):
        func += results[n]*np.sqrt(k) * np.sqrt(k) * j0(k * R)
    func = func * (kmax - kmin) / length
    return func

def propagation_Bessel(results, distance, kmin, kmax, length, sampling, energy, shape, x0=0, y0=0):
    wavelength = energy2wavelength(energy)
    krange = np.linspace(kmin, kmax, length)
    results_shift = results * np.exp(-1j * np.pi * wavelength * distance * krange**2)
    m, n = shape
    wavefunction = inv_Bessel(results_shift, m, n, sampling, kmin, kmax, length, x0, y0)
    return results_shift, wavefunction

def match_parameters_Bessel(ew, sampling, gmin, gmax, length=50, x0=[0, 0, 1, 0], method="BFGS"):
    def Bessel_fit(coeffs):
        rx, ry, bg_amp, bg_pha = coeffs
        bg = bg_amp*np.exp(1j*bg_pha)
        results, components = Bessel_transform(ew-bg, sampling, gmin, gmax, length, rx, ry)
        return np.linalg.norm(ew-bg-components)
    res = minimize(Bessel_fit, x0, method=method)
    return res.x
