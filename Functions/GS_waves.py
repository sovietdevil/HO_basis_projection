from abtem.waves import PlaneWave
import numpy as np
import scipy
from abtem.core.energy import energy2wavelength
import hyperspy.api as hs
from tqdm import tqdm

def select_freq_range(exitwave, gmin, gmax, sampling):
    exitwave = np.array(exitwave)
    m, n = exitwave.shape
    ft_exitwave = scipy.fft.fft2(exitwave)
    freq_gx = np.fft.fftfreq(m, sampling)
    freq_gy = np.fft.fftfreq(n, sampling)
    gx, gy = np.meshgrid(freq_gx, freq_gy)
    g2 = gx ** 2 + gy ** 2
    ft_exitwave[g2 < gmin ** 2] = 0
    ft_exitwave[g2 > gmax ** 2] = 0
    ew = scipy.fft.ifft2(ft_exitwave)
    return ew

def multislice_from_structure(structure, sampling, energy, gmin=0, gmax=2):
    '''
    Conduct multislice simulation.
    '''
    pw = PlaneWave(sampling=sampling, energy=energy)
    ew = pw.multislice(structure)
    return select_freq_range(np.array(ew.array), gmin, gmax, sampling)

def from_files(fname, fpath='./', reverse=False):
    '''
    Extract the exit wave function from .tif files.
    '''
    #Check for the path name
    if fpath[-1] == '/':
        fullname = fpath+fname
    else:
        fullname = fpath+'/'+fname
    #Extract data from the file
    loaded = np.array(hs.load(fullname))
    real = loaded[:,:,0]
    imag = loaded[:,:,1]
    if reverse==True:
        wavefunction = real - 1j*imag
    else:
        wavefunction = real + 1j*imag
    return wavefunction

def propagation_ew(waves, distance, sampling, energy):
    wavelength = energy2wavelength(energy)
    waves = np.array(waves)
    m, n = waves.shape
    kx = np.fft.fftfreq(m, sampling)
    ky = np.fft.fftfreq(n, sampling)
    Kx, Ky = np.meshgrid(kx, ky)
    k2 = Kx ** 2 + Ky ** 2
    kernel = np.exp(- 1.j * k2 * np.pi * wavelength * distance)
    waves = scipy.fft.ifft2(scipy.fft.fft2(waves)*kernel)
    return waves

def real_space_filter(wavefunction, sampling, center=None, coeff=1, sigma=1):
    wavefunction = np.array(wavefunction)
    m, n = wavefunction.shape
    a = m*sampling/2
    if center is None:
        center = (m//2, n//2)
    x = np.linspace((1-center[0])*sampling, (m-center[0])*sampling, m)
    y = np.linspace((1-center[1])*sampling, (n-center[1])*sampling, n)
    x, y = np.meshgrid(x, y)
    r2 = x**2 + y**2
    mask = coeff * np.exp(-r2/(2*sigma))
    wave_filtered = wavefunction * mask
    wave_remained = wavefunction * (1-mask)
    return wave_filtered, wave_remained    

def Fourier_space_filter(wavefunction, sampling, coeff=1, sigma=1):
    wavefunction = np.array(wavefunction)
    m, n = wavefunction.shape
    kx = np.fft.fftfreq(m, sampling)
    ky = np.fft.fftfreq(n, sampling)
    interval = 1/(sampling*m)
    Kx, Ky = np.meshgrid(kx, ky)
    k2 = Kx ** 2 + Ky ** 2
    mask = coeff * np.exp(-k2/(2*sigma))
    wavefunc_filtered = scipy.fft.ifft2(scipy.fft.fft2(wavefunction)*mask)
    wavefunc_remained = scipy.fft.ifft2(scipy.fft.fft2(wavefunction)*(coeff-mask))
    return wavefunc_filtered, wavefunc_remained

def sel_central(ew, n_period, shift=[0, 0, 0, 0]):
    n_row, n_col = ew.shape
    upper = int((1/2-1/(4*n_period))*n_col + shift[0])
    lower = int((1/2+1/(4*n_period))*n_col + shift[1])
    left = int((1/2-1/(4*n_period))*n_row + shift[2])
    right = int((1/2+1/(4*n_period))*n_row + shift[3])
    ew_sel = ew[upper:lower, left:right]
    return ew_sel

def defocus_stack(wavefunction, energy, sampling, start, end, step, progress=True):
    defocus = np.arange(start, end, step)
    stack = []
    if progress:
        for deltaf in tqdm(defocus):
            stack.append(propagation_ew(wavefunction, deltaf, sampling, energy))
    else:
        for deltaf in defocus:
            stack.append(propagation_ew(wavefunction, deltaf, sampling, energy))
    return np.array(stack)

def select_defocus_profile(ew_stack, n_row, n_col, width=2):
    lx, ux = n_col - width, n_col + width
    ly, uy = n_row - width, n_row + width
    stack_sel = ew_stack[:,ly:uy, lx:ux]
    return np.sum(np.sum(stack_sel,axis=1),axis=1)/(2*width)**2

def fit_circle(circle_profile, background_profile):
    xc, yc, r, sigma = cf.least_squares_circle(circle_profile - background_profile)
    return xc, yc, r, sigma