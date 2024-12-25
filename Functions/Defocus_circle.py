from abtem.waves import PlaneWave
import numpy as np
import scipy
from abtem.core.energy import energy2wavelength
import hyperspy.api as hs
from tqdm import tqdm
import circle_fit as cf

def generate_continious_phase(complex_series):
    output_phase = np.angle(complex_series)
    for n in range(len(complex_series)-1):
        if (output_phase[n]-output_phase[n+1])>np.pi:
            output_phase[n+1:] = output_phase[n+1:] + 2*np.pi
        elif (output_phase[n+1]-output_phase[n])>np.pi:
            output_phase[n+1:] = output_phase[n+1:] - 2*np.pi
    return output_phase

def find_angle(points, xc, yc):
    rx = points[:,0] - xc
    ry = points[:,1] - yc
    r_complex = rx + 1j*ry
    return generate_continious_phase(r_complex)