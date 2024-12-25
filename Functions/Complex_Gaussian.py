import numpy as np
from sklearn.linear_model import LinearRegression
from lmfit import Model
from scipy.optimize import minimize

def __init__():
    super().__init__()

def Gaussian_function(k, a_re, a_im, b_re, b_im, c_re, c_im):
    theta = -b_im * k**2
    y_re = c_re + np.exp(-b_re * k**2) * (a_re*np.cos(theta) - a_im*np.sin(theta))
    y_im = c_im + np.exp(-b_re * k**2) * (a_re*np.sin(theta) + a_im*np.cos(theta))
    return y_re + 1j*y_im

def Gaussian_in_focus(k, a_re, a_im, b, c_re, c_im):
    y_re = c_re + a_re * np.exp(-b * k**2)
    y_im = c_im + a_im * np.exp(-b * k**2)
    return y_re + 1j*y_im

def Gaussian_fitting(profile, freq_data, freq_range, defocus=True):
    index_sel = np.where(np.logical_and(freq_data>=freq_range[0],
                                        freq_data<=freq_range[1]))
    #Step 1: fit line for the complex component
    fit_x = freq_data[index_sel]; fit_y = profile[index_sel]
    peak_complex = fit_y[0]
    re = np.real(fit_y).ravel().reshape(-1, 1)
    im = np.imag(fit_y).ravel()
    res = LinearRegression().fit(re, im)
    slope = res.coef_[0]; intercept = res.intercept_
    #Step 2: fit Gaussian function for the line
    X = (fit_x**2).reshape(-1, 1)
    y = np.log(np.abs(fit_y-intercept))
    res = LinearRegression().fit(X, y); k_init = -res.coef_[0]
    if defocus:
        model = Model(Gaussian_function)
        params = model.make_params(a_re=np.real(peak_complex),
                                a_im=np.imag(peak_complex)-intercept,
                                b_re=intercept, b_im=0,
                                c_re=k_init, c_im=0)
        result = model.fit(fit_y, params, k=fit_x)
    else:
        model = Model(Gaussian_in_focus)
        params = model.make_params(a_re=np.real(peak_complex),
                                a_im=np.imag(peak_complex)-intercept,
                                b=0,
                                c_re=k_init, c_im=0)
        result = model.fit(fit_y, params, k=fit_x)
    return result

def Gaussian_function_partition(coeffs, krange):
    k_thres = coeffs[0]
    factor = coeffs[1]
    coeffs = coeffs[2:].reshape(2, -1)
    a_re, a_im, b_re, b_im, c_re, c_im = coeffs[0,:]
    theta = -b_im * krange**2
    y1_re = c_re + np.exp(-b_re * krange**2) * (a_re*np.cos(theta) - a_im*np.sin(theta))
    y1_im = c_im + np.exp(-b_re * krange**2) * (a_re*np.sin(theta) + a_im*np.cos(theta))
    a_re, a_im, b_re, b_im, c_re, c_im = coeffs[1,:]
    y1 = y1_re + 1j*y1_im
    theta = -b_im * krange**2
    y2_re = c_re + np.exp(-b_re * krange**2) * (a_re*np.cos(theta) - a_im*np.sin(theta))
    y2_im = c_im + np.exp(-b_re * krange**2) * (a_re*np.sin(theta) + a_im*np.cos(theta))
    y2 = y2_re + 1j*y2_im
    #weight = 1/(1+np.exp(-factor*(k_thres**2-krange**2)))
    weight = 0.5 * (1 + np.tanh(-factor*(k_thres**2-krange**2)))
    return y2*weight + y1*(1-weight)

def Modified_Gaussian_fitting(results, krange):
    gmin = krange[0]; gmax = krange[-1]; length=len(krange)
    #Step 1: Initial fit
    fitted = Gaussian_fitting(results, krange, [gmin+0.01, 1], defocus=True)
    coeff = np.array(list(fitted.best_values.values()))
    bg = coeff[-2] + 1j*coeff[-1]
    recover = Gaussian_function(krange, *coeff)
    residue = results - recover + bg
    fit_res = Gaussian_fitting(residue, krange, [1, gmax], defocus=True)
    coeff_res = np.array(list(fit_res.best_values.values()))
    bg_res = coeff_res[-2] + 1j*coeff_res[-1]
    recover_res = Gaussian_function(krange, *coeff_res)
    def partitioned_Gaussian_fit(coeffs):
        k_thres = coeffs[0]
        factor = coeffs[1]
        coeffs = coeffs[2:].reshape(2, -1)
        a_re, a_im, b_re, b_im, c_re, c_im = coeffs[0,:]
        theta = -b_im * krange**2
        y1_re = c_re + np.exp(-b_re * krange**2) * (a_re*np.cos(theta) - a_im*np.sin(theta))
        y1_im = c_im + np.exp(-b_re * krange**2) * (a_re*np.sin(theta) + a_im*np.cos(theta))
        a_re, a_im, b_re, b_im, c_re, c_im = coeffs[1,:]
        y1 = y1_re + 1j*y1_im
        theta = -b_im * krange**2
        y2_re = c_re + np.exp(-b_re * krange**2) * (a_re*np.cos(theta) - a_im*np.sin(theta))
        y2_im = c_im + np.exp(-b_re * krange**2) * (a_re*np.sin(theta) + a_im*np.cos(theta))
        y2 = y2_re + 1j*y2_im
        #weight = 1/(1+np.exp(-factor*(k_thres**2-krange**2)))
        weight = 0.5 * (1 + np.tanh(-factor*(k_thres**2-krange**2)))
        recover = y2*weight + y1*(1-weight)
        return np.linalg.norm(np.abs(recover-results))
    coeffs=np.hstack((1, 50, coeff, coeff_res)).ravel()
    res = minimize(partitioned_Gaussian_fit, coeffs, method="BFGS")
    return res.x