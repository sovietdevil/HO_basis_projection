import torch
import numpy as np
from torch.optim import LBFGS, Adam

initial_parameters = torch.tensor([0.6593729, -0.04567545,
                                   0.7216003, 0.40138984,
                                   -0.8390287, -0.56562936,
                                   0.35573682, 0.10767087,
                                   0.06771779, 0.10599555])

class Gaussian_Fitting(torch.nn.Module):
    def __init__(self, n_grids, sampling, range_constrain=False, range_list=None, periodic_boundary=True):
        super().__init__()
        self.n_grids = n_grids
        self.n_row = int(np.sqrt(n_grids))
        self.sampling = sampling
        self.length = sampling * self.n_row
        self.coeffs = torch.nn.parameter.Parameter(torch.rand(10), requires_grad=True)
        if range_list==None:
            range_list = torch.tensor([[-0.5, 0.5], [0, 1], [0, 1], [-2, 2],
                                   [-2, 2], [-2, 2], [-np.pi, np.pi], [-1, 1], [-1, 1]])
        self.range_constrain = range_constrain
        self.range_list = range_list
        self.periodic_boundary = periodic_boundary
    def forward(self, input):
        total = input
        a = self.length
        x = np.linspace(-a/2, a/2, self.n_row)
        sampling = a/self.n_row
        X, Y = np.meshgrid(x, x)
        X = torch.tensor(X); Y = torch.tensor(Y)
        if self.range_constrain:
            coeffs = torch.clamp(self.coeffs, min=self.range_list[:,0], max=self.range_list[:,1])
        else:
            coeffs = self.coeffs
        bg_re, bg_im, amp, var, kx, ky, k2, phi0, rx, ry = coeffs
        if self.periodic_boundary:
            Xp = ((X - rx) + a/2)%a - a/2
            Yp = ((Y - ry) + a/2)%a - a/2
            r2 = Xp**2 + Yp**2
            phase_shift = k2*r2 + kx*Xp + ky*Yp + phi0
            result_re = bg_re + amp * torch.exp(-r2/(2*var)) * torch.cos(phase_shift)
            result_im = bg_im + amp * torch.exp(-r2/(2*var)) * torch.sin(phase_shift)
        else:
            r2 = (X - rx)**2 + (Y - ry)**2
            phase_shift = k2*X**2 + kx*X + k2*Y**2 + ky*Y + phi0
            result_re = bg_re + amp * torch.exp(-r2/(2*var)) * torch.cos(phase_shift)
            result_im = bg_im + amp * torch.exp(-r2/(2*var)) * torch.sin(phase_shift)
        result = (result_re + 1j*result_im).ravel()
        total = total.clone() + result
        return total
    
def optimize_columns(ew, sampling, n_iter=100, range_constrain=False, range_list=None, periodic_boundary=True):
    n_grids = ew.shape[0]*ew.shape[1]
    model = Gaussian_Fitting(n_grids, sampling, range_constrain, range_list, periodic_boundary)
    state_dict = model.state_dict()
    state_dict["coeffs"] = initial_parameters
    model.load_state_dict(state_dict)
    x = torch.tensor(np.zeros((n_grids,)), dtype=torch.complex32)
    y = model(x).detach().numpy().reshape(ew.shape)
    Y = torch.tensor(ew.ravel())
    optimizer = LBFGS(model.parameters(), lr=0.001, line_search_fn="strong_wolfe")
    def closure():
        optimizer.zero_grad()
        loss = torch.norm(Y - model(x))
        loss.backward()
        return loss
    model.train()
    for step in range(n_iter):
        loss = optimizer.step(closure)
    return model.state_dict()["coeffs"].detach().numpy(), loss