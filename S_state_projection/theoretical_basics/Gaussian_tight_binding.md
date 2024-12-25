#### Functions
The General solutions of the Schrödinger's equation:
$$
\left(\frac{\hat{p}^2}{2m}+\frac{1}{2}m\omega^2x^2\right)\psi_n(x)=E_n\psi_n(x)
$$
are:
$$
\psi_n(x)=\frac{1}{\sqrt{2^nn!}}\left(\frac{m\omega}{\pi\hbar}\right)^{1/4}e^{-\frac{m\omega x^2}{2\hbar}}H_n\left(\sqrt{\frac{m\omega}{\hbar}}x\right)
$$
with $H_n(x)$ the Hermite polynomial, generated from
$$
H_n(x)=(-1)^ne^{x^2}\frac{d^n}{dx^n}\left(e^{-x^2}\right)
$$
The corresponding eigen-energies are:
$$
E_n=\hbar\omega(n+1/2)
$$
Convert using $\sigma=\frac{\hbar}{m\omega}$ derives
$$
\psi_n(x)=\frac{1}{\sqrt{2^n\sigma n!\sqrt{\pi}}}H_n\left(\frac{x}{\sigma}\right)e^{-\frac{x^2}{2\sigma^2}}
$$
#### Estimations
Based on the Gaussian potential distribution
$$
\exp\left(-\frac{r^2}{2\sigma^2}\right)=1-\frac{r^2}{2\sigma^2}+\frac{r^4}{8\sigma^4}+...
$$
The first two term is $1-r^2/2\sigma^2$, which implies the quadratic approximation of the potential distribution for coordinates close to the column peak.
#### Gaussian Periodic Potential Decomposition
Given the potential $V(\mathbf{r})$ in the form:
$$
V(\mathbf{r}) = V_0 e^{-\alpha \mathbf{r}^2}
$$
For a periodic structure, this potential can be expressed as:
$$
V(\mathbf{r}) = \sum_{\mathbf{G}} V_{\mathbf{G}} e^{i \mathbf{G} \cdot \mathbf{r}}
$$
where $V_{\mathbf{G}}$ is the Fourier coefficient of the potential, given as $V_{\mathbf{G}} = V_0 e^{-\alpha \mathbf{G}^2}$.
#### Separation into Harmonic Oscillator and Perturbation Term
1. **Harmonic Oscillator Term**:
The leading term of the Gaussian potential can be approximated by a harmonic oscillator potential:
$$
V_{\text{HO}}(\mathbf{r}) = \frac{1}{2} m \omega^2 \mathbf{r}^2
$$
where $m \omega^2$ is a fitting parameter to match the curvature of the Gaussian potential near the minimum.
2. **Perturbation Term**:
The remaining terms after the harmonic oscillator approximation form the perturbation:
$$
V_{\text{pert}}(\mathbf{r}) = V(\mathbf{r}) - V_{\text{HO}}(\mathbf{r})
$$
#### Schrödinger Equation with Perturbation
The total Hamiltonian can be written as:
$$

H = H_{\text{HO}} + H_{\text{pert}}

$$
where $H_{\text{HO}}$ is the Hamiltonian for the harmonic oscillator and $H_{\text{pert}}$ is the perturbative Hamiltonian.

The Schrödinger equation becomes:
$$
H \psi(\mathbf{r}) = \left( H_{\text{HO}} + H_{\text{pert}} \right) \psi(\mathbf{r}) = E \psi(\mathbf{r})
$$
#### Solutions using Perturbation Theory
1. **Solve the Harmonic Oscillator**:
First, solve the Schrödinger equation for the harmonic oscillator term:
$$
H_{\text{HO}} \psi_n^{(0)}(\mathbf{r}) = E_n^{(0)} \psi_n^{(0)}(\mathbf{r})
$$
The solutions $\psi_n^{(0)}(\mathbf{r})$ are well-known Hermite polynomials $H_n$ multiplied by a Gaussian function.

2. **Apply Perturbation Theory**:
Treat $H_{\text{pert}}$ as a perturbation to find the corrections to the energy levels and wavefunctions.

The first-order correction to the energy is:
$$
E_n^{(1)} = \left\langle \psi_n^{(0)} \left| H_{\text{pert}} \right| \psi_n^{(0)} \right\rangle
$$
The first-order correction to the wavefunction is:
$$
\psi_n^{(1)} = \sum_{m \neq n} \frac{\left\langle \psi_m^{(0)} \left| H_{\text{pert}} \right| \psi_n^{(0)} \right\rangle}{E_n^{(0)} - E_m^{(0)}} \psi_m^{(0)}
$$
#### Gaussian Fourier Coefficients
For the Gaussian potential:
$$
V_{\mathbf{G}} = V_0 e^{-\alpha \mathbf{G}^2}
$$
The perturbation term in reciprocal space:
$$
H_{\text{pert}} = \sum_{\mathbf{G} \neq 0} V_{\mathbf{G}} e^{i \mathbf{G} \cdot \mathbf{r}}
$$
#### Summary
- **Harmonic Oscillator**: $V_{\text{HO}}(\mathbf{r}) = \frac{1}{2} m \omega^2 \mathbf{r}^2$
- **Perturbation Term**: $V_{\text{pert}}(\mathbf{r}) = V(\mathbf{r}) - V_{\text{HO}}(\mathbf{r})$
- **Energy Corrections**:
	- First-order: $E_n^{(1)} = \left\langle \psi_n^{(0)} \left| H_{\text{pert}} \right| \psi_n^{(0)} \right\rangle$
	- Wavefunction corrections: $\psi_n^{(1)}$

- The periodicity of the potential $V(\mathbf{r})$ ensures that both the harmonic oscillator approximation and the perturbation term $V_{\text{pert}}(\mathbf{r})$ are periodic.

- Bloch's theorem applies, allowing the wavefunctions to be written as a product of a plane wave and a periodic function.

- Perturbation theory utilizes these periodic conditions to calculate corrections to energy levels and wavefunctions.