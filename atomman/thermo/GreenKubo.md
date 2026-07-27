# Green-Kubo methods

The Green-Kubo (GK) methods are a family of equilibrium simulation analyses that can be used to estimate different transport properties (diffusion, viscosity, and thermal conductivity) of a material. At the center of all of the GK methods is an integral of an auto-correlation function computed from a measured property. Evaluating the integral is the trickest part of the GK methods.

## Auto-correlation functions

The auto-correlation function, ACF, of some property $X(t)$ is $\left<X_0 \cdot X_t\right>$, where $X_0$ and $X_t$ are measurements of $X(t)$ taken at different times separated by some $t$.  ACF can then be evaluated by varying the $t$ from which $\left<X_0 \cdot X_t\right>$ is evaluated. As a note of caution, ACF can be evaluated differently across different fields, so it is important that the evaluation is done in the correct way for GK calculations!  Notably, ACF should not be normalized and have units that are consistent with $X^2$, and measurements should be divided by the number of samples, i.e. be sure it is $\left<X_0 \cdot X_t\right>$ and not just $X_0 \cdot X_t$!

All Green-Kubo methods are based on the integral of an auto-correlation function with respect to time,

$$ I = \int_0^\infty{\left<X_0 \cdot X_t\right>dt}. $$

The particular property $X$ that is evaluated depends on which materials property is being estimated; diffuision ($D$), viscosity ($\mu$), or thermal conductivity ($\kappa$)

$$ D = \frac{1}{d}\int_0^\infty{\left<\mathbf{v}_0 ⋅ \mathbf{v}_t\right>dt} $$

$$ \mu = \frac{V}{d k_B T} \int_0^{\infty}{\left<\mathbf{P}_0 ⋅ \mathbf{P}_t\right>dt} $$

$$ \kappa = \frac{V}{d k_B T^2} \int_0^{\infty}{\left< \mathbf{J}_0 ⋅ \mathbf{J}_t\right>dt}$$

In these equations, $d$ is the number of dimensions included in the estimate, $V$ is volume, $k_B$ is the Boltzmann constant, $T$ is temperature, $\mathbf{v}$ are atomic velocities, $\mathbf{P}$ are shear pressures, and $\mathbf{J}$ are heat fluxes.

As the ACF measures the correlation of $X$ at a moment $t$ to itself at an earlier moment in time, for complex chaotic systems we can expect that there will be some value of $t$ beyond which $\left<X_0 \cdot X_t\right> = 0$. Thus, the integral $I$ should be a finite constant. In practice, however, the measured ACF data tends to contain considerable fluctuations and may converge around a non-zero constant, both of which can lead to the integral diverging as time increases rather than converging. Two ways to improve the integral calculation are to increase the number of samples of $X_0 \cdot X_t$ to smooth the mean, and to select a cutoff time for the integral, $t_c$, to ignore values beyond where the ACF is expected to be zero.

### Increasing the number of samples

Since the ACF involves averaging $X_0 \cdot X_t$ values, the fluctuations in the computed ACF curves can be reduced by increasing the number of measurements of $X_0 \cdot X_t$ for a given $t$. The three ways of going about this are to increase the number of atoms in the system, increase how long the simulations that measure $X(t)$ are ran, and increase how many simulations are performed. All three increase the total computational time and have their own strengths and weaknesses.

Increasing the system size, i.e. the number of atoms, affects the calculations differently based on what property $X$ is being evaluated. For $\mathbf{v}$, increasing the number of atoms directly increases the number of $\mathbf{v}_0 \cdot \mathbf{v}_t$ samples as it can be independently evaluated for each atom. In contrast, single values of $\mathbf{P}$ and $\mathbf{J}$ are measured for the whole system at each $t$.  Increasing the number of atoms improves the measurements of $\mathbf{P}$ and $\mathbf{J}$, but it does not increase the number of $\mathbf{P}_0 \cdot \mathbf{P}_t$ or $\mathbf{J}_0 \cdot \mathbf{J}_t$ samples. As such, increasing the number of atoms likely improves the results but the relationship in measurement quality is not as clear as it is for $\mathbf{v}$. One additional consideration for system sizes is that smaller systems will be more constrained than larger systems, for example $\mathbf{J}$ depends on phonon fluctuations which may be constrained if the system dimensions are too small.  Be sure to try different system sizes to see if there is a size affect!

How long each simulation used to evaluate the properties $X$ also influences the measurement quality. At the bare minimum, the simulation time should be large enough that the ACF does approach and remain at/near zero to capture the full integral space. For the properties $\mathbf{P}$ and $\mathbf{J}$, running the calculation for longer times results in a continuous improvement of the ACF as it increases the number of $X_0 \cdot X_t$ samples for each $t$. This is because $\mathbf{P}$ and $\mathbf{J}$ are evaluated and saved for each $t$ allowing for the $X_0$ value to roll through the table rather than being fixed to the first sample only. Thus, if you run for a total of 1000 steps of some $Δt$, you would end up with 999 samples of $X_0 \cdot X_{Δt}$, 998 samples of $X_0 \cdot X_{2Δt}$, and so on. 

This, however, is not the case for $\mathbf{v}$ as the velocity auto-correlation function (VACF) calculation method *does* fix $\mathbf{v}_0$ to be for only the first time measurement. Since $\mathbf{v}$ is computed per atom, keeping track of each atom's velocity for each time step can quickly become memory intensive. By fixing $\mathbf{v}_0$, then only the $t=0$ reference state has to be stored in order to compute $\left<\mathbf{v}_0 ⋅ \mathbf{v}_t\right>$ for each subsequent time $t$. Assuming the max time involved is beyond the VACF tail, increasing the calculation time will not improve the integral calculation in any way.

Running multiple independent simulations is also a beneficial means for improving the ACF and the resulting materials property estimate. With independent simulation results, the ACF measurements for each $t$ can be averaged together to produce a smoother ACF curve, or the final computed property estimates can be averaged together. Performing separate simulations also provides a simple means of evaluating the calculation error from the standard deviation of the outputted values. 

Interestingly, the final property value estimated from the combined smooth ACF will be identical to the averaged property estimate if the simulations were performed for the same number of atoms for the same number of time steps, and use the same time cutoff for the integral. This happens because averages of the same sample sizes can be averaged together without weighting, and the integral is itself a summation.

For $\mathbf{v}$, increasing the number of simulations offers an improvement in the ACF on par with increasing the number of atoms. For $\mathbf{P}$ and $\mathbf{J}$, it is actually better to increase the simulation run time rather than running separate simulations. This is due to the number of $X_0 \cdot X_t$ samples taken for subsequent $t$ values decreasing by the number of separate simulations. Thus, running one long simulation provides much better estimates for moderate $t$ values than averaging across independent simulations.

### Selecting a cutoff time for the integral

Alternatively, the integral calculation can be greatly improved through selection of a $t_c$ that appropriately excludes the noise from the ACF function after it ideally should be zero. The challenge is in finding an appropriate method that can freely identify $t_c$ from the ACF curve alone. The method adapted here is to evaluate a fluctuation f defined as

$$ f(t) = \frac{\sigma(ACF)}{\left<ACF\right>}, $$

where the standard deviation $\sigma$ and mean of ACF are evaluated for rolling sample windows of roughly 1000 samples. $t_c$ is then automatically selected based on $f(t)$ exceeding some critical value, which for simplicity is often taken as 1. Thus, $t_c$ is found at the $t$ where the fluctuations in ACF exceed the mean value. 

Initial tests have shown that the fluctuation method does provide a decent automated means of identifying $t_c$ values resulting in mostly consistent integral values. However, there are some downsides in that $f(t)$ can be sensitive to the rolling window size and setting the cutoff to $f(t)=1$ is arbitrary. In practice, I've so far found that a better cutoff is $\max(2f(0), 1)$ as the standard deviation tends to be large at the beginning of the ACF curves due to a rapid initial drop-off.  With the rolling window size, 1000 seems decent for $\mathbf{J}$, while a much smaller 100 is better for $\mathbf{P}$.

As for $\mathbf{v}$, 