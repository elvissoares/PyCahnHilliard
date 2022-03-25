# PyCahnHilliard
This program solves the Cahn Hilliard equation using the an implicit pseudospectral method. 

The Cahn-Hilliard equation is defined as

<img src="https://render.githubusercontent.com/render/math?math=\frac{\partial c}{\partial t} = M \nabla^2\left[ \frac{\delta F}{\delta c}\right] ">

with M being a mobility and the functional of free-energy given by 

<img src="https://render.githubusercontent.com/render/math?math=F[c] = \int \text{d}{\mathbf{r}} \left[ \frac{\kappa}{2} (\nabla c(\mathbf{r} ))^2 %2B  f(c)\right]">

where <img src="https://render.githubusercontent.com/render/math?math=\kappa"> is a parameter related to the interfacial energy and <img src="https://render.githubusercontent.com/render/math?math=f"> is the bulk free-energy density given by 

<img src="https://render.githubusercontent.com/render/math?math=f(c) = W c^2(1-c)^2">

where W is the height ot the thermodynamic barrier. 

## Pseudo-spectral method

The concentration field can be expanded as a Fourier series in the form 

<img src="https://latex.codecogs.com/svg.image?\inline&space;\displaystyle&space;c(\boldsymbol{r},t)&space;=&space;\frac{1}{L^2}&space;\sum_{\boldsymbol{k}}&space;\widehat{c}_{\boldsymbol{k}}(t)&space;e^{i&space;\boldsymbol{k}&space;\cdot&space;\boldsymbol{r}&space;}" title="https://latex.codecogs.com/svg.image?\inline \displaystyle c(\boldsymbol{r},t) = \frac{1}{L^2} \sum_{\boldsymbol{k}} \widehat{c}_{\boldsymbol{k}}(t) e^{i \boldsymbol{k} \cdot \boldsymbol{r} }" />

where the Fourier coefficients are given by 

<img src="https://latex.codecogs.com/svg.image?\inline&space;\displaystyle&space;\widehat{c}_{\boldsymbol{k}}(t)&space;=&space;\mathcal{FT}\{c(\boldsymbol{r},t)&space;\}&space;=&space;\int&space;d&space;\boldsymbol{r}\&space;c(\boldsymbol{r},t)e^{-i&space;\boldsymbol{k}&space;\cdot&space;\boldsymbol{r}&space;}" title="https://latex.codecogs.com/svg.image?\inline \displaystyle \widehat{c}_{\boldsymbol{k}}(t) = \mathcal{FT}\{c(\boldsymbol{r},t) \} = \int d \boldsymbol{r}\ c(\boldsymbol{r},t)e^{-i \boldsymbol{k} \cdot \boldsymbol{r} }" />

and <img src="https://latex.codecogs.com/svg.image?\inline&space;k_i&space;=&space;\{-\pi&space;N_i/L_i,&space;-\pi(N_i-1)/L_i,&space;\ldots,&space;\pi(N_i-1)/L_i,\pi&space;N_i/L_i\}" title="https://latex.codecogs.com/svg.image?\inline k_i = \{-\pi N_i/L_i, -\pi(N_i-1)/L_i, \ldots, \pi(N_i-1)/L_i,\pi N_i/L_i\}" /> where <img src="https://latex.codecogs.com/svg.image?\inline&space;N_i&space;=&space;L_i/\Delta_i" title="https://latex.codecogs.com/svg.image?\inline N_i = L_i/\Delta_i" /> and <img src="https://latex.codecogs.com/svg.image?\inline&space;\Delta_i" title="https://latex.codecogs.com/svg.image?\inline \Delta_i" /> is the stepsize of the meshgrid on the i direction.

The Fourier transform of the dynamical equation is 

<img src="https://latex.codecogs.com/svg.image?\frac{\partial&space;\widehat{c}_{\boldsymbol{k}}}{\partial&space;t}&space;=&space;M&space;\left[&space;-k^2&space;\mathcal{FT}\{f'\}&space;-&space;\kappa&space;k^4&space;\widehat{c}_{\boldsymbol{k}}\right]&space;">

and using an *implicit* Euler integration, we have

<img src="https://latex.codecogs.com/svg.image?\frac{\widehat{c}_{\boldsymbol{k}}^{n&plus;1}-\widehat{c}_{\boldsymbol{k}}^{n}}{\Delta&space;t}&space;=&space;M&space;\left[&space;-k^2&space;\mathcal{FT}\{f'(c^n)\}&space;-&space;\kappa&space;k^4&space;\widehat{c}_{\boldsymbol{k}}^{n&plus;1}&space;\right]&space;">

such that 

<img src="https://latex.codecogs.com/svg.image?\widehat{c}_{\boldsymbol{k}}^{n&plus;1}=&space;\frac{&space;\widehat{c}_{\boldsymbol{k}}-&space;\Delta&space;t&space;M&space;k^2&space;\mathcal{FT}\{f'(c^n)\}}{1&plus;&space;\Delta&space;t&space;\kappa&space;k^4&space;}" title="https://latex.codecogs.com/svg.image?\widehat{c}_{\boldsymbol{k}}^{n&plus;1}=&space;\frac{&space;\widehat{c}_{\boldsymbol{k}}-&space;\Delta&space;t&space;M&space;k^2&space;\mathcal{FT}\{f'(c^n)\}}{1&plus;&space;\Delta&space;t&space;\kappa&space;k^4&space;}" />

where <img src="https://latex.codecogs.com/svg.image?\inline&space;\Delta&space;t&space;"> is the time step value. 

## Example

The following figure is a result for the system with M=1.0, W=2.0, <img src="https://render.githubusercontent.com/render/math?math=\kappa=0.5">, dx=0.1, dt=0.01. The initial condition is given by a normal distribution 

<img src="https://latex.codecogs.com/svg.image?c(\boldsymbol{r},t=0)&space;=&space;c_0&space;&plus;&space;0.1&space;\mathcal{N}(0,1)" title="https://latex.codecogs.com/svg.image?c(\boldsymbol{r},t=0) = c_0 + 0.1 \mathcal{N}(0,1)" />

And the system is evolved until N = 1000 steps. 

### C0 = 0.3

![GIF](https://github.com/elvissoares/PyCahnHilliard/blob/master/ch-c0%3D0.3.gif)

### C0 = 0.5

![GIF](https://github.com/elvissoares/PyCahnHilliard/blob/master/ch-c0%3D0.5.gif)

### C0 = 0.7

![GIF](https://github.com/elvissoares/PyCahnHilliard/blob/master/ch-c0%3D0.7.gif)
