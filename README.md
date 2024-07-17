# PyCahnHilliard
This program solves the Cahn-Hilliard equation using an implicit pseudospectral method. 

The Cahn-Hilliard equation is defined as

$$\frac{\partial c}{\partial t} = M \nabla^2\left[ \frac{\delta F}{\delta c}\right] = M \left[-\kappa \nabla^4 c + \nabla^2 f'(c)\right]$$

with M being a mobility and the functional of free-energy given by 

$$F[c] = \int \left[ \frac{\kappa}{2} (\nabla c(\boldsymbol{r} ))^2 + f(c)\right] \text{d}{\boldsymbol{r}} $$

where $\kappa$ is a parameter related to the interfacial energy and $f$ is the bulk free-energy density given by 

$$f(c) = W c^2(1-c)^2$$

where $W$ is the height of the thermodynamic barrier. The next Figure presents this bulk free-energy.

![Bulk](https://github.com/elvissoares/PyCahnHilliard/blob/master/ch-bulk-free-energy.png)

## Pseudo-spectral method

The concentration field can be expanded as a Fourier series in the form 

$$\displaystyle c(\boldsymbol{r},t) = \frac{1}{L^2} \sum_{\boldsymbol{k}} \widehat{c}_{\boldsymbol{k}}(t) e^{i \boldsymbol{k} \cdot \boldsymbol{r} }$$

where the Fourier coefficients are given by 

$$\widehat{c}_{\boldsymbol{k}}(t) = \mathcal{FT}[c(\boldsymbol{r},t) ] = \int_V  c(\boldsymbol{r},t)e^{-i \boldsymbol{k} \cdot \boldsymbol{r} }\text{d}{\boldsymbol{r}} $$

and $k_i = \{-\pi N_i/L_i, -\pi(N_i-1)/L_i, \ldots, \pi(N_i-1)/L_i,\pi N_i/L_i\}$ where $\Delta_i$ is the gridsize of the meshgrid on the $i$ direction.

The Fourier transform of the dynamical equation is 

```math
\frac{\partial}{\partial t} \widehat{c}_{\boldsymbol{k}}  = M [ - k^2 \mathcal{FT}[f']-\kappa k^4 \widehat{c}_{\boldsymbol{k}} ]
``` 

and using an *implicit* Euler time integration, we have

$$ \frac{\widehat{c}_{\boldsymbol{k}}^{n+1} -\widehat{c}_{\boldsymbol{k}}^{n} }{\Delta t}=M\left [-k^2\mathcal{FT}[f'(c^n)]-\kappa k^4 \widehat{c}_{\boldsymbol{k}}^{n+1} \right ] $$

such that 

$$ \widehat{c}_{\boldsymbol{k}}^{n+1} =\frac{\widehat{c}_{\boldsymbol{k}}^n -\Delta t M k^2 \mathcal{FT}[f'(c^n)]}{1 +\Delta t \kappa k^4} $$

where $\Delta t$ is the time stepsize. 

## Example

The following figures are results from the CH equations for a system with M=1.0, W=2.0, $\kappa=0.5$ and three different initial conditions $c_0 = 0.3, 0.5, 0.7$. The gridsize is $L = 64\pi$ with the number of gridpoints $N = 2^9 = 512$. The initial condition is given by a normal distribution 

$$c(\boldsymbol{r},t=0) = c_0 + 0.1 \mathcal{N}(0,1),$$

and our system is evolved during 10000 steps with stepsize of dt=0.1.

### C0 = 0.3

![GIF](https://github.com/elvissoares/PyCahnHilliard/blob/master/ch-c0%3D0.3.gif)

### C0 = 0.5

![GIF](https://github.com/elvissoares/PyCahnHilliard/blob/master/ch-c0%3D0.5.gif)

### C0 = 0.7

![GIF](https://github.com/elvissoares/PyCahnHilliard/blob/master/ch-c0%3D0.7.gif)

## Options and Dependences

In *cahnhilliard.py* we use just the Numpy package to do the fft. In *cahnhilliard-pytorch.py* we use the torch package to do the fft using the CUDA-capable NVIDIA GPU. 

* [NumPy](https://numpy.org) is the fundamental package for scientific computing with Python.
* [PyTorch](https://pytorch.org/) is a high-level library for machine learning, with multidimensional tensors that can also be operated on a CUDA-capable NVIDIA GPU. 
* [Matplotlib](https://matplotlib.org/stable/index.html) is a comprehensive library for creating static, animated, and interactive visualizations in Python.


# Cite My work

If you use *cahnhilliard.py* or *cahnhilliard-pytorch.py* in your work, please consider to cite it using the following reference:

Soares, E. do A., Barreto, A. G. & Tavares, F. W. *Exponential Integrators for Phase-Field Equations using Pseudo-spectral Methods: A Python Implementation.* 1–12 (2023). ArXiv: [2305.08998](https://arxiv.org/abs/2305.08998)

Bibtex:

    @article{Soares2023,
    archivePrefix = {arXiv},
    arxivId = {2305.08998},
    author = {Soares, Elvis do A. and Barreto, Amaro G. and Tavares, Frederico W},
    eprint = {2305.08998},
    month = {may},
    pages = {1--12},
    title = {{Exponential Integrators for Phase-Field Equations using Pseudo-spectral Methods: A Python Implementation}},
    url = {http://arxiv.org/abs/2305.08998},
    year = {2023}
    }


# Contact
Elvis Soares: elvis.asoares@gmail.com

Universidade Federal do Rio de janeiro

School of Chemistry
