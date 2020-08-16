# PyCahnHilliard
This program solves the Cahn Hilliard equation using the an implicit pseudospectral method. 

The Cahn-Hilliard equation is given by 

<img src="https://render.githubusercontent.com/render/math?math=\frac{\partial c}{\partial t} = M \nabla^2\left[ \frac{\delta F}{\delta \rho}\right] ">

with the free-energy written as 

<img src="https://render.githubusercontent.com/render/math?math=F[c] = \int \text{d}{\mathbf{r}} \left[ \frac{\kappa}{2} (\nabla c(\mathbf{r} ))^2 %2B  f(c)\right]">

where <img src="https://render.githubusercontent.com/render/math?math=\kappa"> is a parameter related to the interfacial energy and <img src="https://render.githubusercontent.com/render/math?math=f"> is the bulk free-energy density given by 

<img src="https://render.githubusercontent.com/render/math?math=f(c) = W c^2(1-c)^2">

where W is the height ot the thermodynamic barrier. 

The following figure is a result for the system with M=1.0, W=2.0, <img src="https://render.githubusercontent.com/render/math?math=\kappa=0.5">, dx=0.1, dt=0.01. The initial condition is given by a normal distribution 

<img src="https://render.githubusercontent.com/render/math?math=c_0 = 0.5%2B 0.1 \mathcal{N}(0,1)">

And the system is evolved until N = 2000 steps. 

![Output](https://github.com/elvissoares/PyCahnHilliard/blob/master/cahn-hilliard-c0-0.5.png)