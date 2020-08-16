# PyCahnHilliard
This program solves the Cahn Hilliard equation using the an implicit pseudospectral method. 

The Cahn-Hilliard equation is given by 
$$ \frac{\partial c}{\partial t} = M \nabla^2\left[ \frac{\delta F}{\delta \rho}\right] $$
white the free-energy written as 
$$ F[c] = \int \text{d}{\mathbf{r}} \left[ \frac{\kappa}{2} (\nabla c(\mathbf{r} ))^2 + f(c)\right] $$
where $\kappa$ is a parameter related to the interfacial energy and $f$ is the bulk free-energy density given by 
$$f(c) = W c^2(1-c)^2$$
where $W$ is the height ot the thermodynamic barrier. 