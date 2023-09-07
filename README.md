# CahnHilliard_PhaseSepration
The python script to solve the Cahn-Hilliard equation using  an implicit pseudospectral algorithm

Author, date: Wentao Gong, On Sep 04, 2023

References:

[Lecture 19 Part VI: Dynamics, the Cahn-Hilliard Equation](https://www.youtube.com/watch?v=O4_6eagE4Gc)

[Phase Field methods: From fundamentals to applications](https://www.youtube.com/watch?v=FTiBq1o-8e4)

[Multi-Phase Fluid Flows: The Cahn-Hilliard Navier-Stokes Framework by Rahul Pandit](https://www.youtube.com/watch?v=0Rmh2OSxayc)
Papers in tutorial folder


## Theory and Concepts
__Cahn-Hilliard Equation:__ This equation models the time evolution of a scalar field $\phi$ which is a fourth-order partial differential equation that models phase separation. Its basic form in real space can be written as:


$$\frac{\partial \phi(x, t)}{\partial t} = \nabla^2 \left( \nabla^2 \phi(x, t) - f'(\phi) \right)$$


Where $f'(\phi)$ is the derivative of the free energy density with respect to $\phi$), the order parameter.

Phase-field model for liquid-liquid phase separation

$$ \frac{\partial c}{\partial t} = M \nabla^2\left[ \frac{\delta F}{\delta c}\right] $$

with **M** being a mobility and the functional of free-energy given by 


$$F(c) = \int [\frac{\kappa}{2}(\nabla c(\mathbf{r})^2+f(c))]d\mathbf{r}$$

Here, $f(c)$ is the bulk free energy density, $\kappa$ is a gradient penalty term (related to interfacial tension), and $c$ is the order parameter (often related to concentration or density). The first term $\frac{\kappa}{2} (\nabla c)^2$ represents the interfacial free energy density. The higher $\kappa$, the higher the energy cost to form an interface.


For the above functional, the functional derivative $\delta F/\delta c$  is commonly found to be:

$$\frac{\delta F}{\delta c} =\frac{\delta f}{\delta c}-\kappa \nabla^2c$$




This expression provides the "force" that acts to change the concentration field $c(\boldsymbol{r},t)$ to minimize the free energy of the system



where **f** is the bulk free-energy density given by 
$$f(c) = Wc^2(1-c)^2$$
The derivative of bulk free energy density
$$\frac{\delta f}{\delta c}=2W[c(1-c)^2+(c-1)c^2]$$
where **W** is the height ot the thermodynamic barrier which must be overcome for a reaction or phase transition to occur. This barrier is essentially the energy difference between the initial state and the highest-energy transition state on the pathway to the final state. For chemical reactions, this is often referred to as the "activation energy". 
Barrier Height is given as
$$W=E_{transition state}- E_{initial state}$$
Higher values of W increase the energetic penalty for mixing, leading to more pronounced phase separation.
### Limitations
1. **Simplicity**:  The concept simplifies complex landscapes of potential energy, especially for systems with multiple possible pathways.
2. **Non-Equilibrium Systems**: The thermodynamic barrier is a concept mainly applicable to systems near or at equilibrium. For systems far from equilibrium, the concept might not be directly applicable.


The next Figure presents this bulk free-energy.

![Bulk](https://github.com/wentaogong111/CahnHilliard_PhaseSepration/blob/main/Report/ch-bulk-free-energy.png)

## Choices of Numerical methods #


In the study of partial differential equations (PDEs) like the Cahn-Hilliard equation, various numerical methods can be employed for solving them, including Fourier methods, Finite Difference Methods (FDM), Finite Volume Methods (FVM), and Finite Element Methods (FEM).

### Spectral Method

### Pros

1. **High Accuracy**: The method approximates the solution in the entire domain simultaneously, offering high-order accuracy for smooth solutions.
2. **Efficiency**: By utilizing Fast Fourier Transforms (FFT), the pseudo-spectral method can be computationally more efficient than finite-difference methods for certain types of problems.
3. **Deal with Nonlinearity**: The Cahn-Hilliard equation often contains nonlinear terms. The pseudo-spectral method is well-equipped to handle such nonlinearities through techniques like operator splitting.
4. **Periodic Boundary Conditions**: For problems where periodic boundary conditions are relevant, the pseudo-spectral method is particularly advantageous.
5. **Reduced Numerical Diffusion**: Compared to finite difference or finite volume methods, pseudo-spectral methods often demonstrate reduced numerical diffusion, which can be crucial for capturing sharp gradients and preserving interface integrity.
6. **Spectral Convergence**: For smooth solutions, spectral methods can offer exponential convergence, making them more accurate than other methods for a given computational cost.

  
### Cons

1. **Limited to Smooth Solutions**: If the solution has discontinuities or sharp gradients, the method might not be effective due to Gibbs phenomenon.
2. **Complexity**: Requires a more sophisticated understanding of the Fourier or other basis functions being used.
3. **Boundary Conditions**: The method is most naturally suited for periodic problems and can be challenging to implement for more complex boundary conditions.
4. **Time-Stepping**: Due to high accuracy in spatial discretization, time-stepping methods must also be chosen carefully to maintain stability.
### Finite Difference Method (FDM)

#### Pros

1. **Simplicity**: Easier to implement and understand.
2. **Flexibility**: Works well for complex geometries and diverse boundary conditions.
3. **Stability**: Known for numerical stability in time-dependent problems.

#### Cons

1. **Accuracy**: Generally less accurate than spectral methods.
2. **Stiff Equations**: May require implicit solvers, increasing complexity.

### Finite Volume Method (FVM)

#### Pros

1. **Conservation Laws**: Excellent for problems where conservation properties need to be preserved.
2. **Robustness**: Good for handling discontinuities.
3. **Unstructured Grids**: Can be used on unstructured, complex geometries.

#### Cons

1. **Computational Cost**: Can be computationally intensive.
2. **Complexity**: Requires the understanding of flux and control volumes.

### Finite Element Method (FEM)

#### Pros

1. **Generalizability**: Excellent for a wide variety of problems, including non-linear and complex geometries.
2. **High Accuracy**: Allows for high-order accuracy.
3. **Adaptive Meshing**: Allows mesh refinement in regions of interest.

#### Cons

1. **Implementation Difficulty**: Complex to implement and requires a deep understanding of the method.
2. **Computational Cost**: High computational expense for three-dimensional problems.

## Conclusion

The choice of method depends on several factors, such as the type of problem you're dealing with, the properties you want to preserve, and the computational resources available.

The spectral method offers a powerful tool for solving the Cahn-Hilliard equation efficiently and accurately under the right conditions. Its ability to deal with complex and nonlinear terms efficiently makes it an appealing choice for researchers working on phase separation problems.

However, if you're dealing with non-periodic boundaries, complex geometries, or need to capture discontinuities accurately, FEM, FDM or FVM might be more appropriate.


## Case1: Pseudo-spectral method

To simplify solving this equation, Fourier transform can be applied. The Fourier Transform of a function **g(x)** is defined as:
$$\widehat{g}(k)=\int e^{ikx}g(x)dx
$$
The concentration field can be expanded as a Fourier series in the form 

$$ 
c(\mathbf{r},t) = \frac{1}{L^2}\sum_{\boldsymbol{k}}\widehat{c}_{\boldsymbol{k}}(t)e^{i\boldsymbol{k}\boldsymbol{r}}
$$

where $\widehat{c}_{\boldsymbol{k}}(t)$ are the Fourier coefficients, and $\boldsymbol{k}= \frac{2\pi N}{L}$ is the wavenumber for a domain of length L.

The Fourier coefficients are given by 
$$
\widehat{c}_{\boldsymbol{k}}(t)=\mathcal {FT}{c(\boldsymbol{r},t)}=\int c(\boldsymbol{r},t)e^{i\boldsymbol{k}\boldsymbol{r}}d\boldsymbol{r}
$$
and $k_i= \{-\pi N_i/L_i, -\pi(N_i-1)/L_i, \ldots, \pi(N_i-1)/L_i,\pi N_i/L_i\}$, where $N_i = L_i/\Delta_i$  and $\Delta_i$ is the stepsize of the meshgrid on the i direction.

The Fourier transform of the dynamical equation is 
$$
\frac{\partial\widehat{c}_{\boldsymbol{k}}}{\partial t}= M[-\kappa k^4\widehat{c}_{\boldsymbol{k}}-k^2\mathcal {FT}\{f'\}]
$$


### Explanation of Terms

1. **$\frac{\partial \hat{c}(k, t)}{\partial t}$**: Rate of change of the Fourier component of $c$ with respect to time $t$.
  
2. **$-k^4 \hat{c}_{\boldsymbol{k}}$**: The Laplacian term in Fourier space. $k^4$ is a consequence of two spatial derivatives $\nabla^2$.

3. **$k^2 \hat{f'}$**: Fourier transform of the derivative of free energy density. $k^2$ arises from one spatial derivative $\nabla^2$.

and using an *implicit* Euler integration, we have
$$
\frac{\widehat{c}_{\boldsymbol{k}}^{n+1}-\widehat{c}_{\boldsymbol{k}}^{n}}{\Delta t}=M[-\kappa k^4 \widehat{c}_{\boldsymbol{k}}^{n+1}-k^2 \mathcal{FT}\{f'(c^n)\}]
$$
Solver for $\widehat{c}_{\boldsymbol{k}}^{n+1}$: isolate $\widehat{c}_{\boldsymbol{k}}^{n+1}$ to find its value at the next time step:
$$
\widehat{c}_{\boldsymbol{k}}^{n+1}=\frac{\widehat{c}_{\boldsymbol{k}}-\Delta tMk^2\mathcal{FT}\{f'(c^n)\} }{1+\Delta t\kappa M k^4}
$$
where $\kappa$ is the gradient coefficient, **k** is the wave vector, M is mobility, $\Delta t$ is the time step value

1. **Inverse Fourier Transform**:Take the inverse Fourier transform of $\widehat{c}_{\boldsymbol{k}}^{n+1}$ to find ${c}_{\boldsymbol{k}}^{n+1}$ the concentration at the next time step in real space.
2. **Repeat**: Go back to  *implicit* Euler scheme and repeat for all desired time steps.

##  Example
The following figure is a result for the system with M=1.0, W=2.0, $\kappa=0.5$, dx=1.0, dt=0.1. The initial condition is given by a normal distribution 
$$
c(\boldsymbol{r},t=0) = c_0 + 0.1 \mathcal{N}(0,1)
$$


And the system is evolved until N = 1000 steps. 
### c0 = 0.2
![GIF](https://github.com/wentaogong111/CahnHilliard_PhaseSepration/blob/main/Results/cahn-hilliard-c0-0.2.png)
### C0 = 0.4

![GIF](https://github.com/wentaogong111/CahnHilliard_PhaseSepration/blob/main/Results/ch-c0%3D0.4.gif)

### C0 = 0.5

![GIF](https://github.com/wentaogong111/CahnHilliard_PhaseSepration/blob/main/Results/ch-c0%3D0.5.gif)

### C0 = 0.7

![GIF](https://github.com/elvissoares/PyCahnHilliard/blob/master/ch-c0%3D0.7.gif)
