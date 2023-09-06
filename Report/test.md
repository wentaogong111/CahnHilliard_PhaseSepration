# Interfacial Free Energy Density: Equation and Theory

## Theory

The interfacial free energy density is a measure of the energetic cost of creating an interface between two phases. This is a crucial concept in understanding phenomena such as phase separation, domain formation, and other complex morphological behaviors in systems like polymers, colloids, or binary alloys.

## The Landau-Ginzburg Functional

One of the common theoretical frameworks to describe this is the Landau-Ginzburg functional, usually denoted \( F \), which for a one-component system can be expressed as:

$$
F = \int d^3r \left[ \frac{\kappa}{2} (\nabla \phi)^2 + f(\phi) \right]
$$

Here, $f(\phi)$ is the bulk free energy density, $\kappa$ is a gradient penalty term (related to interfacial tension), and $\phi$ is the order parameter (often related to concentration or density). The first term $\frac{\kappa}{2} (\nabla \phi)^2$ represents the interfacial free energy density. The higher $\kappa$, the higher the energy cost to form an interface.

## Equations in Fourier Space

In Fourier space, the interfacial term can be represented as:

$$
F_{\text{interf}} = \frac{\kappa}{2} \int \text{d}^3k \, |\phi(k)|^2 k^2
$$
Where $k$ is the wave vector, and $\phi(k)$ is the Fourier transform of $\phi(r)$.

Phase-field model for liquid-liquid phase separation

$$ \frac{\partial c}{\partial t} = M \nabla^2\left[ \frac{\delta F}{\delta c}\right] $$
the functional of free-energy given by
$$F(c) = \int [f(c))+\frac{\kappa}{2}(\nabla c(\mathbf{r})^2]d\mathbf{r}$$
and
$$\frac{\delta F}{\delta c} =\frac{\delta f}{\delta c}-\kappa \nabla^2c$$

the bulk free-energy density given by 

$$f(c) = Wc^2(1-c)^2$$
The derivative of bulk free energy density
$$\frac{\delta f}{\delta c}=2W[c(1-c)^2+(c-1)c^2]$$

The Fourier coefficients are given by 
$$
\widehat{c}_{\boldsymbol{k}}(t)=\mathcal {FT}{c(\boldsymbol{r},t)}=\int c(\boldsymbol{r},t)e^{i\boldsymbol{k}\boldsymbol{r}}d\boldsymbol{r}
$$
and $k_i= \{-\pi N_i/L_i, -\pi(N_i-1)/L_i, \ldots, \pi(N_i-1)/L_i,\pi N_i/L_i\}$, where $N_i = L_i/\Delta_i$  and $\Delta_i$ is the stepsize of the meshgrid on the i direction.

The Fourier transform of the dynamical equation is 
$$
\frac{\partial\widehat{c}_{\boldsymbol{k}}}{\partial t}= M[-\kappa k^4\widehat{c}_{\boldsymbol{k}}-k^2\mathcal {FT}\{\frac{\delta f}{\delta c}\}]
$$

The derivative of bulk free energy density
$$\frac{\delta f}{\delta c}=2W[c(1-c)^2+(c-1)c^2]$$
$\frac{\kappa}{2} (\nabla c)^2$ represents the interfacial free energy density.

$$
F_{\text{interf}} = \frac{\kappa}{2} \int |\widehat{c}_k|^2 k^2\text{d}^3k \
$$
The interfacial energy density 
$$
f_{\text{interf}} = \kappa|\widehat{c}_k|^2 k
$$
