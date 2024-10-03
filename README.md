<div id="top"></div>

# Matrix method in Python for acoustic levitation simulations
<a href="https://www.codefactor.io/repository/github/slkiser/acousticlevitation"><img src="https://www.codefactor.io/repository/github/slkiser/acousticlevitation/badge" alt="CodeFactor" /></a> <img alt="GitHub" src="https://img.shields.io/github/license/slkiser/acousticLevitation"> <img src="https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2Fslkiser%2FacousticLevitation&count_bg=%23FFB031&title_bg=%23555555&icon=&icon_color=%23E7E7E7&title=hits&edge_flat=false"/> <a href="https://www.linkedin.com/in/shawn-kiser/"> <img src="https://img.shields.io/badge/LinkedIn-blue?style=flat&logo=linkedin&labelColor=grey"></a>

This repository implements the matrix method used in the publication ["Matrix Method for Acoustic Levitation Simulation"](https://www.researchgate.net/publication/224254694_Matrix_Method_for_Acoustic_Levitation_Simulation). Specifically, the results of Section II are written in Python, which describes an acoustic levitator composed of only one circular flat transducer (and with a hole if needed) and one planar reflector.

This [README](README.md) consists of an introduction and theory which derives the **Rayleigh integral equation** used in the publication. This integral is then numerically resolved using the matrix method used in the publication. The matrix method offers an advantage over traditional [numerical integration](https://en.wikipedia.org/wiki/Numerical_integration) for two-dimensional (2D) problems, e.g., the [trapezoidal rule](https://en.wikipedia.org/wiki/Trapezoidal_rule) (1st order approximation) & [Simpson's rule](https://en.wikipedia.org/wiki/Simpson%27s_rule) (2nd order approximation), by obtaining a higher order of accuracy at a lower computation cost.

If you use this code, please cite this work:
> [1] S. L. Kiser, Matrix method in Python for acoustic levitation simulations. [Online]. Available: https://github.com/slkiser/acousticLevitation

## Introduction
Acoustics is a field that has ties to mechanics and fluid dynamics. For solid mechanics, a phenomenological approach to vibration in solids is adequate to describe modal shapes and eigenfrequencies. However, it is not adequate for relating mechanical vibrations to noise radiation and transmission. Therefore, the interaction of sound waves and the vibration of solids requires a fundamental wave approach. Consider the following:

- Solids can store potential energy in compression (longitudinal), flexural (bending/transverse), shear, and torsional waves.
- Fluids can only store energy in compression, thus they can only sustain compressional (longitudinal) waves.

We are interested in an approach is required that can relate the wave-particle velocities of the vibrating solid surface to the volume velocities in the fluid. In other words, we seek an approach that relates the energy exchange between the solid and the fluid. 

At a fundamental level, the radiation of sound from a vibrating boundary (surface) can be formulated in terms of an integral equation involving [Green's functions](https://en.wikipedia.org/wiki/Green%27s_function) with imposed radiation constraints, i.e., the radiation condition ensures that the integral equation for the radiated sound pressure represents outward-traveling sound waves.

In its most general form, the integral equation is attributable to Kirchhoff, although Helmholtz modified it for single-frequency (harmonic) applications. This integral, [Kirchhoff-Helmholtz integral equation](https://en.wikipedia.org/wiki/Kirchhoff%E2%80%93Helmholtz_integral), relates the harmonic vibrational motion of the surface to the radiated sound pressure field in the enclosed fluid. It can be interpreted as representing the sound pressure field of a vibrating surface by a distribution of volume velocity sources and forces on the surface. The velocity sources and the forces are related to normal surface velocity and surface pressure respectively. It is important to note that the surface pressure and the normal surface vibrational velocity interact with one another. 

Rayleigh modified the [Kirchhoff-Helmholtz integral equation](https://en.wikipedia.org/wiki/Kirchhoff%E2%80%93Helmholtz_integral) for the specific case of a planar source located in an infinite baffle and illustrated that it is equivalent to a distribution of point sources. This modification named the  **Rayleigh integral equation** is outlined in the following subsection.

### Theory (need to view on white background)
The nature of **steady-state sound fields** in presence of boundaries (surfaces) is described by the [Helmholtz equation](https://en.wikipedia.org/wiki/Helmholtz_equation):

```math
\large \nabla^2 p + k^2 p = 0
```

To describe the sound pressure  $$p(\mathbf{r},\omega)$$ in [polar coordinates](https://ocw.mit.edu/courses/aeronautics-and-astronautics/16-07-dynamics-fall-2009/lecture-notes/MIT16_07F09_Lec05.pdf) due to the boundary vibration and reflection, a [Green's function](https://en.wikipedia.org/wiki/Green%27s_function) is used to satisfy the following equation:

```math
\nabla^{2} G\left(\mathbf{r} \mid \mathbf{r}'\right)+k^{2} G\left(\mathbf{r} \mid \mathbf{r}'\right)=-\delta\left(\mathbf{r}-\mathbf{r}'\right)
```

Kirchhoff obtained the solution by substituting both previous equations for sound pressure in the integral form (the [Kirchhoff-Helmholtz integral](https://en.wikipedia.org/wiki/Kirchhoff%E2%80%93Helmholtz_integral)):

```math
p(\mathbf{r}, \omega)=\oint_{s}\left[G\left(\mathbf{r} \mid \mathbf{r}'\right) \frac{\partial p}{\partial n'}-p\left(\mathbf{r}', \omega\right) \frac{\partial G\left(\mathbf{r} \mid \mathbf{r}'\right)}{\partial n'}\right] \mathrm{d} s'
```

where $$\partial p / \partial n'$$ represents the gradient on the boundary surface (outward positive from the sound generation).

The boundary conditions about the plate, $$S_1$$, and the rigid reflector, $$S_2$$, are:

```math
\mathbf{v}_{B} = 
\begin{cases} 
\mathbf{v}_{n}\left(\mathbf{r}_{p}, \omega\right) \mathrm{e}^{j \omega t} & \text{if } \mathbf{r}_{p}=\left(x', y'\right) \in S_{1} \\ 
0 & \text{if } \mathbf{r}_{p}=\left(x', y'\right) \in S_{2} 
\end{cases}
```

where the polar velocity is $$\mathbf{v}=d\mathbf{r}/dt$$.

By applying [Euler's equation](https://en.wikipedia.org/wiki/Euler_equations_(fluid_dynamics)), or a **balance of momentum** ([Newton's second law](https://en.wikipedia.org/wiki/Newton%27s_laws_of_motion#Newton's_second_law)), the pressure gradient on the boundary surface is equal to the boundary vibration velocity:

```math
\frac{\partial p}{\partial n'}=j \rho_{0} \omega \mathbf{v}_{B}
```

For this boundary, a [Green's function](https://en.wikipedia.org/wiki/Green%27s_function) is needed to satisfy:

```math
\frac{\partial G\left(\mathbf{r} \mid \mathbf{r}'\right)}{\partial n'}=0
```

at the boundary surface (when $$z=0$$ ). The [Green's function](https://en.wikipedia.org/wiki/Green%27s_function) selected by Rayleigh is:

```math
G\left(\mathbf{r} \mid \mathbf{r}'\right)=\frac{1}{4 \pi}\left(\frac{\mathrm{e}^{-j k R_{+}}}{R_{+}}+\frac{\mathrm{e}^{-j k R_{-}}}{R_{-}}\right)
```
where:

```math
\begin{array}{l}
R_{+}=\sqrt{\left(x-x'\right)^{2}+\left(y-y'\right)^{2}+\left(z-z'_{+}\right)^{2}} \\
R_{-}=\sqrt{\left(x-x'\right)^{2}+\left(y-y'\right)^{2}+\left(z+z'_{-}\right)^{2}}
\end{array}
```

with $$z'_{+}$$ and $$z'_{-}$$ being small positive and negative perturbations in the $$Z$$ direction. This [Green's function](https://en.wikipedia.org/wiki/Green%27s_function) has a zero gradient on the boundary of $$z=0$$:

```math
\left.\frac{\partial G\left(\mathbf{r} \mid \mathbf{r}'\right)}{\partial z'}\right|_{\{z'_{+}, z'_{-}\}\rightarrow 0}=0
```

and when the sources are limited to the boundary surface $$\mathbf{r}'\rightarrow \mathbf{r}_{p}$$, Rayleigh gets the second Green's function:

```math
G\left(\mathbf{r} \mid \mathbf{r}_{p}\right)=\frac{1}{2 \pi} \frac{\mathrm{e}^{-j k R}}{R}
```

where:

```math
R=\left|\mathbf{r}-\mathbf{r}_{p}\right|=\sqrt{\left(x-x'\right)^{2}+\left(y-y'\right)^{2}+z^{2}}
```

Substituting the past 3 equations into the [Kirchhoff-Helmholtz integral](https://en.wikipedia.org/wiki/Kirchhoff%E2%80%93Helmholtz_integral) results in the **Rayleigh integral** for the sound pressure field of a vibration source and a planar reflector:
```math
\boxed{p(\mathbf{r}, \omega)=\frac{j \rho_{0} \omega}{2 \pi} \int_{S_{1}} \frac{v_{n}\left(\mathbf{r}_{p}, \omega\right)}{R} \mathrm{e}^{-j k R} \mathrm{d} s}
```

The **Rayleigh integral** is then solved using the matrix method outlined in the publication ["Matrix Method for Acoustic Levitation Simulation"](https://www.researchgate.net/publication/224254694_Matrix_Method_for_Acoustic_Levitation_Simulation). Numerically, the **Rayleigh integral** is determined by discretizing the transducer and the reflector surfaces in small area cells, and multiplication of the transfer matrices gives the name of the matrix method.


## Setup and usage
I recommend the usage of the [Anaconda distribution](https://www.anaconda.com/products/individual) and the included [Spyder IDE](https://www.spyder-ide.org/). Spyder is a unique Python IDE that offers a [MATLAB](https://www.mathworks.com/products/matlab.html)-like experience with an integrated variable viewer.

### Package installation
For required packages, refer to [`requirements.txt`](requirements.txt). The packages can be installed quickly using pip and/or in a virtual environment:
```
pip install -r requirements.txt
```
Or if using Anaconda the requirements are already included.

### Validation
[`validation.py`](validation.py) is a script that recreated Fig. 2 (a) in [1]. It creates a grid of 18mm by 40mm, with a flat circular transducer of radius `R = 5*mm`, and a flat reflector along the bottom of the grid. A side by side comparison is presented:
Original plot, Fig 2 (a) of [1]            |  [`validation.py`](validation.py)
:-------------------------:|:-------------------------:
![Image of original plot](https://github.com/slkiser/acousticLevitation/blob/main/fig%202a.png)  |  ![Image of recreated plot](https://github.com/slkiser/acousticLevitation/blob/main/validation.png)

### How to use

*Make sure to make the main folder the console's working directory.*

[`main.py`](main.py) is a script that generates a grid of 50mm by 50mm, with a flat circular transducer of radius `R = 15*mm` with hole `R2=2*mm`, and a flat reflector along the bottom of the grid. All modifiable parameters are indicated on lines 51-76:

```
#%% Definition of grid and problem (modify for your needs)
x_min = -50*mm
x_max = 50*mm
z_min = 0                               # Location of the reflector
z_max = 50*mm                           # Location of the transducer
disc = 0.5*mm 	                        # Discretization step size 

# Define air constants
rho = complex(1.214, 0)			        # Density of air in Raleigh, kg/m^3
c = complex(340.1, 0)			        # Speed of sound in air, m/s

# Define transducer parameters
R = 15*mm 		                        # Transducer radius, m
A_trans = complex(np.pi*R**2, 0)		# Transducer area, m^2
R2 = 2*mm                               # Trandsucer hole, m
A_hole = complex(np.pi*R2**2, 0)        # Trandsucer hole area, m^2
f = complex(56000, 0)			        # Frequency, Hz
U_0 = 0.0000060                         # Displacement amplitude, m

# Define reflector parameters
A_refl = complex(np.pi*(x_max/2)**2, 0) # Reflector area, m^2

# Constants
wavelength = complex(c/f, 0)		    # Wavelength, m
omega = complex(2*np.pi*f, 0)		    # Angular frequency, rad/s
wavenumber = omega/c			        # Wavenumber - k
E = complex(0, 1)/wavelength		    # Constant used for reflected waves
D = (omega*rho*c)/wavelength		    # Constant used for transmitted wave
```

The script calculates the distance arrays shown in [figure 1 of the publication](https://www.researchgate.net/publication/224254694_Matrix_Method_for_Acoustic_Levitation_Simulation) `r_nm`, `r_im`, and `r_in` on lines 92-116. The transfer matrices of [eq. (2) and eqs.(4-6)](https://www.researchgate.net/publication/224254694_Matrix_Method_for_Acoustic_Levitation_Simulation) `T_TM`,`T_TR`,`T_RT`, and `T_RM` are calculated on lines 125-158. The pressure calculation of [eq. (3)](https://www.researchgate.net/publication/224254694_Matrix_Method_for_Acoustic_Levitation_Simulation) `P` is transcribed identically (up to the 4th reflected wave) on lines 166-170.

In the end, two plots are created showing the acoustic pressure as a function of spatial Z versus spatial X. 

![Image of continuous plot](https://github.com/slkiser/acousticLevitation/blob/main/continuous.png)

![Image of contour plot](https://github.com/slkiser/acousticLevitation/blob/main/contour.png)

___

## Acknowledgments

This research was part of a PhD thesis funded by [Arts et Métiers](https://artsetmetiers.fr/) (École nationale supérieure d'arts et métiers). Laboratory equipment was provided by H2020 FastMat (fast determination of fatigue properties of materials beyond one billion cycles) project under the European Research Council (ERC) (grant agreement No 725142) at [PIMM laboratory](https://pimm.artsetmetiers.fr/).

![Image of logos](https://github.com/slkiser/lineSpectraVibration/blob/main/logo.png)

## License

Distributed under the MIT License. See [`LICENSE`](LICENSE) for more information.

<p align="right">(<a href="#top">back to top</a>)</p>
