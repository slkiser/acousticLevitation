<div id="top"></div>

# Matrix method in Python for acoustic levitation simulations
This implements the matrix method used in the publication ["Matrix Method for Acoustic Levitation Simulation"](https://www.researchgate.net/publication/224254694_Matrix_Method_for_Acoustic_Levitation_Simulation). Specifically, the results of Section II are written in Python, which describes an acoustic levitator composed of only one circular flat transducer and one planar reflector.

This [README](README.md) consists of an introduction and theory which derives the **Rayleigh integral equation** used in the publication. This integral is then numerically resolved using the matrix method used in the publication. The matrix method offers an advantage over traditional [numerical integration](https://en.wikipedia.org/wiki/Numerical_integration) for two-dimensional (2D) problems, e.g., the [trapezoidal rule](https://en.wikipedia.org/wiki/Trapezoidal_rule) (1st order approximation) & [Simpson's rule](https://en.wikipedia.org/wiki/Simpson%27s_rule) (2nd order approximation), by obtaining a higher order of accuracy at a lower computation cost.

If you use this code, please cite this work:
> [1] S. L. Kiser, Matrix method in Python for acoustic levitation simulations. [Online]. Available: https://github.com/slkiser/acousticLevitation

## Introduction
Acoustics is a field that has ties to mechanics and fluid dynamics. For solid mechanics, a phenomenological approach to vibration is solids is adequate enough to describe modal shapes and eigenfrequencies. However, it is not adequate for relating mechanical vibrations to noise radiation and transmission. Therefore, the interaction of sound waves and the vibration of solids requires a fundamental wave approach. Consider the following:

- Solids can store potential energy in compression (longitudinal), flexural (bending/transverse), shear, and torsional waves.
- Fluids can only store energy in compression, thus they can only sustain compressional (longitudinal) waves.

We are interested in an approach is required that can relate the wave particle velocities of the vibrating solid surface to the volume velocities in the fluid. In other words, we seek an approach that relates the energy exchange between the solid and the fluid. 

At a fundamental level, the radiation of sound from a vibrating boundary (surface) can be formulated in terms of an integral equation involving [Green's functions](https://en.wikipedia.org/wiki/Green%27s_function) with imposed radiation constraints, i.e., the radiation condition ensures that the integral equation for the radiated sound pressure represents outward-traveling sound waves.

In its most general form, the integral equation is attributable to Kirchhoff, although Helmholtz modified it for single-frequency (harmonic) applications. This integral, [Kirchhoff-Helmholtz integral equation](https://en.wikipedia.org/wiki/Kirchhoff%E2%80%93Helmholtz_integral), relates harmonic vibrational motion of the surface to the radiated sound pressure field in the enclosed fluid. It can be interpreted as representing the sound pressure field of a vibrating surface by a distribution of volume velocity sources and forces on the surface. The velocity sources and the forces are related to normal surface velocity and surface pressure respectively. It is important to note that the surface pressure and the normal surface vibrational velocity interact with one another. 

Rayleigh modified the [Kirchhoff-Helmholtz integral equation](https://en.wikipedia.org/wiki/Kirchhoff%E2%80%93Helmholtz_integral) for the specific case of a planar source located in an infinite baffle and illustrated that it is equivalent to a distribution of point sources. This modification, named the  **Rayleigh integral equation** is outlined in the following subsection.

### Theory
The nature of **steady-state sound fields** in presence of boundaries (surfaces) is described by the [Helmholtz equation](https://en.wikipedia.org/wiki/Helmholtz_equation):

<div align=center> <img src="https://render.githubusercontent.com/render/math?math=%5Clarge+%5Cdisplaystyle+%5Cnabla%5E%7B2%7D+p%2Bk%5E%7B2%7D+p%3D0"> </div>

To describe the sound pressure  <img src=
"https://render.githubusercontent.com/render/math?math=%5Ctextstyle+p%28%5Cmathbf%7Br%7D%2C%5Comega%29" 
alt="p(\mathbf{r},\omega)"> in [polar coordinates](https://ocw.mit.edu/courses/aeronautics-and-astronautics/16-07-dynamics-fall-2009/lecture-notes/MIT16_07F09_Lec05.pdf) due to the boundary vibration and reflection, a [Green's function](https://en.wikipedia.org/wiki/Green%27s_function) is used to satisfy the following equation:

<div align=center> <img src=
"https://render.githubusercontent.com/render/math?math=%5Clarge+%5Cdisplaystyle+%5Cnabla%5E%7B2%7D+G%5Cleft%28%5Cmathbf%7Br%7D+%5Cmid+%5Cmathbf%7Br%7D%27%5Cright%29%2Bk%5E%7B2%7D+G%5Cleft%28%5Cmathbf%7Br%7D+%5Cmid+%5Cmathbf%7Br%7D%27%5Cright%29%3D-%5Cdelta%5Cleft%28%5Cmathbf%7Br%7D-%5Cmathbf%7Br%7D%27%5Cright%29" 
alt="\nabla^{2} G\left(\mathbf{r} \mid \mathbf{r}'\right)+k^{2} G\left(\mathbf{r} \mid \mathbf{r}'\right)=-\delta\left(\mathbf{r}-\mathbf{r}'\right)"></div>

Kirchhoff obtained the solution by substituting both previous equations for sound pressure in the integral form (the [Kirchhoff-Helmholtz integral](https://en.wikipedia.org/wiki/Kirchhoff%E2%80%93Helmholtz_integral)):

<div align=center>
<img src=
"https://render.githubusercontent.com/render/math?math=%5Clarge+%5Cdisplaystyle+p%28%5Cmathbf%7Br%7D%2C+%5Comega%29%3D%5Coint_%7Bs%7D%5Cleft%5BG%5Cleft%28%5Cmathbf%7Br%7D+%5Cmid+%5Cmathbf%7Br%7D%27%5Cright%29+%5Cfrac%7B%5Cpartial+p%7D%7B%5Cpartial+n%27%7D-p%5Cleft%28%5Cmathbf%7Br%7D%27%2C+%5Comega%5Cright%29+%5Cfrac%7B%5Cpartial+G%5Cleft%28%5Cmathbf%7Br%7D+%5Cmid+%5Cmathbf%7Br%7D%27%5Cright%29%7D%7B%5Cpartial+n%27%7D%5Cright%5D+%5Cmathrm%7Bd%7D+s%27" 
alt="p(\mathbf{r}, \omega)=\oint_{s}\left[G\left(\mathbf{r} \mid \mathbf{r}'\right) \frac{\partial p}{\partial n'}-p\left(\mathbf{r}', \omega\right) \frac{\partial G\left(\mathbf{r} \mid \mathbf{r}'\right)}{\partial n'}\right] \mathrm{d} s'">
</div>

where <img src="https://render.githubusercontent.com/render/math?math=%5Ctextstyle+%5Cpartial+p+%2F+%5Cpartial+n%27" alt="\partial p / \partial n'">represents the gradient on the boundary surface (outward positive from the sound generation).

The boundary conditions about the plate, <img src=
"https://render.githubusercontent.com/render/math?math=%5Ctextstyle+S_1" 
alt="S_1">, and the rigid reflector, <img src=
"https://render.githubusercontent.com/render/math?math=%5Ctextstyle+S_2" 
alt="S_2">, are:
<div align=center>
<img src=
"https://render.githubusercontent.com/render/math?math=%5Clarge+%5Cdisplaystyle+%5Cmathbf%7Bv%7D_%7BB%7D%3D%5Cleft%5C%7B%5Cbegin%7Barray%7D%7Bcl%7D%0A%5Cmathbf%7Bv%7D_%7Bn%7D%5Cleft%28%5Cmathbf%7Br%7D_%7Bp%7D%2C+%5Comega%5Cright%29+%5Cmathrm%7Be%7D%5E%7Bj+%5Comega+t%7D+%26+%5Cmathbf%7Br%7D_%7Bp%7D%3D%5Cleft%28x%27%2C+y%27%5Cright%29+%5Cin+S_%7B1%7D+%5C%5C%0A0+%26+%5Cmathbf%7Br%7D_%7Bp%7D%3D%5Cleft%28x%27%2C+y%27%5Cright%29+%5Cin+S_%7B2%7D%0A%5Cend%7Barray%7D%5Cright." 
alt="\mathbf{v}_{B}=\left\{\begin{array}{cl}
\mathbf{v}_{n}\left(\mathbf{r}_{p}, \omega\right) \mathrm{e}^{j \omega t} & \mathbf{r}_{p}=\left(x', y'\right) \in S_{1} \\
0 & \mathbf{r}_{p}=\left(x', y'\right) \in S_{2}
\end{array}\right.">
</div>

where the polar velocity is <img src="https://render.githubusercontent.com/render/math?math=%5Ctextstyle+%5Cmathbf%7Bv%7D%3Dd%5Cmathbf%7Br%7D%2Fdt" alt="\mathbf{v}=d\mathbf{r}/dt">.

By applying [Euler's equation](https://en.wikipedia.org/wiki/Euler_equations_(fluid_dynamics)), or a **balance of momentum** ([Newton's second law](https://en.wikipedia.org/wiki/Newton%27s_laws_of_motion#Newton's_second_law)), the pressure gradient on the boundary surface is equal to the boundary vibration velocity:

<div align=center>
<img src=
"https://render.githubusercontent.com/render/math?math=%5Clarge+%5Cdisplaystyle+%5Cfrac%7B%5Cpartial+p%7D%7B%5Cpartial+n%27%7D%3Dj+%5Crho_%7B0%7D+%5Comega+%5Cmathbf%7Bv%7D_%7BB%7D" 
alt="\frac{\partial p}{\partial n'}=j \rho_{0} \omega \mathbf{v}_{B}">
</div>

For this boundary, a [Green's function](https://en.wikipedia.org/wiki/Green%27s_function) is needed to satisfy:

<div align=center>
<img src=
"https://render.githubusercontent.com/render/math?math=%5Clarge+%5Cdisplaystyle+%5Cfrac%7B%5Cpartial+G%5Cleft%28%5Cmathbf%7Br%7D+%5Cmid+%5Cmathbf%7Br%7D%27%5Cright%29%7D%7B%5Cpartial+n%27%7D%3D0" 
alt="\frac{\partial G\left(\mathbf{r} \mid \mathbf{r}'\right)}{\partial n'}=0">
</div>

at the boundary surface (when <img src="https://render.githubusercontent.com/render/math?math=%5Ctextstyle+z%3D0" alt="z=0"> ). The [Green's function](https://en.wikipedia.org/wiki/Green%27s_function) selected by Rayleigh is:

<div align=center>
<img src=
"https://render.githubusercontent.com/render/math?math=%5Clarge+%5Cdisplaystyle+G%5Cleft%28%5Cmathbf%7Br%7D+%5Cmid+%5Cmathbf%7Br%7D%5E%7B%5Cprime%7D%5Cright%29%3D%5Cfrac%7B1%7D%7B4+%5Cpi%7D%5Cleft%28%5Cfrac%7B%5Cmathrm%7Be%7D%5E%7B-j+k+R_%7B%2B%7D%7D%7D%7BR_%7B%2B%7D%7D%2B%5Cfrac%7B%5Cmathrm%7Be%7D%5E%7B-j+k+R_%7B-%7D%7D%7D%7BR_%7B-%7D%7D%5Cright%29" 
alt="G\left(\mathbf{r} \mid \mathbf{r}^{\prime}\right)=\frac{1}{4 \pi}\left(\frac{\mathrm{e}^{-j k R_{+}}}{R_{+}}+\frac{\mathrm{e}^{-j k R_{-}}}{R_{-}}\right)">
</div>
where:

<div align=center>
<img src=
"https://render.githubusercontent.com/render/math?math=%5Clarge+%5Cdisplaystyle+%5Cbegin%7Barray%7D%7Bl%7D%0AR_%7B%2B%7D%3D%5Csqrt%7B%5Cleft%28x-x%27%5Cright%29%5E%7B2%7D%2B%5Cleft%28y-y%27%5Cright%29%5E%7B2%7D%2B%5Cleft%28z-z%27_%7B%2B%7D%5Cright%29%5E%7B2%7D%7D+%5C%5C%0AR_%7B-%7D%3D%5Csqrt%7B%5Cleft%28x-x%27%5Cright%29%5E%7B2%7D%2B%5Cleft%28y-y%27%5Cright%29%5E%7B2%7D%2B%5Cleft%28z%2Bz%27_%7B-%7D%5Cright%29%5E%7B2%7D%7D%0A%5Cend%7Barray%7D" 
alt="\begin{array}{l}
R_{+}=\sqrt{\left(x-x'\right)^{2}+\left(y-y'\right)^{2}+\left(z-z'_{+}\right)^{2}} \\
R_{-}=\sqrt{\left(x-x'\right)^{2}+\left(y-y'\right)^{2}+\left(z+z'_{-}\right)^{2}}
\end{array}">
</div>

with <img src="https://render.githubusercontent.com/render/math?math=%5Ctextstyle+z%27_%7B%2B%7D+" alt="z'_{+} "> and <img src="https://render.githubusercontent.com/render/math?math=%5Ctextstyle++z%27_%7B-%7D" alt=" z'_{-}"> being small positive and negative perturbations in the Z direction. This [Green's function](https://en.wikipedia.org/wiki/Green%27s_function) has a zero gradient on the boundary of <img src="https://render.githubusercontent.com/render/math?math=%5Ctextstyle+z%3D0" alt="z=0">:

<div align=center>
<img src="https://render.githubusercontent.com/render/math?math=%5Clarge+%5Cdisplaystyle+%5Cleft.%5Cfrac%7B%5Cpartial+G%5Cleft%28%5Cmathbf%7Br%7D+%5Cmid+%5Cmathbf%7Br%7D%27%5Cright%29%7D%7B%5Cpartial+z%27%7D%5Cright%7C_%7B%5C%7Bz%27_%7B%2B%7D%2C+z%27_%7B-%7D%5C%7D%5Crightarrow+0%7D%3D0" alt="\left.\frac{\partial G\left(\mathbf{r} \mid \mathbf{r}'\right)}{\partial z'}\right|_{\{z'_{+}, z'_{-}\}\rightarrow 0}=0"></div>

and when the sources are limited to the boundary surface <img src= "https://render.githubusercontent.com/render/math?math=%5Ctextstyle+%5Cmathbf%7Br%7D%27%5Crightarrow+%5Cmathbf%7Br%7D_%7Bp%7D" alt="\mathbf{r}'\rightarrow \mathbf{r}_{p}">, Rayleigh gets the second Green's function:

<div align=center>
<img src=
"https://render.githubusercontent.com/render/math?math=%5Clarge+%5Cdisplaystyle+G%5Cleft%28%5Cmathbf%7Br%7D+%5Cmid+%5Cmathbf%7Br%7D_%7Bp%7D%5Cright%29%3D%5Cfrac%7B1%7D%7B2+%5Cpi%7D+%5Cfrac%7B%5Cmathrm%7Be%7D%5E%7B-j+k+R%7D%7D%7BR%7D" 
alt="G\left(\mathbf{r} \mid \mathbf{r}_{p}\right)=\frac{1}{2 \pi} \frac{\mathrm{e}^{-j k R}}{R}">
</div>

where:

<div align=center>
<img src=
"https://render.githubusercontent.com/render/math?math=%5Clarge+%5Cdisplaystyle+R%3D%5Cleft%7C%5Cmathbf%7Br%7D-%5Cmathbf%7Br%7D_%7Bp%7D%5Cright%7C%3D%5Csqrt%7B%5Cleft%28x-x%27%5Cright%29%5E%7B2%7D%2B%5Cleft%28y-y%27%5Cright%29%5E%7B2%7D%2Bz%5E%7B2%7D%7D" 
alt="R=\left|\mathbf{r}-\mathbf{r}_{p}\right|=\sqrt{\left(x-x'\right)^{2}+\left(y-y'\right)^{2}+z^{2}}">
</div>

Substituting the past 3 equations into the [Kirchhoff-Helmholtz integral](https://en.wikipedia.org/wiki/Kirchhoff%E2%80%93Helmholtz_integral) results in the **Rayleigh integral** for the sound pressure field of a vibration source and a planar reflector:
<div align=center>
<img src=
"https://render.githubusercontent.com/render/math?math=%5Clarge+%5Cdisplaystyle+%5Cboxed%7Bp%28%5Cmathbf%7Br%7D%2C+%5Comega%29%3D%5Cfrac%7Bj+%5Crho_%7B0%7D+%5Comega%7D%7B2+%5Cpi%7D+%5Cint_%7BS_%7B1%7D%7D+%5Cfrac%7Bv_%7Bn%7D%5Cleft%28%5Cmathbf%7Br%7D_%7Bp%7D%2C+%5Comega%5Cright%29%7D%7BR%7D+%5Cmathrm%7Be%7D%5E%7B-j+k+R%7D+%5Cmathrm%7B%7Ed%7D+s%7D" 
alt="\boxed{p(\mathbf{r}, \omega)=\frac{j \rho_{0} \omega}{2 \pi} \int_{S_{1}} \frac{v_{n}\left(\mathbf{r}_{p}, \omega\right)}{R} \mathrm{e}^{-j k R} \mathrm{~d} s}">
</div>

The **Rayleigh integral** is then solved using the matrix method outlined in the publication ["Matrix Method for Acoustic Levitation Simulation"](https://www.researchgate.net/publication/224254694_Matrix_Method_for_Acoustic_Levitation_Simulation). Numerically, the **Rayleigh integral** is determined by discretizing the transducer and the reflector surfaces in small area cells, and multiplication of the transfer matrices gives the name of the matrix method.


## Setup and usage
We recommend the usage of the [Anaconda distribution](https://www.anaconda.com/products/individual) and the included [Spyder IDE](https://www.spyder-ide.org/). Spyder is a unique Python IDE that offers a [MATLAB](https://www.mathworks.com/products/matlab.html)-like experience with an integrated variable viewer.

### Package installation
For required packages, refer to [`requirements.txt`](requirements.txt). The packages can be installed quickly using pip and/or in a virtual environment:
```
pip install -r requirements.txt
```
Or if using Anaconda:
```
conda install pip
pip install -r requirements.txt
```
