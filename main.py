"""
%=============================================================================
% DESCRIPTION:
% This is the Python implementation of the implements the matrix method for an
% acoustic levitator composed of only one circular flat transducer and one 
% planar reflector.
%
%=============================================================================
% Version 1.0, Authored by:
% Shawn L. KISER (Msc) @ https://www.linkedin.com/in/shawn-kiser/
%   Laboratoire PIMM, Arts et Metiers Institute of Technology, CNRS, Cnam,
%   HESAM Universite, 151 boulevard de l’Hopital, 75013 Paris (France)
%
% Based on:
% [1] M. A. B. Andrade, N. Perez, F. Buiochi, and J. C. Adamowski, “Matrix 
%     method for acoustic levitation simulation,” IEEE Trans. Ultrason., 
%     Ferroelect., Freq. Contr., vol. 58, no. 8, pp. 1674–1683, Aug. 2011, 
%     doi: 10.1109/TUFFC.2011.199.

%
%=============================================================================
% The MIT License (MIT)
% 
% Copyright (c) 2021 Shawn L. KISER
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.

%=============================================================================
"""
#%% Imports
import numpy as np
import matplotlib.pylab as plt
mm = 1e-3

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

#%% Creation of grid, indices, and distance arrays ((r) in paper)
x = np.linspace(x_min, x_max, int((x_max-x_min)/disc))
z = np.linspace(z_min, z_max, int((z_max-z_min)/disc))
transducer = np.linspace(-R, R, int((2*R)/disc))

# Array indices
N = int((2*R)/disc)
M = int(len(x)*len(z))
L = int(len(x))
Z = int(len(z))

# Create zeroed distance arrays in memory
r_nm = np.zeros((N, M))
r_im = np.zeros((L, M))
r_in = np.zeros((N, L))

# Calculate distance array r_nm
for i in range(N):
    q = 0
    for j in range(L):
        for k in range(len(z)):
            r_nm[i, q] = np.sqrt((transducer[i]-x[j])**2+(z_max-z[k])**2)
            q = q+1

# Calculate distance array r_im
for i in range(L):
    q = 0
    for j in range(L):
        for k in range(len(z)):
            r_im[i, q] = np.sqrt((x[i]-x[j])**2+(z_min-z[k])**2)
            q = q + 1

# Calculate distance array r_in
for i in range(N):
    for j in range(L):
        r_in[i, j] = np.sqrt((transducer[i]-x[j])**2+(z_max-z_min)**2)
r_ni = r_in.T

#%% Creation of cell discretization and transfer matrices calculations
cells = complex(100, 0)			        # Number of discrete cells
sn = A_trans/cells			            # Unit cell area for transducer
si = A_refl/(cells*4)                   # Unit cell area for reflector
sh = A_hole/cells                       # Unit cell area for transducer hole

# Create zeroed transfer matrices in memory
T_TM = np.matrix(np.zeros((N, M)), dtype=complex)
T_TR = np.matrix(np.zeros((N, L)), dtype=complex)
T_RT = np.matrix(np.zeros((L, N)), dtype=complex)
T_RM = np.matrix(np.zeros((L, M)), dtype=complex)

# Calculate transfer matrices
for i in range(N):
    for j in range(M):
        T_TM[i, j] = ((sn*np.exp(complex(0, -1) *
                                wavenumber*r_nm[i, j]))/r_nm[i, j]) - \
                     ((sh*np.exp(complex(0, -1) *
                                wavenumber*r_nm[i, j]))/r_nm[i, j])
    for k in range(L):
        T_TR[i, k] = ((sn*np.exp(complex(0, -1) *
                                wavenumber*r_in[i, k]))/r_in[i, k]) - \
                     ((sh*np.exp(complex(0, -1) *
                                wavenumber*r_in[i, k]))/r_in[i, k])

for i in range(L):
    for j in range(N):
        T_RT[i, j] = (si*np.exp(complex(0, -1) *
                                wavenumber*r_ni[i, j]))/r_ni[i, j]

for i in range(L):
    for j in range(M):
        if r_im[i, j] == 0.0:
            continue
        T_RM[i, j] = (si*np.exp(complex(0, -1) * wavenumber*r_im[i, j]))/r_im[i, j]

# Rotation for column notation
T_TM = T_TM.T
T_TR = T_TR.T
T_RT = T_RT.T
T_RM = T_RM.T

#%% Boundary conditions of transducer
U = np.array(np.zeros((N, 1)), dtype=complex)
for i in range(N):
    U[i] = U_0*np.exp(complex(0, omega))

#%% Calculation of pressure, where each line is an order of approximation
P = D * T_TM * U + \
    (D*E) * (T_RM*T_TR) * U + \
        (D*E**2) * (T_TM*T_RT) * (T_TR*U) + \
            (D*E**3) * T_RM * (T_TR*T_RT) * (T_TR*U)+\
                (D*E**4)*(T_TM*T_RT)*(T_TR*T_RT)*(T_TR*U)
    
# Reshape of numpy array with vertical Z and horizontal X
P2 = np.array(np.reshape(P, (len(x), len(z))))
P2 = np.transpose(P2)

#%% Plots
# Plot contour plots
levels = np.linspace(P2.real.min()/2, -P2.real.min()/2, 10)
cntr  = plt.contourf(x,z,P2.real,levels,cmap=plt.cm.jet,extend="both")
plt.colorbar().set_label('Acoustic pressure (Pa)')
plt.contour(x,z,P2.real,levels,colors='black', linewidths=0.5, linestyles='solid')
plt.title('Pressure distribution, 10 contours')
plt.xlabel('X (m)')
plt.ylabel('Z (m)')
plt.savefig('contour.png',dpi=300)
plt.show()

# Plot continuous plots
plt.pcolormesh(x, z, P2.real, cmap=plt.cm.jet, shading='gouraud')
plt.title('Pressure distribution, continuous')
plt.xlabel('X (m)')
plt.ylabel('Z (m)')
plt.clim(min(P.real)/2, -min(P.real)/2)
plt.colorbar().set_label('Acoustic pressure (Pa)')
plt.savefig('continuous.png',dpi=300)
plt.show()