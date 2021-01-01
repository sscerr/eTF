# README #

Here the 2.5D version of the eTF code to evolve the "extended Two-Fluid" plasma model equations 

The eTF model is presented in: Cerri et al., Phys. Plasmas 20, 112112 (2013).
Refer to this paper if using this code.

See eTFdocumentation.pdf for installation and compilation guidelines.

### Version 1.0 ###

* parallel (MPI) code
* "2.5D" code: two-dimensional (x,y) simulation domani (no variation along z, d/dz = 0), fully 3D vectors (e.g., B = (Bx, By, Bz))  
* simulation boundaries: open conditions along x, periodic conditions along y  

* generalized Ohm’s law with: Hall term (JxB), thermo-electric effect (div[P_e]) 
* anisotropic pressures for both species (protons and electrons)  
* 1st-order finite-Larmor-radius (FLR) corrections of the ions (assumes B mainly along z) 
* double-adiabatic closure (zero heat fluxes) 
* massless electrons (no electron-inertia effects)  

* default initial condition: FLR-corrected shear-flow layer (e.g., Kelvin-Helmholtz instability)

