# photoNs-2.0
A simplified version of Cosmological Simulation Code for flat Lambda-CDM models

The first version of photoNs [Q. Wang, et al. RAA, (2018) 18, 62] is a TreePM code designed for massively parallel simulations.

In the second veriosn, we combine the Particle-Mesh (PM) and Fast Multipole Method (FMM) to compute the gravitational interactions. PM method needs a FFT lib of two-dimensional pencil decomposition [http://www.2decomp.org/] and FMM is carried out in Cartesian coordinates [e.g, Dehnen, W. Comput. Astrophys. (2014) 1: 1.] with a new domain decomposition, MPI communications, acceptance crietiria, etc. It contains a basic implementation but main functions for a flat Lamdba-CDM universe, as a benchmark for optimizations and algorithm testification.

The initial condition in Gadget2 format and parameter files for demonstrations can be found in the folder demo/ and 'make demo' will output a snapshot into the default folder run/



