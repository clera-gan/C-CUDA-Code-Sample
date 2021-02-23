# Sample C++ codes 
C++/CUDA C codes for GPU-based DEM for a ship loading process developed by Dr. Gan using C/C++ for CPU part and CUDA C for GPU part. GPU codes are not included in the repo. Please refer to papers[1-2] for more details:
[1].	J. Gan, T. Evans, A. Yu, Impact energy dissipation study in a simulated ship loading process, Powder Technology, 354 (2019) 476-484.
[2].	J. Gan, T. Evans, A. Yu, Application of GPU-DEM simulation on large-scale granular handling and processing in ironmaking related industries, Powder Technology, 361 (2020) 258-273.

Code composition: 
1) input folder
2) src folder for source codes
3) output folder for simulation results
4) Executive code --“GDEM” file, to run the simulation
5) job script --“case” file, for job management on linux system

input:
feed.dat to set feeding stream/mode/mass flowrate/feed material
hopp3d.inp: general information
material.dat:  to set material properties
meshinfo.dat: mesh number/name, offset, scale factor
movements.dat to set the movement of the meshes (from mesh file folder), rotation/translation/vibration
restart.dat file to  continue to run the simulation from certain time
output:  time/statistic results of simulations

src:
source code: to generate execute file, input make to compile the code
AllocateArrays.cpp: to allocate cpu, gpu arrays
Boundary.cpp : to set boundary functions
Feed.cpp : to set feed functions
FreeArrays.cpp : free cpu, gpu arrays
GPUSet.cpp: set GPU functions
Initialization.cpp: initialize parameters and arrays
Materials.cpp: to set material related functions
MemcpyHostDevice.cpp : to realized data copy between cpu and gpu
Movement.cpp  : set movement functions
mpiFunctions.cpp : to set mpi functions
Particle.cpp  : to set particle related functions, suggested to merge into other functions
ReadData.cpp : read input data
WriteData.cpp : to output data
dempacking.cpp : main function to run the DEM simulations including CPU and GPU functions.
Makefile: to compile the *.cpp/*.cu codes, and generate the target file (e.g. objective file /dynamic library)
