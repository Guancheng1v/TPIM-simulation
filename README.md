# Temporal Photonic Ising Machine (TPIM) Simulation
This repository contains the MATLAB simulation code and implementation scripts for the paper: "Temporal photonic Ising machine based on time-wavelength multiplexing and dispersive Fourier transform".
The code simulates a temporal photonic Ising machine (TPIM) that maps Ising spins to optical pulses modulated from an optical frequency comb and stretched by a chirped fiber Bragg grating (CFBG). It includes the numerical physical simulation of pulse propagation, the implementation of the stream-style annealing algorithm, and visualization tools to produce the main figures in the manuscript.
## **Repository Structure**
The files in this repository are organized as follows:
### **1. Physical Simulation Core**
`TPIM_simulation.m`: A MATLAB class definition that models the optical components, including the optical frequency comb, amplitude/phase modulation via electro-optic modulators (EOM), and dispersion propagation through the CFBG.

`JtoKsi.m`: Implements the Cholesky-like decomposition algorithm to decompose a interaction matrix $J$ into Decoupling vectors $\xi$ (Mattis-type spin interactions).
### **2. Main Algorithms & Solvers**
`weight_3N.m`: The main physical simulation script. It performs the stream-style optical propagation and intensity-based Hamiltonian evaluation.

`verify_softmax_annealing.m`: An algorithmic comparison script that compares the proposed algorithm against the classical single-spin-flip Simulated Annealing (SA) benchmark.

`HIcheck.m`: Verifies the linear relationship between the computed physical optical intensity ($I$) and the theoretical mathematical Hamiltonian ($H$).

### **3. Problem Generator Scripts**
`generate_mattis_matrix.m`: Generates specific rank Mattis interaction matrices with discrete values {−1,0,1}.

`generatemaxcut.m`: Generates standard benchmark instances for the 100-spin Max-Cut problem: a Möbius ladder graph, an unweighted random graph, and a frustrated random graph.

### **4. Visualization & Plotting Scripts**
`plotcircle.m`: Plots the circular spin-interaction network diagram (as seen in insets of Fig.3 and Fig. 5a,d,g).

`plotrank1.m` & `plotrankmorethan1.m`: Visualize the decomposed vectors for different matrix ranks (as seen in Fig.4a,d,g).

`plotruns.m`: Plots multiple annealing trajectories showing energy decay alongside the temperature schedule (as seen in Fig.4c,f,i).

`plotmaxcutiteration.m`: Plots the convergence curve showing both the Ising Energy and the Max-Cut size (as seen in Fig.5c,f,i).

`plothistoIsing.m` & `plothistomaxcut.m`: Generate histograms of the final ground state success distributions comparing TPIM and SA (as seen in Fig4.b,e,h and Fig5.b,e,h).

## Requirements and Environment
The code was developed and tested under the following environment:
Software: MATLAB (R2021a or later recommended).
Hardware: A CUDA-enabled GPU is recommended, as several core physical simulation steps utilize gpuArray for accelerated parallel operations.
Note: If a compatible GPU is not available, you can modify gpuArray calls in the scripts to standard CPU arrays, though execution time may increase.

## How to Run the Simulations
### 1: Verify the Physical Model (H-I Calibration)
To verify that the simulated optical intensity from the TS-DFT matches the mathematical Ising energy:
Ensure you have an interaction matrix loaded (e.g., J11 or generate one using generate_mattis_matrix).
Run `HIcheck.m`. This will output a scatter plot correlating physical intensity with the mathematical Hamiltonian.

### 2: Compare Algorithms (Pure Algorithmic Level)
To test the convergence behavior of the stream-style algorithm against classical SA without simulating the full physical pulse propagation:
Run `verify_softmax_annealing.m`.
This script evaluates a selected benchmark problem (e.g., Möbius ladder).

### 3: Run the Complete Physical Optoelectronic Solver
To execute the physical-simulation-based simulated annealing loop:
Run `weight_3N.m`.
This script simulates the physical EOM modulation, fiber dispersion, and intensity sampling within the iterative annealing loop.

## Note on Code Comments
Some inline comments and documentation within the .m files are written in Chinese. However, the variables, function names, and execution flow are written in standard English and align directly with the formulas presented in the manuscript.
