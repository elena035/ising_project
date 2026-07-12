# 2D Ising Model Simulation

This repository contains a parallelized simulation of the 2D Ising Model performed for different lattice sizes ($L \times L$), specifically $L = 16, 32, 64,$ and $128$.

The simulation aims to study the phase transition occurring approximately at the critical inverse temperature $\beta_c \approx 0.4406$. This is achieved by computing the intensive instantaneous magnetization, defined as:

$$ m = \frac{1}{L^2} \sum_{i=1}^{L^2} s_i$$

where $s_i \in \{-1, +1\}$ is the value of the $i$-th spin.


The simulation code is written in C (parallelized with MPI) and is fully available in the parallel_code directory. The program implements the Metropolis algorithm to sample the model's configurations for fixed values of $\beta$ and $L$.

## Metropolis algorithm

The algorithm evolves the system through the following steps:

1. Initialization: Start from a random configuration for the $L \times L$ lattice, where each spin is randomly assigned a value of $+1$ or $-1$.

2. Proposal: Propose a new energy configuration by picking a random spin and flipping its sign (single spin-flip).

3. Energy Evaluation: Compute the energy difference between the new and old configurations:

$$ \Delta E = E_{new} - E_{old} = 2 s_i \sum_{\langle i,j \rangle} s_j$$

where $s_i$ is the value of the flipped spin (before the flip) and the sum runs over its 4 nearest neighbors. Since the lattice has a coordination number of 4, $\Delta E$ is a discrete variable that can only assume 5 specific values: $\{-8J, -4J, 0, +4J, +8J\}$.

4. Acceptance Step: According to the Metropolis criterion, the probability of accepting the proposed configuration is given by:

$$  p= \min(1, e^{-\beta \Delta E}) $$


- If $\Delta E \le 0$, the algorithm always accepts the move.

- If $\Delta E > 0$, the exponential $e^{-\beta \Delta E}$ is calculated and compared with a random number $r \in [0, 1)$. The new configuration is accepted only if $r < e^{-\beta \Delta E}$.

5. Iteration: Repeat these steps to evolve the system over time.


## Termalization and measurement

To avoid initial configuration bias (which would corrupt the macroscopic measurements), the measurement phase begins only after an initial thermalization phase.

The system is updated with $N_{eq}$ full iterations (sweeps) of the Metropolis algorithm for each fixed value of $\beta$. This process erases the initialization memory of the system, allowing it to reach dynamic equilibrium for that specific temperature.

Once equilibrium is reached, the measurement phase begins. To obtain solid and reliable statistics, physical observables are recorded only every $N_{skip}$ iterations. This spacing guarantees that the measurements are not autocorrelated, ensuring that all sampled configurations can be considered independent statistical realizations of the same thermodynamic process.

## Repository Structure:

/parallel_code/ : Contains the C source code and the Makefile.

/scripts/ : Bash scripts for production runs

/input_parameters/ : Contains the input file /.txt/ used in the simulation

/ising_analysis/ : Notebook with python scripts for results visualization and results' analysis

/assets/ : Generated plots, images, and GIFs of the lattice evolution.