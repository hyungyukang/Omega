# Omega1 Split Time Stepping Design

# 1 Overview

This document describes the time-stepping methods used in Omega1, particularly emphasizing the split explicit scheme adapted for p-level vertical coordinates, following the methods described by Higdon (2005). The primary objective is to achieve stable, efficient, and accurate numerical integration of ocean circulation model equations.

# 2 Requirements

## 2.1 Accurate Discretization

Temporal discretization in ocean modeling is crucial because inaccuracies can accumulate and degrade long-term simulations. Therefore, Omega1 must maintain at least second-order accuracy in time to ensure high fidelity in representing ocean dynamics, especially when simulating long-term climate scenarios or capturing fine-scale oceanographic processes.

## 2.2 Stability

Stability in numerical ocean modeling dictates the largest feasible timestep, directly influencing computational cost. The implemented time-stepping methods must allow for reasonably long timesteps while preventing numerical instabilities, such as those arising from internal gravity waves and barotropic modes, particularly important in global-scale and high-resolution ocean modeling.

## 2.3 Performance

Ocean models like Omega1 are typically run on high-performance computing (HPC) platforms. Hence, the implemented scheme must efficiently scale across:

- Single CPU performance to enable rapid development and debugging.
- Parallel CPU architectures to utilize multiple nodes efficiently, reducing wall-clock time for operational and research simulations.
- GPU acceleration for massively parallel architectures, significantly enhancing computational speed and efficiency.

## 2.4 Modularity

Modularity ensures ease of testing and future-proofing of the Omega1 codebase. Implementing a modular design enables straightforward integration of alternative time-stepping schemes, facilitates easier maintenance, and encourages contributions from a broader scientific community, thereby enhancing innovation and flexibility.

## 2.5 Conservation

Ensuring conservation of physical properties like mass and tracers is essential for the physical reliability of ocean model simulations. Omega1 must rigorously conserve total mass (volume-integrated layer thickness, $\(h\)$) and total tracers (area-integrated $\(\phi\)$), thereby ensuring accurate long-term climate simulations and realistic modeling of biogeochemical processes.

## 2.6 Explicit Time Argument for RHS

Explicitly including a time argument in the RHS functions is necessary to accurately handle time-dependent forcing terms (e.g., tides, seasonal wind stresses, and surface buoyancy fluxes). This allows precise alignment of model forcing with observational datasets or climate model outputs, significantly improving model realism and predictability.

# 3 Algorithmic Formulation

## 3.1 Split Explicit Scheme

The split explicit method separates ocean velocity into depth-integrated barotropic and depth-dependent baroclinic components. This separation allows computationally expensive baroclinic modes to run at longer timesteps and computationally efficient barotropic modes to run at shorter timesteps, enhancing computational efficiency and accuracy.

- **Barotropic Equations:**
$$
  \[\frac{\partial \zeta}{\partial t} + \nabla \cdot \left(u \sum_k h_k^{\text{edge}}\right) = 0\]
  \[\frac{\partial u}{\partial t} + fu^\perp = -g\nabla\zeta + G\]
$$

- **Baroclinic Equations:**
$$
  \[\frac{\partial u'_k}{\partial t} = -fu'^\perp_k + T(u_k, w_k, p_k) + g\nabla\zeta - G\]
$$

## 3.2 Runge-Kutta 4th Order (RK4) Scheme

The RK4 scheme provides a high-order, explicit alternative useful for benchmarking and validating simpler schemes. Due to its computational cost, RK4 typically serves as a reference solution to assess the accuracy and convergence of other numerical methods in idealized test cases.

## 3.3 Unsplit Explicit Scheme

The unsplit explicit scheme directly solves the ocean model equations without decomposing the velocity field, providing a straightforward method ideal for benchmarking, verification, and simpler configurations or coarse-resolution experiments.

# 4 Design and Implementation

## 4.1 Implementation Stages

### Stage 1: Baroclinic (3D)
- Integrate baroclinic velocities over large timesteps, efficiently resolving internal ocean dynamics.
- Calculate the barotropic forcing terms based on updated baroclinic velocities, ensuring consistency between model components.

### Stage 2: Barotropic (2D)
- Explicitly subcycle the faster barotropic velocities and sea surface height, ensuring numerical stability via a predictor-corrector approach.
- Employ multiple smaller timesteps for the barotropic mode, capturing rapid changes such as surface gravity waves and ensuring dynamic stability.

### Stage 3: Tracers, Thickness, and Diagnostics
- Update passive tracers, oceanic layer thickness, density fields, and pressure diagnostics using velocities from previous stages, ensuring accurate representation of ocean transport and mixing processes.

## 4.2 Configuration Parameters
- `config_time_integration`: (`'RK4'`, `'split_explicit'`, `'unsplit'`)
- `config_n_ts_iter`: Number of baroclinic iterations per timestep (default: 2)
- `config_n_btr_subcycles`: Number of barotropic subcycles per baroclinic timestep
- `config_compute_tr_midstage`: Toggle for computing tracer updates mid-stage

## 4.3 Variables
- `uBtr`: Barotropic horizontal velocity
- `ssh`: Sea surface height
- `uBtrSubcycle`: Temporally subcycled barotropic velocity
- `sshSubcycle`: Temporally subcycled sea surface height
- `FBtr`: Barotropic mass flux
- `GBtrForcing`: Forcing term derived from baroclinic mode
- `uBcl`: Baroclinic horizontal velocity

# 5 Verification and Testing

Comprehensive testing ensures robustness, accuracy, and scalability:

- **Time-only convergence tests** to evaluate temporal accuracy against analytical solutions.
- **Comparisons among RK4, unsplit, and split explicit methods** to benchmark accuracy and stability.
- **Validation of conservation properties** to ensure numerical solutions adhere to physical laws of mass and tracer conservation.
- **Performance assessment** across diverse computational architectures (single CPU, parallel CPUs, and GPUs).

Verification equations:

$$
\[\frac{\partial u}{\partial t} = -R u\]
\[\frac{\partial h}{\partial t} = 0\]
\[\frac{\partial h\phi}{\partial t} = \frac{h}{\tau}(\phi - \phi_0)\]
$$

These verification tests validate that the implemented schemes accurately represent known dynamics such as exponential decay processes and tracer restoration, common in ocean biogeochemical modeling.

# 6 Summary

This enhanced design document details a structured approach for robust, efficient, and accurate ocean model simulations, ensuring Omega1's continued development as a leading ocean modeling tool for climate research and operational forecasting.

