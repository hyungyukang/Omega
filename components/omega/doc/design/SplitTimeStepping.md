# Omega1 Split Time Stepping Design

# 1 Overview

To enhance computational efficiency by allowing longer timesteps, ocean models require split barotropic–baroclinic time stepping methods. The implementation described here is based on the approach of Higdon (2005) and was employed in MPAS-Ocean, where the approach explicitly subcycles the barotropic terms. The split-explicit time stepping scheme in Omega1 is nearly identical to the `split_explicit` scheme in MPAS-Ocean, but in the tilted pressure $\(p^\*\)$ coordinate. Therefore, this document focuses specifically on split explicit timestepping in the $p^\*$ coordinate. The implicit treatment of the barotropic terms is planned for a future stage.

The split explicit method involves the following sequence:
 1) Decompose velocity into barotropic and baroclinic components.
 2) Advance the baroclinic velocities using a large timestep, and compute the vertically averaged forcing $\(\boldsymbol{\overline{G}})$
 3) Subcycle the barotropic velocity using small explicit timesteps.
 4) Recombine velocities and update other relevant variables.

    
# 2 Requirements

## 2.1 Requirement: A split time-stepping in the tilted pressure $\(p^*\)$ coordinate

The algorithm should be based on Section 2.3 of Higdon (2005), with modifications to accommodate the $p^\*$ coordinate variables. It should accept input parameters for the time-stepping scheme, the number of split-explicit iterations, the number of barotropic subcycles, and the number of baroclinic Coriolis iterations. An unsplit variant will also be included, which mirrors the split-explicit approach except that the full velocity is solved during the baroclinic stage, with no operations performed in the barotropic stage.

## 2.2 Requirement: Stable time integration for long-term and high-resolution simulations

Stability constrains the maximum allowable timestep, which in turn affects computational cost. The implemented split-explicit time stepping methods must allow for reasonably long timesteps while preventing numerical instabilities, such as those arising from internal gravity waves and barotropic modes, particularly important in global-scale and high-resolution ocean modeling. At a minimum, the time-stepping approach used in Omega1 should accommodate the same timestep sizes as MPAS-Ocean for both the baroclinic and barotropic subsystems.

## 2.3 Requirement: Modularization of baroclinic and barotropic time-stepping methods

Modularity ensures ease of testing and future-proofing of the Omega1 codebase. Implementing a modular design enables mix and match time-stepping schemes of the baroclinic and barotropic subsystems, straightforward integration of alternative time-stepping schemes, facilitates easier maintenance by separating the baroclinic and barotropic time stepping codes,

thereby enhancing innovation and flexibility.



# 3 Algorithmic Formulation

## 3.1 Split Explicit Scheme

The split explicit method separates ocean velocity into depth-integrated barotropic and depth-dependent baroclinic components. This separation allows computationally expensive baroclinic modes to run at longer timesteps and computationally efficient barotropic modes to run at shorter timesteps, enhancing computational efficiency and accuracy.

- **Barotropic Equations:**

$$
  \frac{\partial \zeta}{\partial t} + \nabla \cdot \left(u \sum_k h_k^{\text{edge}}\right) = 0
  \frac{\partial u}{\partial t} + fu^\perp = -g\nabla\zeta + G
$$

- **Baroclinic Equations:**

$$
  \frac{\partial u'_k}{\partial t} = -fu'^\perp_k + T(u_k, w_k, p_k) + g\nabla\zeta - G
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
\frac{\partial u}{\partial t} = -R u
\frac{\partial h}{\partial t} = 0
\frac{\partial h\phi}{\partial t} = \frac{h}{\tau}(\phi - \phi_0)
$$

These verification tests validate that the implemented schemes accurately represent known dynamics such as exponential decay processes and tracer restoration, common in ocean biogeochemical modeling.

