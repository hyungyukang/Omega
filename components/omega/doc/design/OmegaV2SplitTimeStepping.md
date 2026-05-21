(omega-v2-split-time-stepping)=
# Omega V2: Split Time Stepping

<!--
Add 4.design
Add 5.tests
Add SE-AB2
-->

<!--- use table of contents if desired for longer documents  -->
**Table of Contents**
1. [Overview](#1-overview)
2. [Requirements](#2-requirements)
3. [Algorithmic Formulation](#3-algorithmic-formulation)
4. [Design](#4-design)
5. [Verification and Testing](#5-verification-and-testing)


## 1. Overview

To enhance computational efficiency by allowing longer timesteps, ocean models require split barotropic-baroclinic time stepping methods. The implementation described here is based on the approach of [Higdon (2005)](https://www.sciencedirect.com/science/article/pii/S0021999104005236) and was employed in MPAS-Ocean, where the approach explicitly subcycles the barotropic terms. The split-explicit time stepping scheme in Omega V2 is nearly identical to the `split_explicit` scheme in MPAS-Ocean, but in the non-Boussinesq pseudo-height $\tilde{z}$ vertical coordinate. Therefore, this document focuses specifically on split-explicit time stepping in the $\tilde{z}$ coordinate. The implicit treatment of the barotropic terms is planned for a future stage.

The split-explicit method involves the following sequence:

- Decompose velocity into barotropic and baroclinic components.
- Advance the baroclinic velocities using a large timestep, and compute the vertically averaged forcing, the $\overline{G}$ term.
- Subcycle the barotropic velocity using small explicit timesteps.
- Recombine velocities and update other relevant variables.

## 2. Requirements

### 2.1 Requirement: A split time-stepping method in the pseudo-height coordinate

The algorithm is based on Section 2.3 of [Higdon (2005)](https://www.sciencedirect.com/science/article/pii/S0021999104005236), with modifications to accommodate the $\tilde{z}$-coordinate variables in the non-Boussinesq framework. It should accept input parameters for the time-stepping scheme, the number of split-explicit iterations, the number of barotropic subcycles, and the number of baroclinic Coriolis iterations. An unsplit variant will also be included, which mirrors the split-explicit approach except that the full velocity is solved during the baroclinic stage, with no operations performed in the barotropic stage.

### 2.2 Requirement: Stable time integration for long-term high-resolution simulations

Stability constrains the maximum allowable timestep, which in turn affects computational cost. The implemented split-explicit time stepping methods must allow for reasonably long timesteps while preventing numerical instabilities, such as those arising from internal gravity waves and barotropic modes, which are particularly important in global-scale and high-resolution ocean modeling. At a minimum, the time-stepping approach used in Omega V2 should accommodate the same timestep sizes as MPAS-Ocean for both the baroclinic and barotropic subsystems since Omega V2 is non-Boussinesq but hydrostatic.

### 2.3 Requirement: Modularization of baroclinic and barotropic time-stepping methods

Modularity ensures ease of testing and future-proofing of the Omega V2 codebase. Implementing a modular design enables mix-and-match time-stepping schemes of the baroclinic and barotropic subsystems, straightforward integration of alternative time-stepping schemes, and easier maintenance by separating the baroclinic and barotropic time-stepping codes, thereby enhancing flexibility.

## 3. Algorithmic Formulation

### 3.1 Barotropic (external) and baroclinic (internal) mode splitting

The split-explicit method separates ocean velocity into depth-integrated barotropic and depth-dependent baroclinic components. This separation allows computationally expensive baroclinic modes to run at longer timesteps and computationally efficient barotropic modes to run at shorter timesteps, enhancing computational efficiency and accuracy.

The layered discrete governing equations for Omega V2 are described in the {ref}`Omega V1 governing equations <omega-design-governing-eqns-omega1>` design document. The mass, tracer, and velocity equations used as the starting point for the mode split are summarized below.

**Mass:**

$$
\frac{\partial \tilde{h}_{i,k}}{\partial t}
+ \nabla \cdot \left( [\tilde{h}_{k}]_e u_{e,k} \right)
+ \left[ \tilde{W}_{tr} \right]^{\text{top}}_{k}
- \left[ \tilde{W}_{tr} \right]^{\text{top}}_{k+1}
= 0 .
$$ (split-discrete-mass)

**Tracer:**

$$
\frac{\partial \tilde{h}_{i,k}\varphi_{i,k}}{\partial t}
+ \nabla \cdot \left( [\tilde{h}_{i,k}\varphi_{i,k}]_e u_{e,k} \right)
+ \left\{
\left[\varphi \tilde{W}_{tr}\right]^{\text{top}}_{k}
-
\left[\varphi \tilde{W}_{tr}\right]^{\text{top}}_{k+1}
\right\}
= [D_h^\varphi]_{i,k} - [D_v^\varphi]_{i,k},
$$ (split-discrete-tracer)

where

$$
[D_h^\varphi]_{i,k}
= \kappa_{2,e}\nabla^2\varphi_{i,k}
- \kappa_{4,e}\nabla^4\varphi_{i,k},
$$ (split-tracer-horizontal-diffusion)

and

$$
[D_v^\varphi]_{i,k}
= [\tilde{\kappa}_{v}]_{i,k}
\left[\frac{\partial h\varphi}{\partial \tilde{z}}\right]_{i,k}
-
[\tilde{\kappa}_{v}]_{i,k+1}
\left[\frac{\partial h\varphi}{\partial \tilde{z}}\right]_{i,k} .
$$ (split-tracer-vertical-diffusion)

**Velocity:**

$$
\frac{\partial u_{e,k}}{\partial t}
+ \left[ {\bf k} \cdot \nabla \times u_{e,k} + f_v \right]_e u^{\perp}_{e,k}
+ [\nabla K]_e
+ \frac{1}{[\tilde{h}_{i,k}]_e}
\left[\tilde{W}_{tr}\frac{\partial U}{\partial \tilde{z}}\right]_{e,k}
= - [(\alpha \nabla p + \nabla \Phi)]_{e,k}
+ [D_h^{\bf u}]_{e,k}
- [D_v^{\bf u}]_{e,k},
$$ (split-discrete-velocity)

where

$$
[D_h^{\bf u}]_{e,k}
= \nu_{2,e}\nabla^2 {\bf u}_{e,k}
- \nu_{4,e}\nabla^4 {\bf u}_{e,k},
$$ (split-velocity-horizontal-diffusion)

and

$$
[D_v^{\bf u}]_{e,k}
= \frac{1}{[\tilde{h}_{i,k}]_e}
\left\{
[\tilde{\nu}_{v}]^{\text{top}}_{e,k}
\left[\frac{\partial {\bf u}}{\partial \tilde{z}}\right]_{e,k}
-
[\tilde{\nu}_{v}]^{\text{top}}_{e,k+1}
\left[\frac{\partial {\bf u}}{\partial \tilde{z}}\right]_{e,k+1}
\right\} .
$$ (split-velocity-vertical-diffusion)

Define the barotropic velocity, baroclinic velocity, barotropic pressure, and barotropic pressure anomaly as follows.

**Barotropic velocity:**

$$
\overline{u}
\equiv
\frac{1}{\tilde{H}}
\sum_{k=0}^{K_{\max}} \tilde{h}_{k}u_{k}.
$$ (split-barotropic-velocity)

**Baroclinic velocity:**

$$
u'_k \equiv u_k - \overline{u}.
$$ (split-baroclinic-velocity)

**Barotropic pressure:**

$$
B \equiv p^{\text{floor}} - p^{\text{surf}}.
$$ (split-barotropic-pressure)

**Barotropic pressure anomaly:**

$$
B' \equiv B - \rho_0 g b,
$$ (split-barotropic-pressure-anomaly)

so that

$$
B = B' + \rho_0 g b,
$$ (split-barotropic-pressure-reconstruction)

and

$$
p^{\text{floor}}
= p^{\text{surf}} + B' + \rho_0 g b,
$$ (split-floor-pressure-reconstruction)

where $b \equiv -z^{\text{floor}}$ is the geometric bottom depth.

Define the column-integrated pseudo thickness and geometric thickness as follows.

**Column-integrated pseudo thickness**, or mass-equivalent column depth:

$$
\tilde{H}
\equiv
\sum_{k=0}^{K_{\max}} \tilde{h}_{k}
=
\frac{1}{\rho_0 g}\sum_{k=0}^{K_{\max}}
\left(p_{k+1}^{\text{top}} - p_k^{\text{top}}\right)
=
\frac{1}{\rho_0 g}\left(p^{\text{floor}}-p^{\text{surf}}\right)
=
\frac{B}{\rho_0 g}
=
\frac{1}{\rho_0 g}\left(B' + \rho_0 g b\right).
$$ (split-column-pseudo-thickness)

**Column-integrated geometric thickness**, or geometric column depth:

$$
H
=
\sum_{k=0}^{K_{\max}} h_k
=
\sum_{k=0}^{K_{\max}} \frac{\rho_0}{\rho_k}\tilde{h}_k
=
\rho_0 \sum_{k=0}^{K_{\max}} \alpha_k \tilde{h}_k
=
\rho_0 S,
$$ (split-column-geometric-thickness)

where

$$
S = \sum_{k=0}^{K_{\max}} \alpha_k \tilde{h}_k
$$ (split-column-specific-volume)

is the column-integrated specific volume.

Sea-surface height can be diagnosed by

$$
\eta = \rho_0 S - b.
$$ (split-diagnostic-ssh)

The sea-surface height is decomposed into two components:

$$
\eta = \eta^{\text{mass}} + \eta^{\text{steric}},
$$ (split-ssh-decomposition)

where

$$
\eta^{\text{mass}} = \frac{B'}{\rho_0 g},
\qquad
\eta^{\text{steric}} = \rho_0 S - \frac{B}{\rho_0 g}.
$$ (split-mass-and-steric-ssh)

This expression is equivalent to Eq. (10.60) in [Griffies (2012)](https://mom-ocean.github.io/assets/pdfs/MOM5_manual.pdf) and [Madec et al. (2015)](https://epic.awi.de/id/eprint/39698/1/NEMO_book_v6039.pdf), but formulated using $S$ and $B$. The term $\eta^{\text{steric}}$ represents the change in sea-surface height resulting from the non-Boussinesq steric effect, while $\eta^{\text{mass}}$ is the mass-related component of sea-surface height that is advanced in time in the barotropic system as the fast process.

Under the Boussinesq approximation, the steric height $\eta^{\text{steric}}$ reduces to zero:

$$
\eta^{\text{steric}} = H - \tilde{H};
\qquad
\rho \rightarrow \rho_0;
\qquad
H - H = 0.
$$ (split-boussinesq-steric-limit)

A common expression for steric height, following [Griffies (2012)](https://mom-ocean.github.io/assets/pdfs/MOM5_manual.pdf) and [Madec et al. (2015)](https://epic.awi.de/id/eprint/39698/1/NEMO_book_v6039.pdf), is

$$
\eta^{\text{steric}}
= -\int_{-H}^{\eta}
\left(\frac{\rho - \rho_0}{\rho_0}\right) dz .
$$ (split-steric-height-griffies-madec)

The derivation below demonstrates that this definition of $\eta^{\text{steric}}$ is identical to the definition employed in the split system of Omega:

$$
\begin{aligned}
\eta^{\text{steric}}
&= -\int_{-H}^{\eta}
\left(\frac{\rho - \rho_0}{\rho_0}\right) dz \\
&= \int_{-H}^{\eta}
\left(1 - \frac{\rho}{\rho_0}\right) dz \\
&= \int_{-H}^{\eta} dz
- \int_{-H}^{\eta}\frac{\rho}{\rho_0} dz \\
&= \int_{-H}^{\eta} dz
- \int_{-H}^{\eta} d\tilde{z} \\
&= H - \tilde{H} \\
&= \rho_0 S - \frac{B}{\rho_0 g} .
\end{aligned}
$$ (split-steric-height-derivation)


The barotropic continuity and momentum equations are written as follows.

**Barotropic continuity equation:**

$$
\frac{\partial B_i'}{\partial t}
+ \left[
\nabla \cdot \left( [(B_i' + \rho_0 g b_i)]_e \overline{u}_e \right)
\right]_i
= -\rho_0 g Q_i .
$$ (split-barotropic-continuity)

**Barotropic momentum equation:**

$$
\frac{\partial \overline{u}_e}{\partial t}
+ f_e \overline{u}^{\perp}_e
= -[\overline{\alpha}_i]_e [\nabla B_i']_e
+ \overline{G}_e .
$$ (split-barotropic-momentum)

**Baroclinic momentum equation:**

$$
\frac{\partial u'_{e,k}}{\partial t}
= -[f_v]_e u_{e,k}^{\prime\perp}
+ \Gamma_{e,k}
+ [\overline{\alpha}_i]_e [\nabla B_i']_e
- \overline{G}_e,
$$ (split-baroclinic-momentum)

where

$$
\Gamma_{e,k}
\equiv
-[\nabla K]_e
-(\alpha \nabla p + \nabla \Phi)_{e,k}
-
[ {\bf k} \cdot \nabla \times u_{e,k}]_e u^{\perp}_{e,k}
-
\frac{1}{[\tilde{h}_{i,k}]_e}
\left[\tilde{W}_{tr}\frac{\partial U}{\partial \tilde{z}}\right]_{e,k}
+ [D_h^{\bf u}]_{e,k}
- [D_v^{\bf u}]_{e,k}.
$$ (split-gamma-definition)

Here, $\overline{G}_e$ includes all remaining terms in the barotropic equation.

### 3.2 Split-explicit time stepping algorithm

The mode-splitting time-stepping algorithm in Omega-V2 follows the MPAS-Ocean split-explicit framework. The algorithm first advances the baroclinic velocity over the large time step, explicitly subcycles the barotropic mode, and then uses the resulting transport velocity to update pseudo thickness and tracers. This structure keeps the barotropic–baroclinic coupling explicit and modular: the baroclinic update provides the vertically averaged forcing to the barotropic solver, while the barotropic subcycling returns time-averaged transports used for the pseudo thickness and tracer updates. Omega-V2 will provide two split-explicit options: `SE-RK2` and `SE-AB2`.

The `SE-RK2` option, corresponding to the split-explicit second-order Runge–Kutta scheme used in MPAS-Ocean (`split_explicit`), can be described as a split-explicit RK2-like predictor–corrector scheme. In this approach, provisional end-of-step values are first estimated during the first time step iteration. Midpoint states are then constructed by averaging the old and provisional new states, and the tendencies are recomputed using these midpoint estimates during the second time step iteration. Because the baroclinic velocity advance, barotropic subcycling, and pseudo thickness/tracer updates are performed sequentially rather than as a fully synchronized RK update of all prognostic variables, the method is RK2-like rather than a fully stage-synchronous RK2 scheme.

The SE-AB2 option, corresponding to the split-explicit second-order Adams–Bashforth scheme used in MPAS-Ocean (`split_explicit_ab2`), follows the same split-explicit mode-splitting framework but performs only one time step iteration. Instead of constructing midpoint states through a second RK2-like correction, the baroclinic tendencies are advanced using a second-order Adams–Bashforth extrapolation based on the current and previous tendencies. Compared with `SE-RK2`, `SE-AB2` is computationally less expensive because it requires only one time step iteration, but it requires additional storage for previous-step tendencies or tendency-related forcing terms.

For the initial implementation, `SE-RK2` is selected as the baseline split-explicit scheme. `SE-AB2` can be implemented later by extending the same framework with Adams–Bashforth extrapolation of the baroclinic tendencies. Therefore, this section first describes the `SE-RK2` algorithm.


#### 3.2.1 Initialization

Define `NTimeStepIteration` as

$$
\text{NTimeStepIteration} = 2
\quad \text{if} \quad
\text{TimeStepper} = \texttt{SE-RK2}.
$$ (split-ntimestepiter-rk2)

and

$$
\text{NTimeStepIteration} = 1
\quad \text{if} \quad
\text{TimeStepper} = \texttt{SE-AB2},
$$ (split-ntimestepiter-ab2)

If the model is not restarting, `SE-AB2` uses `SE-RK2` for the first time step with `NTimeStepIteration = 2`. From the second time step, `NTimeStepIteration = 1`.

Compute `NBtrSubcycle` as

$$
\text{NBtrSubcycle} = \frac{\text{TimeStep}}{\text{BtrTimeStep}} .
$$ (split-nbtrsubcycle)

Compute the barotropic velocity:

$$
\overline{u}_e
\equiv
\frac{1}{[\tilde{H}_i]_e}
\sum_{k=0}^{K_{\max}}
[\tilde{h}_{i,k}]_e u_{e,k}.
$$ (split-initial-barotropic-velocity)

If the model is restarting, $\overline{u}_e$ is read from the previous time step instead.

Compute the baroclinic velocity:

$$
u'_{e,k} = u_{e,k} - \overline{u}_e.
$$ (split-initial-baroclinic-velocity)

Compute the pressure $p$.

Compute the barotropic pressure:

$$
B_i = p_{i,K_{\max}+1} - p_{i,0}.
$$ (split-initial-barotropic-pressure)

Compute the barotropic pressure anomaly:

$$
B_i' = B_i - \rho_0 g b_i.
$$ (split-initial-barotropic-pressure-anomaly)

Prepare variables before the first iteration:

$$
u^{*}_{e,k} = u^n_{e,k},
\qquad
\tilde{W}^{*}_{i,k} = \tilde{W}^{n}_{i,k},
\qquad
\tilde{h}^{*}_{i,k} = \tilde{h}^{n}_{i,k},
\qquad
\varphi^{*}_{i,k} = \varphi^{n}_{i,k},
\qquad
p^{*}_{i,k} = p^{n}_{i,k},
\qquad
\tilde{z}^{*}_{i,k} = \tilde{z}^{n}_{i,k},
\quad \text{etc.}
$$ (split-initial-provisional-variables)

#### 3.2.2 Stage 1: Baroclinic velocity advance with long time step

This stage advances the baroclinic velocity $u'$ with the long time step and computes the barotropic forcing term $\overline{G}$.

Compute the baroclinic forcing plus the barotropic pressure-gradient contribution:

$$
\Gamma^*_{e,k} + [\overline{\alpha}^{*}_{i}]_e [\nabla B_i^{\prime *}]_e .
$$ (split-stage1-baroclinic-forcing)

Compute the column-integrated pseudo thickness:

$$
\tilde{H}^{*}_i
= \sum_{k=0}^{K} \tilde{h}^{*}_{i,k}.
$$ (split-stage1-column-pseudo-thickness)

Compute the Coriolis term using a centered treatment with two iterations. For
$j = 0, \ldots, \text{NBclIter}-1$, with the default value `NBclIter = 2`, compute $f_e u_{e,k}^{\prime\perp *}$ from $u_{e,k}^{\prime *}$:

$$
u_{e,k}^{\prime\perp *}
=
\sum_{e'\in ECP(e)} \tilde{E}_{e,e'} f_{e'} u_{e',k}^{\prime *}.
$$ (split-stage1-baroclinic-coriolis)

Advance the baroclinic velocity:

$$
u_{e,k}^{\prime n+1}
=
u_{e,k}^{\prime n}
+ \Delta t
\left(
-[f_v]_e u_{e,k}^{\prime\perp *}
+ \Gamma^*_{e,k}
+ [\overline{\alpha}^{*}_{i}]_e [\nabla B_i^{\prime *}]_e
\right).
$$ (split-stage1-baroclinic-advance)

Compute $\overline{G}^{*}_e$:

$$
\overline{G}^{*}_e
=
\frac{1}{[\tilde{H}^{*}_{i}]_e\Delta t}
\sum_{k=0}^{K} \tilde{h}^{*}_{i,k}u_{e,k}^{\prime n+1}.
$$ (split-stage1-gbar)

Set

$$
\overline{G}^{*}_e = 0
\quad \text{for the unsplit algorithm}.
$$ (split-stage1-gbar-unsplit)

Compute the midpoint baroclinic velocity:

$$
u_{e,k}^{\prime n+0.5}
= \frac{1}{2}
\left(
 u_{e,k}^{\prime n}
+ u_{e,k}^{\prime n+1}
- \Delta t\,\overline{G}^{*}_e
\right).
$$ (split-stage1-midpoint-baroclinic-velocity)

#### 3.2.3 Stage 2: Barotropic velocity advance, explicitly subcycled

This stage advances $B'$ and $\overline{u}$ as a coupled system through $2M$ subcycles using the predictor-corrector scheme, ending at time $t+2\Delta t$. For the unsplit algorithm, this stage is skipped; that is, $\overline{u}=0$ and $u=u'$.

The discrete barotropic continuity equation is

$$
B_i^{\prime n+1}
= B_i^{\prime n}
- \Delta t
\left[
\nabla \cdot \left( [(B_i^{\prime n}+\rho_0 g b_i)]_e \overline{u}_e^n \right)
\right]_i
- \Delta t\,\rho_0 g Q_i^n.
$$ (split-stage2-discrete-continuity)

The discrete barotropic momentum update is

$$
\overline{u}^{n+1}_e
= \overline{u}^{n}_e
+ \Delta t
\left(
- f_e \overline{u}^{\perp n}_e
- [\overline{\alpha}_i^{n}\nabla B_i^{\prime n}]_e
+ \overline{G}^{*}_e
\right).
$$ (split-stage2-discrete-momentum)

Initialize the barotropic subcycling variables:

$$
\hat{\overline{u}}^{n}_e = \overline{u}^{n}_e,
\qquad
\hat{B}^{\prime n}_i = B_i^{\prime n}.
$$ (split-stage2-initialization)

For each predictor-corrector subcycle, $m=0,\ldots,2M-1$, use the following steps.

**$\overline{u}$ predictor:**

$$
[\overline{u}^{*}_e]^{n+(m+1)/M}
=
[\overline{u}_e]^{n+m/M}
+
\frac{\Delta t}{M}
\left(
-f_e[\overline{u}_e^{\perp}]^{n+m/M}
-[\overline{\alpha}_i]^n_e[\nabla B_i']^{n+m/M}_e
+ \overline{G}_e
\right).
$$ (split-stage2-u-predictor)

**$B'$ predictor:**

$$
[F^{*}_e]^{m+1}
=
\left([B_i']^{n+m/M}+\rho_0 g b_i\right)_e
\left(
(1-\gamma_1)[\overline{u}_e]^{n+m/M}
+
\gamma_1[\overline{u}^{*}_e]^{n+(m+1)/M}
\right),
$$ (split-stage2-b-predictor-flux)

and

$$
[B_i^{\prime *}]^{n+(m+1)/M}
=
[B_i']^{n+m/M}
+
\frac{\Delta t}{M}
\left(
[-\nabla \cdot [F^{*}_e]^{m+1}]_i
-
\rho_0 g Q_i^n
\right).
$$ (split-stage2-b-predictor)

**$\overline{u}$ corrector:**

$$
[\overline{u}_e]^{n+(m+1)/M}
=
[\overline{u}_e]^{n+m/M}
+
\frac{\Delta t}{M}
\left(
-f_e[\overline{u}_e^{\perp *}]^{n+m/M}
-
[\overline{\alpha}_i]^n_e
\nabla
\left(
(1-\gamma_2)[B_i']^{n+m/M}
+
\gamma_2[B_i^{\prime *}]^{n+(m+1)/M}
\right)_e
+
\overline{G}_e
\right).
$$ (split-stage2-u-corrector)

**$B'$ corrector:**

$$
[F_e]^{m+1}
=
\left[
\left(
(1-\gamma_2)[B_i']^{n+m/M}
+
\gamma_2[B_i^{\prime *}]^{n+(m+1)/M}
\right)
+\rho_0 g b_i
\right]_e
\left(
(1-\gamma_3)[\overline{u}_e]^{n+m/M}
+
\gamma_3[\overline{u}^{*}_e]^{n+(m+1)/M}
\right),
$$ (split-stage2-b-corrector-flux)

and

$$
[B_i^{\prime *}]^{n+(m+1)/M}
=
[B_i']^{n+m/M}
+
\frac{\Delta t}{M}
\left(
[-\nabla \cdot [F_e]^{m+1}]_i
-
\rho_0 g Q_i^n
\right).
$$ (split-stage2-b-corrector)

Compute the time average during subcycling:

$$
\overline{u}_e^{\text{bt}}
=
\frac{1}{2M+1}
\sum_{m=0}^{2M}
[\overline{u}_e]^{n+m/M},
$$ (split-stage2-barotropic-velocity-average)

and

$$
\overline{F}_e^{\text{bt}}
=
\frac{1}{2M}
\sum_{m=0}^{2M-1}
[F_e]^{m+1}.
$$ (split-stage2-barotropic-flux-average)

Then perform the boundary update on $\overline{u}_e^{\text{bt}}$ and $\overline{F}_e^{\text{bt}}$.

#### 3.2.4 Barotropic-baroclinic coupling and barotropic pressure consistency

For the mode-split consistency of the barotropic pressure anomaly between $B'$ from the barotropic mode and $\rho_0 g(\tilde{H}-b)$ from the baroclinic mode, Omega follows the scheme from Hallberg and Adcroft (2009), as implemented in MPAS-Ocean.

The barotropic update of $B'$ is given by

$$
\frac{B_i^{\prime n+1}-B_i^{\prime n}}{\Delta t}
+
\left[
\nabla \cdot
\overline{\left([(B_i' + \rho_0 g b_i)]_e \overline{u}_e\right)}^{\text{bt}}
\right]_i
= -\rho_0 g Q_i .
$$ (split-btr-update-consistency)

Here, $\overline{\varphi}^{\text{bt}}$ denotes a time-averaged quantity from the barotropic subcycles, and $n$ indicates the baroclinic time step.

The velocity correction $u^{\text{co}}$ is written as

$$
u_{e,k}^{\text{co}}
=
\left\{
\overline{\left([(B_i' + \rho_0 g b_i)]_e \overline{u}_e\right)}^{\text{bt}}
-
\sum_{k=0}^{K}
[\tilde{h}_i^{*}]_{e,k}
\left(
\overline{u}_e^{\text{bt}}
+u_{e,k}^{\prime n+0.5}
+u_{e,k}^{\text{bolus}*}
\right)
\right\}
\bigg/
[\tilde{H}_i^{*}]_e .
$$ (split-velocity-correction)

The asterisk indicates the provisional variable that is updated during the baroclinic time step iteartion; the most recent available value is always used for forcing terms.

The transport velocity $u^{\text{tr}}$ is defined as

$$
u^{\text{tr}}_{e,k}
=
\overline{u}_e^{\text{bt}}
+u_{e,k}^{\prime n+0.5}
+u_{e,k}^{\text{bolus}*}
+u_{e,k}^{\text{co}}.
$$ (split-transport-velocity)

The transport velocity is used to compute vertical transport velocity and horizontal transport for both pseudo thickness and tracers.

For the unsplit algorithm, the above processes are skipped except that

$$
u_{e,k}^{\text{tr}}
= u_{e,k}^{\prime n+0.5}
+ u_{e,k}^{\text{bolus}*},
$$ (split-unsplit-transport-velocity)

where

$$
u_{e,k}^{\prime n+0.5} = u_{e,k}^{n+0.5}.
$$ (split-unsplit-midpoint-relation)

#### 3.2.5 Stage 3: Update tracers and diagnostics

Compute $\tilde{W}_{i,k}^{*}$ using $u_{e,k}^{\text{tr}}$.

Compute pseudo thickness tendencies using $u_{e,k}^{\text{tr}}$:

$$
\tilde{h}_{i,k}^{n+1}
=
\tilde{h}_{i,k}^{n}
-
\Delta t\,
\nabla \cdot
\left([\tilde{h}_{k}^{*}]_e u_{e,k}^{\text{tr}}\right)
-
\Delta t
\left(
[\tilde{W}_{tr}]^{\text{top}}_{k}
-
[\tilde{W}_{tr}]^{\text{top}}_{k+1}
\right).
$$ (split-stage3-pseudo-thickness-update)

Compute tracer tendencies using $u_{e,k}^{\text{tr}}$:

$$
\varphi_{i,k}^{n+1}
=
\varphi_{i,k}^{n}
-
\Delta t\,
\nabla \cdot
\left([
\tilde{h}_{i,k}^{*}]_e
[\varphi_{i,k}^{*}]_e
u_{e,k}^{\text{tr}}
\right)
-
\Delta t
\left(
[\varphi^{*}\tilde{W}_{tr}]^{\text{top}}_{i,k}
-
[\varphi^{*}\tilde{W}_{tr}]^{\text{top}}_{i,k+1}
\right).
$$ (split-stage3-tracer-update)

#### 3.2.6 Reset variables

If iterating, reset the provisional variables as follows:

$$
u^{\prime *} = u^{\prime n+0.5} \quad \text{from Stage 1},
$$ (split-reset-baroclinic-velocity)

$$
\overline{u}^{*} = \overline{u}^{\text{bt}} \quad \text{from Stage 2},
$$ (split-reset-barotropic-velocity)

$$
u^{*} = \overline{u}^{*} + u^{\prime *},
$$ (split-reset-full-velocity)

$$
\psi^{*}
=
\frac{1}{2}\left(\psi^{n}+\psi^{n+1}\right)
\quad \text{for pseudo thickness and tracers},
$$ (split-reset-provisional-psi)

$$
\tilde{H}^{*} = \sum_{k=0}^{K} \tilde{h}_k^{*},
$$ (split-reset-column-pseudo-thickness)

and

$$
B^{\prime *} = \rho_0 g\left(\tilde{H}^{*}-b\right).
$$ (split-reset-barotropic-pressure-anomaly)

Diagnostic variables are then updated.

After the final iteration,

$$
u^{\prime n+1} \quad \text{is obtained from Stage 1},
$$ (split-final-baroclinic-velocity)

$$
\overline{u}^{n+1} = \overline{u}^{\text{bt}} \quad \text{from Stage 2},
$$ (split-final-barotropic-velocity)

$$
u^{n+1} = \overline{u}^{n+1} + u^{\prime n+1},
$$ (split-final-full-velocity)

$$
\psi^{n+1} \quad \text{is retained for pseudo thickness and tracers},
$$ (split-final-psi)

$$
\tilde{H}^{n+1} = \sum_{k=0}^{K} \tilde{h}_k^{n+1},
$$ (split-final-column-pseudo-thickness)

and

$$
B^{\prime n+1} = \rho_0 g\left(\tilde{H}^{n+1}-b\right).
$$ (split-final-barotropic-pressure-anomaly)

Diagnostic variables are then updated.

### 3.3 Unsplit time stepping algorithm

The unsplit algorithm follows the same overall structure as the split-explicit algorithm, except that the full velocity is advanced directly in Stage 1 and the barotropic subcycling stage is skipped.

#### 3.3.1 Initialization

Compute the pressure $p$.

Prepare variables before the first iteration:

$$
u^{*}_{e,k} = u^n_{e,k},
\qquad
\tilde{W}^{*}_{i,k} = \tilde{W}^{n}_{i,k},
\qquad
\tilde{h}^{*}_{i,k} = \tilde{h}^{n}_{i,k},
\qquad
\varphi^{*}_{i,k} = \varphi^{n}_{i,k},
\qquad
p^{*}_{i,k} = p^{n}_{i,k},
\qquad
\tilde{z}^{*}_{i,k} = \tilde{z}^{n}_{i,k},
\quad \text{etc.}
$$ (unsplit-initial-provisional-variables)

#### 3.3.2 Stage 1: Velocity advance

This stage advances the full velocity $u$.

Compute $\Gamma_{e,k}^{*}$:

$$
\Gamma_{e,k}^{*}
=
-[{\bf k}\cdot\nabla\times u_{e,k}]_e u_{e,k}^{\perp}
-[\nabla K]_e
-
\frac{1}{[\tilde{h}_{i,k}]_e}
\left[\tilde{W}_{tr}\frac{\partial U}{\partial \tilde{z}}\right]_{e,k}
-(\alpha\nabla p+\nabla\Phi)_{e,k}
+[D_h^{\bf u}]_{e,k}
-[D_v^{\bf u}]_{e,k}.
$$ (unsplit-gamma-definition)

Compute the column-integrated pseudo thickness:

$$
\tilde{H}_i^{*}
= \sum_{k=0}^{K} \tilde{h}^{*}_{i,k}.
$$ (unsplit-column-pseudo-thickness)

Compute the Coriolis term using a centered treatment with two iterations. For $j=0,\ldots,\text{NBclIter}-1$, with the default value `NBclIter = 2`, compute $f_e u_{e,k}^{\perp *}$ from $u^{*}_{e,k}$:

$$
u_{e,k}^{\perp *}
=
\sum_{e'\in ECP(e)} \tilde{E}_{e,e'} f_{e'} u_{e',k}^{*}.
$$ (unsplit-coriolis)

Advance the velocity:

$$
u_{e,k}^{n+1}
=
u_{e,k}^{n}
+
\Delta t
\left(
-[f_v]_e u_{e,k}^{\perp *}
+
\Gamma_{e,k}^{*}
\right).
$$ (unsplit-velocity-advance)

Set

$$
\overline{G}_e^{*} = 0.
$$ (unsplit-gbar)

Compute the midpoint velocity:

$$
u_{e,k}^{n+0.5}
=
\frac{1}{2}
\left(u_{e,k}^{n}+u_{e,k}^{n+1}\right).
$$ (unsplit-midpoint-velocity)

#### 3.3.3 Stage 2: Barotropic velocity advance, explicitly subcycled

For the unsplit time stepper, $\overline{u}=0$. This stage is skipped.

#### 3.3.4 Stage 3: Update tracers and diagnostics

Compute $u_{e,k}^{\text{tr}}$:

$$
u_{e,k}^{\text{tr}}
= u_{e,k}^{n+0.5}
+ u_{e,k}^{\text{bolus}*}.
$$ (unsplit-transport-velocity)

Compute $\tilde{W}_{i,k}^{*}$ using $u_{e,k}^{\text{tr}}$.

Compute pseudo thickness tendencies using $u_{e,k}^{\text{tr}}$:

$$
\tilde{h}_{i,k}^{n+1}
=
\tilde{h}_{i,k}^{n}
-
\Delta t\,
\nabla \cdot
\left([\tilde{h}_{k}^{*}]_e u_{e,k}^{\text{tr}}\right)
-
\Delta t
\left(
[\tilde{W}_{tr}]^{\text{top}}_{k}
-
[\tilde{W}_{tr}]^{\text{top}}_{k+1}
\right).
$$ (unsplit-pseudo-thickness-update)

Compute tracer tendencies using $u_{e,k}^{\text{tr}}$:

$$
\varphi_{i,k}^{n+1}
=
\varphi_{i,k}^{n}
-
\Delta t\,
\nabla \cdot
\left([
\tilde{h}_{i,k}^{*}]_e
[\varphi_{i,k}^{*}]_e
u_{e,k}^{\text{tr}}
\right)
-
\Delta t
\left(
[\varphi^{*}\tilde{W}_{tr}]^{\text{top}}_{i,k}
-
[\varphi^{*}\tilde{W}_{tr}]^{\text{top}}_{i,k+1}
\right).
$$ (unsplit-tracer-update)

#### 3.3.5 Reset variables

If iterating, reset the provisional variables as follows:

$$
u^{*} = u^{n+0.5} \quad \text{from Stage 1},
$$ (unsplit-reset-velocity)

$$
\psi^{*}
=
\frac{1}{2}\left(\psi^{n}+\psi^{n+1}\right)
\quad \text{for pseudo thickness and tracers},
$$ (unsplit-reset-provisional-psi)

$$
\tilde{H}^{*} = \sum_{k=0}^{K} \tilde{h}_k^{*}.
$$ (unsplit-reset-column-pseudo-thickness)

Diagnostic variables are then updated.

After the final iteration,

$$
u^{n+1} \quad \text{is obtained from Stage 1},
$$ (unsplit-final-velocity)

$$
\psi^{n+1} \quad \text{is retained for pseudo thickness and tracers},
$$ (unsplit-final-psi)

$$
\tilde{H}^{n+1} = \sum_{k=0}^{K} \tilde{h}_k^{n+1}.
$$ (unsplit-final-column-pseudo-thickness)

Diagnostic variables are then updated.

## 4. Design
### 4.1 Data types and parameters
#### 4.1.1 Parameters
- `NTimeStepIteration`: Number of baroclinic iterations per timestep (default: 2)

#### 4.1.2 Class/structs/data types

### 4.2 Methods

## 5. Verification and testing
### 5.1 Unit testing

### 5.2 Polaris tests
- ** Inertial gravity wave test
- ** Internal tide test
- ** Baroclinic channel test

<!--
## References

- Higdon, R. L. (2005). Reference used for the split-explicit barotropic-baroclinic time-stepping formulation.
- Hallberg, R., & Adcroft, A. (2009). Barotropic-baroclinic coupling and pressure consistency scheme. https://doi.org/10.1016/j.ocemod.2009.02.008
- Griffies (2012) and Madec et al. (2015) are referenced for the steric sea-surface height formulation used in [](#split-steric-height-griffies-madec).
--!>
