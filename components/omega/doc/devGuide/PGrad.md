(omega-dev-pgrad)=

# Pressure Gradient (PGrad)

Omega includes a `PressureGrad` class that computes horizontal pressure gradient
tendencies for the non-Boussinesq momentum equation. The implementation supports
centered and TEOS-10 pressure-integrated schemes, with a placeholder for future
high-order methods. The class follows the same factory pattern used by other Omega
modules.

## PressureGradType enum

An enumeration of the available pressure gradient schemes is defined in `PGrad.h`:

```c++
enum class PressureGradType {
   Centered,
   Integrated,
   HighOrder1,
   HighOrder2
};
```

This is used to select which pressure gradient method is applied at runtime.

## Initialization

An instance of `PressureGrad` requires both a [`HorzMesh`](#omega-dev-horz-mesh) and
a [`VertCoord`](#omega-dev-vert-coord), so these classes and all of their dependencies
must be initialized before `PressureGrad` can be initialized. The static method:

```c++
OMEGA::PressureGrad::init();
```

initializes the default `PressureGrad` instance using the default `HorzMesh` and
`VertCoord` instances and the global Omega configuration. A pointer to the default
instance can be retrieved at any time using:

```c++
OMEGA::PressureGrad* DefPGrad = OMEGA::PressureGrad::getDefault();
```

## Creating additional instances

Additional named instances can be created using:

```c++
OMEGA::PressureGrad* MyPGrad =
    OMEGA::PressureGrad::create("MyPGrad", Mesh, VCoord, Options);
```

where `Mesh` is a pointer to a `HorzMesh`, `VCoord` is a pointer to a `VertCoord`,
and `Options` is a pointer to the `Config`. A named instance can be retrieved later
using:

```c++
OMEGA::PressureGrad* MyPGrad = OMEGA::PressureGrad::get("MyPGrad");
```

## Constructor behavior

The constructor reads the `PressureGrad` section from the configuration, stores
references to mesh and vertical coordinate data needed for computation, and enables
the appropriate functor based on the configured `PressureGradType`. It also allocates
placeholder arrays for tidal potential and self-attraction/loading, which are intended
to be populated by a future tidal forcing module. Currently these arrays are initialized
to zero.

## Computing the pressure gradient

To compute pressure gradient tendencies and accumulate them into a tendency array:

```c++
PGrad->computePressureGrad(Tend, PressureMid, PressureInterface, SpecVol,
                           GeomZInterface, PseudoThick, ConservTemp,
                           AbsSalinity);
```

where:

- `Tend` is a 2D array `(NEdgesAll × NVertLayers)` that the pressure gradient
  tendency is accumulated into
- `PressureMid` and `PressureInterface` contain gauge pressure at cell midpoints
  and interfaces
- `SpecVol` contains layer-centered specific volume
- `GeomZInterface` contains geometric height at layer interfaces
- `PseudoThick` contains layer pseudo-thickness
- `ConservTemp` and `AbsSalinity` contain conservative temperature and absolute
  salinity; the integrated functor uses both arrays directly

The method uses hierarchical Kokkos parallelism: an outer `parallelForOuter` loop
iterates over edges, and an inner `parallelForInner` loop iterates over vertical
chunks. The appropriate functor is dispatched based on `PressureGradChoice`.

## Functors

### PressureGradCentered

This functor implements a centered difference approximation of the pressure gradient
tendency. For each edge, it first computes the layer-invariant tidal and
self-attraction/loading contribution:

```
GradGeoPot = grad(TidalPotential) + grad(SelfAttractionLoading)
```

Then, for each vertical layer `K`, it computes three terms:

1. **Montgomery potential gradient**: The average of the horizontal gradients of the
   Montgomery potential ($\alpha p + g z$) at the top (interface `K`) and bottom
   (interface `K+1`) of the layer. This compactly represents the combined effect
   of the pressure gradient and the geopotential contribution from tilted coordinate
   surfaces.

2. **Specific volume correction**: A correction term equal to the edge-averaged
   pressure at mid-layer multiplied by the horizontal gradient of specific volume.
   This accounts for horizontal density variations that are not captured by the
   Montgomery potential form.

3. **Tidal and geopotential forcing** (`GradGeoPot`): The external geopotential
   contribution from tidal forcing and self-attraction/loading, applied uniformly
   across all layers at an edge.

The tendency update for each layer is:

```
Tend(IEdge, K) += EdgeMask(IEdge, K) * (-GradMontPot + PGradAlpha - GradGeoPot)
```

where `EdgeMask` is applied to enforce land boundary conditions. The functor operator
signature is:

```c++
KOKKOS_FUNCTION void operator()(const Array2DReal &Tend, I4 IEdge, I4 KChunk,
                                const Array2DReal &PressureMid,
                                const Array2DReal &PressureInterface,
                                const Array2DReal &GeomZInterface,
                                const Array1DReal &TidalPotential,
                                const Array1DReal &SelfAttractionLoading,
                                const Array2DReal &SpecVol) const;
```

### PressureGradIntegrated

This functor reconstructs one thermodynamic state per edge and layer by
arithmetic averaging:

```text
CtEdge = 0.5 * (ConservTemp(ICell0, K) + ConservTemp(ICell1, K))
SaEdge = 0.5 * (AbsSalinity(ICell0, K) + AbsSalinity(ICell1, K))
```

At each layer interface, it evaluates

```text
Residual = integral(alpha(SaEdge, CtEdge, p), pCell0, pCell1)
           + Gravity * (zCell1 - zCell0)
```

where `alpha` is the TEOS-10 75-term specific-volume polynomial. The pressure
integral is evaluated analytically in normalized pressure, rather than by
numerical quadrature. The top and bottom interface residuals are averaged and
converted to an edge-normal acceleration:

```text
PGrad = -0.5 * (ResidualTop + ResidualBot) / DcEdge
Tend(IEdge, K) += EdgeMask(IEdge, K) * (PGrad - GradGeoPot)
```

Its functor operator signature is:

```c++
KOKKOS_FUNCTION void operator()(
    const Array2DReal &Tend, I4 IEdge, I4 KChunk,
    const Array2DReal &PressureInterface,
    const Array2DReal &GeomZInterface,
    const Array1DReal &TidalPotential,
    const Array1DReal &SelfAttractionLoading,
    const Array2DReal &ConservTemp,
    const Array2DReal &AbsSalinity) const;
```

### PressureGradHighOrder

This functor is a placeholder for a future high-order pressure gradient implementation
suitable for ice shelf cavities and complex bathymetry. Currently it performs no
computation (a no-op).

## Configuration

The pressure gradient type is selected in the input YAML file:

```yaml
PressureGrad:
   PressureGradType: Centered
```

Valid options for `PressureGradType` are:

- `Centered` or `centered`: centered difference approximation (default)
- `Integrated` or `integrated`: analytic TEOS-10 pressure-integrated method
- `HighOrder1`: first high-order method (placeholder, future implementation)

If an unrecognized value is provided, the implementation falls back to the centered
scheme and logs an informational message.

## Data members

The `PressureGrad` class stores the following key data:

| Member | Type | Description |
| ------ | ---- | ----------- |
| `NEdgesAll` | `I4` | Total number of edges including halo |
| `NVertLayers` | `I4` | Number of vertical layers |
| `NVertLayersP1` | `I4` | Number of vertical layers plus one |
| `MinLayerEdgeBot` | `Array1DI4` | Minimum active layer index for each edge |
| `MaxLayerEdgeTop` | `Array1DI4` | Maximum active layer index for each edge |
| `TidalPotential` | `Array1DReal` | Tidal potential (placeholder, currently zero) |
| `SelfAttractionLoading` | `Array1DReal` | Self-attraction and loading term (placeholder, currently zero) |
| `CenteredPGrad` | `PressureGradCentered` | Centered pressure gradient functor |
| `IntegratedPGrad` | `PressureGradIntegrated` | TEOS-10 pressure-integrated functor |
| `HighOrderPGrad` | `PressureGradHighOrder` | High-order pressure gradient functor |
| `PressureGradChoice` | `PressureGradType` | Selected pressure gradient method |

## Removal

To remove all `PressureGrad` instances:

```c++
OMEGA::PressureGrad::clear();
```

To remove a specific named instance:

```c++
OMEGA::PressureGrad::erase("MyPGrad");
```
