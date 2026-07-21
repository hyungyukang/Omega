(omega-user-pgrad)=

# Pressure Gradient

The pressure gradient term in the momentum equation represents the force per unit
mass due to horizontal variations in pressure and geopotential. This term is
essential for capturing the dynamics of
ocean circulation, including both barotropic and baroclinic motions.

## Physical Background

In the layered non-Boussinesq momentum equation solved in Omega, the pressure
gradient tendency for each edge and layer includes three contributions:

1. **Pressure gradient**: The centered horizontal gradient of layer-midpoint
   pressure multiplied by the edge-averaged specific volume.

2. **Geometric-height gradient**: The centered horizontal gradient of
   layer-midpoint geometric height multiplied by gravity. This accounts for
   tilted layers in a general vertical coordinate.

3. **External geopotential forcing**: Contributions from the tidal potential and
   the self-attraction and loading (SAL) terms. These represent gravitational
   forcing from astronomical tides and the deformation of the solid Earth and ocean
   surface in response to the ocean mass distribution. These terms are currently set
   to zero and will be provided by a future tidal forcing module.

## Configuration Options

The pressure gradient method is configured in the input YAML file under the
`PressureGrad` section:

```yaml
PressureGrad:
   PressureGradType: 'centered'
```

### Available Methods

**Centered Difference** (`'centered'` or `'Centered'`)
- Computes the pressure gradient using a centered finite-difference approximation
  of the pressure and layer-midpoint geometric-height gradients
- Suitable for global ocean simulations without ice shelf cavities
- Default option

**High-Order** (`'HighOrder1'`)
- Uses monotonic piecewise-linear temperature and salinity reconstruction with
  three-point Gauss-Legendre integration through each layer and across each
  edge
- Evaluates TEOS-10 specific volume at the quadrature pressure and reconstructs
  geometric height with the same hydrostatic quadrature
- Intended for simulations with ice shelf cavities and steep bathymetry where the
  centered scheme may be inaccurate

## Dependencies

The pressure gradient calculation requires the following Omega components to be
initialized first:

- [**Horizontal Mesh**](omega-user-horz-mesh): provides mesh geometry including
  distances between cell centers and edge connectivity
- [**Vertical Coordinate**](omega-user-vert-coord): provides pressure at layer
  midpoints and interfaces and geometric midpoint heights ($z$)
- [**Equation of State**](omega-user-eos): provides the specific volume field
- [**Tracers**](omega-user-tracers): provides Conservative Temperature and
  Absolute Salinity
