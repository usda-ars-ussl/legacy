# U.S. Salinity Laboratory Legacy Codes

This repository contains computer programs no longer under active development.
The programs were developed over the years by the the U.S. Salinity Laboratory,
USDA-ARS, Riverside, CA, in collaboration with the Department of Environmental
Sciences of the University of California, Riverside. The programs all deal with
various aspects of water and/or solute movement in the subsurface.

PDF versions of the program manuals are in the “docs” subdirectory.

# License:

U.S.: [Public Domain](https://www.usa.gov/publicdomain/label/1.0)
International: [CC0](https://creativecommons.org/publicdomain/zero/1.0)

# Contents

The following programs are in the repo:

**CFITM**: Fortran 77 program for determining equilibrium transport parameter from
solute displacement experiments. [van Genuchten, 1980]

**CFITIM**: Fortran 77 program for determining nonequilibrium transport parameter
from miscible displacement experiments. [van Genuchten, 1981]

**CHAIN**: Fortran 77 program for solving analytical solutions for one-dimensional
solute transport involving sequential first-order decay chain reactions.
Several examples illustrate the use of this program for problems involving the
migration of radionuclide decay chains, or the simultaneous movement of various
interacting nitrogen or pesticides species. [van Genuchten, 1985]

**CHAIN_2D**: Fortran 77 program for simulating two-dimensional variably-saturated
water flow, heat movement, and the transport of solutes involved in sequential
first-order decay reactions. The flow equation incorporates a sink term to
account for water uptake by plant roots. The water flow part of the model can
deal with prescribed head, gradient, and flux boundaries, as well as boundaries
controlled by atmospheric conditions. The program considers free drainage
boundary conditions as well as a simplified representation of nodal drains
using results of electric analog experiments is also included. The heat
transport equation considers transport due to conduction and convection with
flowing water. The solute transport equations consider convectivedispersive
transport in the liquid phase, as well as diffusion in the gaseous phase. The
transport equations also include provisions for nonlinear nonequilibrium
reactions between the solid and liquid phases, linear equilibrium reactions
between the liquid and gaseous phases, zero-order production, and two
first-order degradation reactions: one which is independent of other solutes,
and one which provides the coupling between solutes involved in the sequential
first-order decay reactions. The program may be used to analyze water and
solute movement in unsaturated, partially saturated, or fully saturated porous
media.  [Šimůnek et al., 1994]

**CXTFIT**: Fortran 77 program for estimating solute transport parameters using
nonlinear least-squares parameter estimation.  The program may also be used to
solve the direct or forward problem to determine the concentration as a
function of time and/or position. Three different one-dimensional transport
models are included: (i) the conventional equilibrium convection dispersion
equation (CDE); (ii) the chemical and physical nonequilibrium CDE's; (iii) a
stochastic stream tube model based upon the local-scale CDE, and assuming
either equilibrium or nonequilibrium transport.  [Toride et al., 1995]

**SALT**: Fortran 77 program for analyzing experimentally derived crop salt
tolerance data. The program uses a non-linear least squares inversion method to
find the unknown parameters in several salt tolerance response functions. One
of the three models inlcuded in the program is the familiar piecewise linear
response function. Application of this function leads to estimates for the
salinity threshold and the slope of the response curve. Two alternative types
of salinity response functions are also considered.  [van Genuchten, 1983]

**SWMS_2D**: Fortran 77 program for simulating two-
dimensional water flow and solute transport in variably saturated media.
The program numerically solves the Richards equation for saturatedunsaturated
water flow and the convection-dispersion equation for solute
transport. The flow equation incorporates a sink term to account for water
uptake by plant roots. The transport equation includes provisions for linear
equilibrium adsorption, zero-order production, and first-order degradation.
The program may be used to analyze water and solute movement in
unsaturated, partially saturated, or fully saturated porous media.
[Šimůnek et al., 1994]

**SWMS_3D**: Three-dimensional version of SWMS_2D. [Simunek et al., 1995]

**3DADE**: Fortran 77 program evaluating analytical solutions for two- and
three-dimensional equilibrium solute transport.  [Leij and Bradford, 1994]

**N3DADE**: Fortran 77 program evaluating analytical solutions for two- and
three-dimensional non-equilibrium solute transport.  [Leij and Toride, 1997]

# References

Leij, F. J., and S. A. Bradford, 3DADE: A computer program for evaluating
three-dimensional equilibrium solute transport in porous media, Research Report
No. 134, U. S. Salinity Laboratory, USDA, ARS, Riverside, CA, 1994.

Leij, F. J., and N. Toride, N3DADE: A computer program for evaluating
nonequilibrium three-dimensional equilibrium solute transport in porous media,
Research Report No. 143, U. S. Salinity Laboratory, USDA, ARS, Riverside, CA, 1997.

Šimůnek, J., T. Vogel, and M. Th. van Genuchten, The SWMS_2D code for
simulating water flow and solute transport in two-dimensional variably
saturated media. Version 1.2, Research Report No. 132, U.S. Salinity
Laboratory, USDA, ARS, Riverside, California, 196 pp., 1994.

Šimůnek, J., and M. Th. van Genuchten, The CHAIN_2D code for simulating
two-dimensional movement of water flow, heat, and multiple solutes in
variably-saturated porous media, Version 1.1, Research Report No 136, U.S.
Salinity Laboratory, USDA, ARS, Riverside, California, 205pp., 1994.

Šimůnek, J., K. Huang, and M. Th. van Genuchten, The SWMS_3D code for
simulating water flow and solute transport in three-dimensional variably
saturated media. Version 1.0, Research Report No. 139, U.S. Salinity
Laboratory, USDA, ARS, Riverside, California, 155 pp., 1995.

Toride, N., F. J. Leij, and M. Th. van Genuchten, The CXTFIT code for
estimating transport parameters from laboratory or field tracer experiments.
Version 2.0, Research Report No. 137, U. S. Salinity Laboratory, USDA, ARS,
Riverside, CA, 1995.

van Genuchten, M. Th., Determining transport parameters from solute
displacement experiments, Research Report No. 118, U. S. Salinity Laboratory,
USDA, ARS, Riverside, CA, 1980.

van Genuchten, M. Th., Non-equilibrium transport parameters from miscible
displacement experiments, Research Report No. 119, U. S. Salinity Laboratory,
USDA, ARS, Riverside, CA, 1981.

van Genuchten, M. Th., Analyzing crop salt tolerance data: Model description
and user's manual, Research Report No. 120, U. S. Salinity Laboratory,
USDA, ARS, Riverside, CA, 1983.

van Genuchten, M. Th., Convective-dispersive transport of solutes involved in
sequential first-order decay reactions, Computers & Geosciences, 11(2),
129-147, 1985.
