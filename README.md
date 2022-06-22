# vortex_tube_CFD_tools
Scripts to mesh Ranque-Hilsch vortex tube geometries, run vortex tube simulations, and post-process vortex tube simulation results

## Meshing
Steps to produce a Ranque-Hilsch vortex tube mesh for CFX or Fluent in an environment with [OpenFOAM](https://openfoam.org/) v8 and [blockmeshbuilder](https://github.com/NauticalMile64/blockmeshbuilder) (>= v0.2.0) installed:
1. `$ cd vortex_tube_OFCase`
2. `$ python3 vortex_tube_OFCase/vortex_tube.py` generates the `blockMeshDict` file
3. `$ blockMesh` generates the hex mesh in the `constant/polyMesh` folder.
4. `$ foamMeshToFluent` converts the polyMesh into a `.msh` file readable by CFX or Fluent 

The other files in this minimal OpenFOAM case are required to keep OpenFOAM's `blockMesh` command happy, though any user is welcome to update the case structure to be able to run using OpenFOAM.

## Simulation Setup
An ANSYS CFX simulation can be initialized using [ModelA1_setup_SASSST.ccl](ModelA1_setup_SASSST.ccl) as a starting point.
