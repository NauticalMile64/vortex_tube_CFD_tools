# vortex_tube_CFD_tools
Scripts to mesh Ranque-Hilsch vortex tube geometries, run vortex tube simulations, and post-process vortex tube simulation results.

## Meshing
Steps to produce a Ranque-Hilsch vortex tube mesh for CFX or Fluent in an environment with [OpenFOAM](https://openfoam.org/) v8 and [blockmeshbuilder](https://github.com/NauticalMile64/blockmeshbuilder) (>= v0.2.0) installed:
1. `$ cd vortex_tube_OFCase`
2. `$ python3 vortex_tube_OFCase/vortex_tube.py` generates the `blockMeshDict` file
3. `$ blockMesh` generates the hex mesh in the `constant/polyMesh` folder.
4. `$ foamMeshToFluent` converts the polyMesh into a `.msh` file readable by CFX or Fluent 

The other files in this minimal OpenFOAM case are required to keep OpenFOAM's `blockMesh` command happy, though any user is welcome to update the case structure to be able to run using OpenFOAM.

## Simulation Setup
An ANSYS CFX simulation can be initialized using [ModelA1_setup_SASSST.ccl](ModelA1_setup_SASSST.ccl) as a starting point.

## Results Post-processing

### Streamline extraction
If CFX or Fluent results files are post-processed using CFD-Post, the [streamline_post_process.cst](streamline_post_process.cst) file will enable export of streamlines and relevant data from the simulation.

This script assumes that a locus of stagnation points are present somewhere along the surface of the hot exit plug. Straight lines are drawn from the plug tip to the rear of the plug at a constant circumferential angle. Then the location where the in-plane velocity is minimized is used as a starting point for the streamline trace (this is assumed to be a stagnation point within the flow field). Streamlines are traced backwards in the flow field towards the inlets.

Users need to make the following adjustments:

1. Change the number of streamlines by adjusting the `$numDivs` parameter on line 1.
2. Update the case name on line 17 according to the naming of the results file.
3. Adjust the variables `$zp1` and `$rp1` if necessary. In the associated publication, it was important to use the other two (commented) values for these variables for the k-omega SST turbulence model results, as the default method extracted streamlines that would re-circulate towards the hot plug.
4. Once satisfied with the visual results from steps 1-3, replace the streamline export path on line 400 with a valid path and file name prefix, uncomment the export command on line 402.
