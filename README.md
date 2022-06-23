[![DOI](https://zenodo.org/badge/506301276.svg)](https://zenodo.org/badge/latestdoi/506301276)

# vortex_tube_CFD_tools
Scripts to mesh Ranque-Hilsch vortex tube geometries, run vortex tube simulations, and post-process vortex tube simulation results.

## Related Publications
The scripts herein were used to pre- and post-process the data in the following journal article:

> N. J. Dyck, M. J. Parker, and A. G. Straatman, “The Impact of Boundary Treatment and Turbulence Model on CFD Simulations of the Ranque-Hilsch Vortex Tube,” Int. J. Refrig. In Press.

More information on the scripts found in this repository can be found in the Data in Brief article which will be linked here after its publication.

## Acquiring Simulation Data
Users wishing to understand how the scripts work before running them on their own simulation results files can download the simulation results (`.cgns`) and streamline data files (`.csv`) from the Zenodo repository at [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6625406.svg)](https://doi.org/10.5281/zenodo.6625406). The data files can be copied into the root of a clone of this repository and the Python scripts should run as expected, or the filepaths in the Python scripts can be adjusted.

## Meshing
Steps to produce a Ranque-Hilsch vortex tube mesh for CFX or Fluent in an environment with [OpenFOAM](https://openfoam.org/) v8 and [blockmeshbuilder](https://github.com/NauticalMile64/blockmeshbuilder) (>= v0.2.0) installed:
1. `$ cd vortex_tube_OFCase`
2. `$ python3 vortex_tube_OFCase/vortex_tube.py` generates the `blockMeshDict` file
3. `$ blockMesh` generates the hex mesh in the `constant/polyMesh` folder.
4. `$ foamMeshToFluent` converts the polyMesh into a `.msh` file readable by CFX or Fluent 

The other files in this minimal OpenFOAM case are required to keep OpenFOAM's `blockMesh` command happy, though any user is welcome to update the case structure to be able to run using OpenFOAM.

The mesh can be inspected using [Paraview](https://www.paraview.org/) by opening the vortex_tube.foam and selecting the `Surface With Edges` representation.

## Simulation Setup
An ANSYS CFX simulation can be initialized using [ModelA1_setup_SASSST.ccl](ModelA1_setup_SASSST.ccl) as a starting point.

## Results Post-processing
This section describes how streamlines can be extracted from a completed vortex tube simulation using CFD-Post, and discusses scripts to analyze the energy transfer across the extracted streamlines.

### Streamline extraction
If CFX or Fluent results files are post-processed using CFD-Post, the [streamline_post_process.cst](streamline_post_process.cst) file will enable export of streamlines and relevant data from the simulation.

This script assumes that a locus of stagnation points are present somewhere along the surface of the hot exit plug. Straight lines are drawn from the plug tip to the rear of the plug at a constant circumferential angle. Then the location where the in-plane velocity is minimized is used as a starting point for the streamline trace (this is assumed to be a stagnation point within the flow field). Streamlines are traced backwards in the flow field towards the inlets.

Users need to make the following adjustments:

1. Change the number of streamlines by adjusting the `$numDivs` parameter on line 1.
2. Update the case name on line 17 according to the naming of the results file.
3. Adjust the variables `$zp1` and `$rp1` if necessary. In the associated publication, it was important to use the other two (commented) values for these variables for the k-omega SST turbulence model results, as the default method extracted streamlines that would re-circulate towards the hot plug.
4. Once satisfied with the visual results from steps 1-3, replace the streamline export path on line 400 with a valid path and file name prefix, uncomment the export command on line 402.

### Streamline Analysis
The streamlines extracted using the methods in the above section (or any streamline really) can be analyzed using the script [compute_SL_energy_transfer.py](compute_SL_energy_transfer.py). For each streamline in the provided set the script will:

- plot the 2D "unravelled" streamline
- plot the temperature variation along the streamline
- plot the local energy transfer due to heat conduction, axial shear work, and circumferential shear work along with the net energy transfer along the streamline, subject to some reasonable assumptions discussed in the associated publication
- Print out the total energy transfer due to each mode, integrated over the entire streamline

The code will also average the total energy transfers across all the given streamlines.

The script [compare_turbulent_vt_results.py](compare_turbulent_vt_results.py) is similar, except it contrasts specific streamlines, and is more useful for comparing different simulations (e.g. results using different turbulence models).

Both of the above scripts use the function `ReadCFXExport` to read the streamline data stored in the `.csv` files exported by CFD-Post.
