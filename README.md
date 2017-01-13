# RadTrans

RadTrans is a neutral-particle, deterministic radiation transport solver in two spatial dimensions.  This repository is a simplified/stripped-down version of a larger code that is kept in a private repository.  The main difference is that private version implements an eigenvalue solver for multiplying media.

## External Dependencies

There are a couple of external dependencies.  The main two are [HDF5](https://support.hdfgroup.org/HDF5/) and [MOAB](http://sigma.mcs.anl.gov/moab-library/) libraries.  HDF5 is the format used to store the significant pieces of data to disk (e.g., geometry, materials, solution vectors, etc.).  The MOAB library provides a lot of the tools used to represent and work with mesh-related data.

## Input

There is a nascent input file format, but I also developed quite a bit of Python to automatically build up all of the input files and data, including the primary input file.  The utility functions are bundled into a `utils` module in the `models` folder.  The `models` folder also has a simple Python script for building a homogenous slab test case.  This test script can generate input files, run RadTrans (through a shell command), perform some simple output processing and visualization, and run parametric studies.

## Output

Output from RadTrans is still evolving.  The two modes I currently rely on are an HDF5 file of the solution vector (`output.h5`) and VTK-formatted file containing the geometry and scalar flux (`output.vtk`).  The VTK file, for example, can be read into [ParaView](http://www.paraview.org/) for visualization.

More documentation to come!
