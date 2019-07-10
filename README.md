# A library for HPC reservoir simulations

## Introduction

This library relies mainly on the power of existing OpenFOAM functionality to build highly capable solvers for (isotherm processes at 
the moment) fluid flow in oil and gas reservoirs.

It should compile seamlessly on most OpenFOAM forks with versions 4.0+ although the native dynamic mesh capabilities of OF6 solvers
are not used (But planned).

> The main technical focus goes to (In This Order):

> - Comprehensive well control models

> - Parallel Execution Efficiency, both on POSIX systems and Windows.

> - Windows-friendly Implementations as it's probable someone
>   will eventually try to run this on a Windows machine.

By design, the library doesn't care so much for user convenience (GUI-wise). Integration with some GUI packages ( Think: 
[Sim-Flow](https://sim-flow.com) ) may be of interest at a later point. Users are expected to plot things using
*gnuplot* or python's *matplotlib* until a capable framework of **functionObjects** is developed.

## Main purpose of this repository

This repository holds the (to-be) open source code for the framework; At the moment this is a private repo; only invited
developers are allowed to see its content. There is:

- A [wiki](https://bitbucket.org/FoamScience/reservoir-simulator/wiki) that documents code development advances.
- An [Issues](https://bitbucket.org/FoamScience/reservoir-simulator/issues?status=new&status=open) page to discuss design decisions
  feature requests/implementations.

The most recent advances in code developments (Including new features and Bug fixes) are tracked in a dedicated 
[ChangeLog](https://bitbucket.org/FoamScience/reservoir-simulator/wiki/ChangeLog) page.

## Installation and usage

### Requirements

- A working installation of (a recent version 4.0+) OpenFOAM.
- Meshing is almost always done using `blockMesh` (A default utility) and [cfMesh](http://www.cfmesh.com/cfMesh)
- Python is occasionally used to perform/automate minor tasks.
- `gnuplot` is often used for plotting/monitoring.
- The post processor of choice is Kitware's **ParaView** (Open Source, python-scriptable).
- Parameter variation and optimization cases should be run with Sandia Labs' open source 
  [DAKOTA](https://dakota.sandia.gov/) package.

### Install on Linux/Unix

First, you need to download the source code of the framework:

```
git clone https://FoamScience@bitbucket.org/FoamScience/reservoir-simulator.git
```

Or just download the .zip file from this page.

Installation on Unix-like systems (assuming a sourced, working installation of OpenFOAM) should go seamlessly if one runs the 
following command in a shell inside the framework's directory:

```bash
./Allwmake
```

The previous command compiles library parts and solvers/utilities and installs them into the user's dedicated directories;
namely:

- `$FOAM_USER_LIBBIN` for library files (.so)
- `$FOAM_USER_APPBIN` for solver and utility binaries.

Use the `echo` command to figure out where these are:

```bash
echo $FOAM_USER_LIBBIN
```

To verify that everything works as expected, cop one of the tutorials and run its `run` script. Mainly
the run script will consist of a series of commands to create the mesh, set initial/boundary fields, run the appropriate 
solver and then post-process results.


### Install on Windows

Porting OpenFOAM to Windows is, relatively, a simple task. Many attempts have created multiple forks of OpenFOAM with
native Windows executable. The recommended package is [blueCFD-Core](http://bluecfd.github.io/Core/Downloads/).

On Windows, `cfMesh` probably won't compile, so all cases using it to create the mesh should carry their mesh files around.

For compilation of Windows, blueCFD uses MSys2 as a development infrastructure; Resulting binaries are of course native to Windows
and can be used independently. So, after installing blueCFD on a Windows machine, the command `./Allwmake` in the framework's 
root directory would work as expected.

#### Troubleshooting Windows compilation

Depending on your system's configuration, compiling with `Allwmake` may fail.
If this happens, here are the steps to follow:

- First, clean any compilation attempts using the command `./Allwclean`
- Use a bash for loop to compile each library part on its own, and then compile solvers/utilities separately:

```bash
# Store current directory (should be framework root directory)
curr=$PWD

# Loop through libs parts and compile them
for direc in libs/{relativePermeabilityModels,capillaryPressureModels,transportModels,wellModels,finiteVolume/fields/fvPatchFields};
do
        cd $direc
        wmake -j libso
        cd $curr
done

# Compile solvers
wmake all solvers
# Compile utilities
wmake all utilities
```

## Library structure

We followed the structure of the original OpenFOAM code:

- **doc**: Configuration for Doxygen documentation, just run the command `doxygen Doxyfile` to generate HTML 
  documentation of all classes/namespaces in the framework.
- **libs**: This is where Source/Header files for library parts are stored.
- **solvers**: Code for individual solvers (impesFoam, blackOilFoam).
- **tutorials**: Verification and usage-proof cases for different solvers.
- **utilities**: Code for pre-processing utilities (Grid conversion ...).


## Extend the framework for your own needs

It's recommended that you keep the original source files intact. OpenFOAM's dynamic library loading is a fantastic option to 
have around at times like this.

Let's assume one needs to add a new relative permeability model to the toolkit, we'll call it "myNewModel". 
All you have to do is to:

- Copy the code of similar model (say, `krBrooksCorey`) to wherever you like:
```bash
cp -r $PATH-TO-ORIGINAL-LIBRARY/libs/relativePermeabilityModels/twPhaseModels/BrooksCorey myNewModel
```
- Rename all files from "krBrooksCorey" to "myNewModel".

- Change every appearance of "BrooksCorey" to "myNewModel" in copied files:
```bash
sed -i 's/BrooksCorey/myNewModel/g' myNewModel/*
```

- The important thing to keep in mind is that every new model should publicly inherit from `relativePermeabilityModel` class
  and thus virtually define a `correct()` member method that calculates relative permeabilities with their derivatives
  as functions of water saturation.
  
- To compile the new model(s), a `Make` directory (consisting of two files: "files,options") should be created.

- The `Make/files` file must point to the compiled source files (.C) and library installation path:
```bash
# File: Make/files
myNewModel/myNewModel.C
LIB = $(FOAM_USER_LIBBIN)/libmyCustomKrModels
```

- The `Make/options` file must show which existing libraries to use and point to included header files (.H)
```bash
# File: Make/options
EXE_INC = \
    -I$(LIB_SRC)/finiteVolume/lnInclude \
	-I/path/to/original/libraries/libs/relativePermeabilityModels/lnInclude
	
LIB_LIBS = \
    -lfiniteVolume \
	-L$(FOAM_USER_LIBBIN) \
	-lrelativePermeabilityModels
```

Running the command `wmake libso` now in the parent directory of `myNewModel` directory should compile and install the new library.
To use it, it's sufficient to introduce the keyword:
```cpp
libs ("libmyNewModel.so");
```
At the end of `system/controlDict` in your OpenFOAM case.

> You should be aware that solvers are compiled with a general interface to boundary conditions; so, to use custom boundary
> conditions like those defined in the `"libreservoirBCs.so"` library, you have to include it in `controlDict`.
