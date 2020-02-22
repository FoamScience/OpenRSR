> This work is not affiliated with OpenCFD, the owner of the OpenFOAM trademark, in any
> way.
> If you need old versions of the code, they are now available in 
> [this Bitbucket repo](https://bitbucket.org/FoamScience/reservoir-simulator/src) 
> although the code is messy.

[![CodeFactor](https://www.codefactor.io/repository/github/foamscience/openrsr/badge/dev)](https://www.codefactor.io/repository/github/foamscience/openrsr/overview/dev) [TravisBuild][https://travis-ci.com/FoamScience/OpenRSR.svg?branch=dev] [![Codacy Badge](https://api.codacy.com/project/badge/Grade/fba58294886e494fa636370dd24cc79e)](https://www.codacy.com/manual/elwardifadeli/OpenRSR?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=FoamScience/OpenRSR&amp;utm_campaign=Badge_Grade)

# Open Reservoir Simulation Research tool

The toolbox makes use of latest GPLed code in Foam-Extend repositories to build capable
(but not yet efficient enough) solvers for BlackOil equations in isotherm porous media.
It's thanks to the OpenFOAM developpers, `henrus`, `bgscaid` and the other 
[contributors](https://github.com/Unofficial-Extend-Project-Mirror/foam-extend-foam-extend-4.0/graphs/contributors)
to the project that this tool exists. We hope this project could contribute to theirs at 
some point.

The use of the extend branch instead of mainline OpenFOAM is imposed upon us because we
are in dire need of the coupling functionality which mainline OpenFOAM lacks.

> This is version 0.1, so you can expect things to break easily.
> - I'll try to keep the master branch clean and functional
> - The whole library is configured to compile in Debugging mode which naturally causes performance issues,
>   you can compile it in Optimized mode by tweaking your `wmake` configuration.
> - Parallel simulations will (virtually) work as long as you don't use `labelToCell` to
>   specify well cells, or split a single well's cells to many processor. Otherwise, it's
>   a bug which you should report. There is no MPI code in well models yet.

## Introduction

> The current state of the project allows only for 2-phase simulations, I know, silly but
> keep reading please

I developed this library originally for my Master Thesis project to show the potential of 
some Open Source software (OpenFOAM) in Petroleum Reservoir Simulation. So, the code is fairly
simple, with no extensive use of templates/advanced C++ stuff: Point is, anyone who can
read, can read it. This would change in the future though.

The whole purpose of the library is to provide open-access for reservoir engineers to go 
beyond classic IMPES simulations and research their own ideas easily, with no third party 
(Giant Oil & Gas Companies) telling them what models to use and what's not. People should 
always have access to the source code of the software they use, especially when scientific
tasks are involved.

> The most simple class is the capillary pressure model: 
> [none](https://github.com/FoamScience/OpenRSR/blob/master/libs/capillaryPressureModels/none/none.H), 
> so it's good for new comers to go 
> there first; The most complex ones are the well classes for now.

The coupled solver is probably the best choice in many cases, so stick with it until a
fully implicit solver (Newton-Raphson) solver is developed; or develop your own. :point_left:

## Installation and Usage Instructions

### Compiling the library

To compile all library parts to binary shared objects and
executables
(`Allwmake` script should have execution permissions):

```sh
./Allwmake
```

Take a look at `Make/options` for each library department to change 
compiled classes and library name, ..., etc. I assume you're familiar with this kind of
situations.

### Use in own solvers

One can dynamically link some libraries in `controlDict` of a simulation 
case by adding the following entry:

```cpp
libs
(
    libmyOwnKrModels.so
);
```

Which should allow for use of custom relative permeability models (for example) 
without modifying existing solvers (Like UDFs in commercial software).

If a custom solver is to instantiate an object with the help
of some library class, the library must be statically linked (At compile time).
IMPES solver shows the simplest example:

```cpp
// Include parent virtual class header file for all two-phase Kr models
#include "relativePermeabilityModel.H"

// Now we can construct the model using static New method
// And read its parameters from constant/transportProperties
autoPtr<relativePermeabilityModel> krModel = 
        relativePermeabilityModel::New("krModel", transportProperties, Sw);

krModel->correct(); // Most models have a correct method which 
                    // calculates or updates model results.

// We can now retrieve values (or refs to values) of interest from the class
volScalarField& krw = krModel->phase1Kr();
// which would return a ref to Krw in current processor region.
```

And make sure the solver's `Make/options` statically links the library in
question:

```sh
EXE_INC = \
    # Some lnIncludes the solver already has
    -I../../libs/permeabilityModels/lnInclude # <-- path to library sources

EXE_LIBS = \
    # Some libraries the solver already loads
    -L$(FOAM_USER_LIBBIN) \ # <-- if the library is located here (it should)
    -lrelativePermeabilityModels  # <-- load the actual library
    # This assumes the shared object is called librelativePermeabilityModels.so
```

> Changing library code will take effect immediately (in most cases). There
> is no need to re-compile solvers.

### Compiling native Windows executable

Never tried it, and probably wouldn't try. But I would be happy to see someone try it.
A hit: follow the same procedure described in the original Foam-Extend INSTALL file.

## Library and Solver Verification

The following benchmark cases belong to my master thesis

- Buckley-Levrett Theory (many scenarios under different circumstances).
- SPE10 Dataset A
- A Reduced model for GRDC reservoir (described in the master thesis) will be published 
  in the repos.


[![](https://codescene.io/projects/6391/status.svg) Get more details at **codescene.io**.](https://codescene.io/projects/6391/jobs/latest-successful/results)
