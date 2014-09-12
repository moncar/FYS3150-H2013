The solar system
=============
A simulation of the solar system.

## Description
The program can take a file containing a setup for a solar
system. This file can for each object contain a name,
initial position, initial velocity and mass. It is arbitrary
how many objects you add. Examples are added for full solar
system and for Sun-Earth system.

When run, the program takes a timestep, a runtime and a file
containing a system layout. Datafiles for each object is
then stored in a path described by the constant 

```c++
OBJECTS_DATA_PATH
```

Furthermore, plotting can be done with the Python script

```bash
$ ./readObjectsData.py
```

Please note that when running this script you must be inside
the `python` directory and the `OBJECTS_DATA_PATH` folder
must be default.

Also, the directory `data/objects/` must exist. This is
where the datafiles are saved while simulation is running.

## Compiling
Before you compile, make sure you've created a
build-directory in the clone root. This directory must be
called `build`.

The make-script `Make.py` can then be run to create an
executable in the build directory. Run

```bash
$ ./Make.py
```

from clone root to compile the project.

Note that build-directory (with path) can be changed in the
file `meta.dat`. Simply edit the string after
`buildDir:<newPathHere>`.
