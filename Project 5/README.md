Open cluster simulation
=============

This project simulates a cold spherical collapse of an open star
cluster using the Leapfrog and RK4 integration schemes.

A more detailed documentation is provided in the
[wiki](https://github.com/benedicteb/FYS3150-Prj05/wiki).

## Compiling

After clone project can be compiled by running the `Make.py` python
script from the root of the project. This will create a build
directory if not already present and place the executable here.

```bash
$ ./Make.py
```

All tags relating to compilation can be set in the `meta.dat` file
which has an easy-to-understand syntax. This will create a build
directory with the executable, which can then be run as described in
the [wiki](https://github.com/benedebr/FYS3150-Prj05/wiki/Executable).
