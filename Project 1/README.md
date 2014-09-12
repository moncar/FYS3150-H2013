# Introduction
This repository contains code for project 1 in computational
physics, autumn 2013.

# Usage
There are mostly two kinds of files.

## C++ files
C++ is used as the main language in this project. All computations
are written in C++. Many relevant functions are found in the
"main.cpp" file. Yet some are also separated, for example the
function containing the tridiagonal matrix algorithm.

The program can be run via a commandline interface. The relevant
flags are:
* -b: For finding the solution to one specific problem. The "-n"
  flag must also be given.
* -c: For performing an error analysis. All three flags "-n",
  "-logA" and "-logB" must be given.
* -d: For calculating a solution to the problem via the
  LU-decomposition. The flag "-n" must also be given.
* -s: For finding the biggest relative error for a specific
  solution. The solution is found via the tridiagonal matrix
  algorithm. The flag "-n" must also be given.
* -n <number of steps>: Gives the number of steps for the
  calculation. If the flag "-d" is given, this is how many steps
  for the error analysis, not the computation.
* -logA <lower log(n)>: The lower limit for error analysis. The
  logarithm of the number of steps.
* -logB <upper log(n)>: The upper limit for error analysis. The
  logarithm of the number of steps.

## Python-files
I used Python to create a convenient make-file. Also, there are
several files used for plotting and reading data.

# Dependencies
The project depends on the Armadillo library.
