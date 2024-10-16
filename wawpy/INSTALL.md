Installation              {#installation}
============

The installation of the library consists of the following steps

-# Install the required/recommended packages on our computer
-# Generating the Makefiles
-# Compiling the library and examples
-# Testing the compiled software (optional)
-# Generating the documentation (optional)
-# Installing the library and examples

Make sure to install all necessary (or optional) software packages as 
described in the [README](@ref readme) file.

Ubuntu 18.04
============

Required Packages
-----------------

Required Packages          | Comment
:------------------------: | :-------
cmake                      | Needed to generate the Makefiles to compile the software
make                       | Needed to execute the Makefiles during compilation
g++                        | GNU C++ compiler
gfortran                   | GNU Fortran compiler
libgsl23                   | GNU Scientific Library
libgsl-dev                 | GNU Scientific Library (headers)
libfftw3-bin               | Fastest Fourier Transform in the West 3 library
libfftw3-dev               | Fastest Fourier Transform in the West 3 library (headers)
libboost-test1.65.1        | Boost Unit test framework library. Only needed for testing the library
libboost-test1.65-dev      | Boost Unit test framework library (headers). Only needed for testing the library 
doxygen                    | Generation of documentation. Only needed if documentation is desired.
graphviz                   | Generation of graphs in documentation. Only needed if documentation is desired.
*python/swig* packages     | Packages to compile the python 3.0 interface to the library.


Recommended Packages       | Comment
:------------------------: | :------
numpy                      | Numerical library for python
scipy                      | Scientific computation library for python
matplotlib                 | Plotting library for python

To install these packages execute the following command in a terminal window:

    $> sudo apt-get install <package-name-1> <package-name-2> ...

where <package-name-x> are the names of the packages listed in the table above.
Once the packages have been installed the installation of the wavpy library
can be done as explained in [Installation document](@ref installation).


Generating the Makefiles
------------------------

Open a terminal window and navigate to the project's main directory.

Enter into the build subdirectory:

    $> cd build

Execute cmake to generate the Makefiles for the whole project inside
the build subdirectory

    $> cmake ..

A list of messages will appear showing the progress of the Makefile generation.
If a necessary package is missing the process will stop. You will need to 
install the missing package as described in the [README](@ref readme) file.

After a successful generation of the Makefiles a set of auxiliary files and 
directories will have been created inside the build subdirectory.

Compiling the library and examples
----------------------------------

While in the same build subdirectory execute the command:

    $> make

This will compile all source files, generate the wavpy library and example programs.

Testing the compiled software
-----------------------------

A unit test suite to check the functionalities of the library can be launched
 (under development). This step is optional and will not add nor remove any
functionality to the wavpy library. It will only reveal if the library works
as expected.

To run the test suite move into the 'test' subdirectory
under the 'build' directory with
    
    $> cd  test

and execute the test suite with

    $> ctest


Generating the documentation
---------------------------

To generate the documentation of the library move back to the 'build' 
directory and then execute

    $> make doc_doxygen

This will generate the HTML documentation of the library, which can be accessed
 from any browser by opening the file

    build/docs/html/index.html


Installing the library and examples
-----------------------------------

To install the library, move again to the build directory. From there
you need to run

  $> make install

This will install the library at the installation path determined at the first
step when creating the Makefiles. In the instructions above, no installation
path has been provided, and so, the default path will be 'usr/local'.
If a different installation path is desired you should rerun all the
installation process above with the instruction

  $> cmake .. -DCMAKE_INSTALL_PREFIX=MyDesiredInstallationPath

where MyDesiredInstallationPath should be a valid path on your computer.


