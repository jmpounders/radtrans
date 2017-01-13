# Installation Instructions #

 INTRODUCTION
==============

RTS is a C++ solver for radiation transport problems.  It is presently built and maintained for UNIX/Linux environments.  This set of instructions describes how to build the RTS code from a UNIX-like command terminal.  It is assumed that you have basic familiarity in working from a terminal (e.g., the commands "cd", "mkdir", "ls", etc.).

These instructions are not comprehensive in the sense that there may be alternative ways to build the code on your system.  This is especially true for the installation of the HDF5 libraries.  What follows is what has worked for me.  Good luck!


 PREREQUISITES
===============

The RTS code is mostly independent of your operating system.  There are a few portions of the code (e.g., timing objects), however, that must interface with the operating system.  These commands have been designed and tested on UNIX/Linux systems.  To run on Windows it is recommended that you use the Cygwin emulator, which can be freely downloaded from the Internet.  If you use Cygwin then also make sure to install the “gcc”, “bison”, “flex”, “cmake”, “git”, and “make” packages (all under the “Devel” grouping of packages).

The build system for RTS uses the CMake software.  This is also free software and can also be downloaded as a Cygwin package on Windows.

In principle, any compiler can be used to build the RTS code (at least at present).  The clang compiler on Mac OS X and the g++ compiler on Linux/Cygwin are known to work.

Additionally, an HDF5 library (built with the same compiler that you are using) must be installed on your system.  Instructions for installing HDF5 can be found on the HDF website.  On UNIX/Linux it is recommended you use the CMake build process.  On Windows (Cygwin) you should download the HDF tar ball, e.g. "hdf5.1.8.14.tar.gz" and unpack it using the command "tar -xzf hdf5.1.8.15.tar.gz".  In the "release_docs" directory that is created you should find a file named "INSTALL_Cygwin.txt" with detailed installation instructions.  When I went through this process (early 2015), this lead to the following procedure:
   1. Run `./configure –enabl-cxx –without-zlib`
   2. Enter “make > make.out”
   3. Enter “make > makecheck.out”
   4. Check that there are no errors in make.out and makecheck.out
   5. Enter “make install > makeinstall.out”


 INSTALLATION 
==============

Before installing you must be granted permission to the RTS Git repository.  Permission can be obtained from Justin Pounders.  To obtain this permission you will need to create a (free) account on bitbucket.org, which is the website that hosts the RTS code repository.

1. Clone the source repository
------------------------------

“Cloning” the repository means that you are downloading all of source files and directories comprising the RTS code.  You need to first decide where in your computer’s file system you want to save and build the code.  From a command terminal (in Cygwin if on Windows), change directory to wherever you want to save the code and enter the command “git clone https://username@bitbucket.org/jmpounders/rts.git” where username is your BitBucket user name.  This will create an "rts" directory with the following structure

    rts
     |---include
     |---postProcess
     |---src

Next you need to build the code.  Because the code uses the HDF5 library that you previously installed, you need to tell RTS where to find HDF5.  This is done by creating a file named "LocalConfig.cmake" in the rts/src directory.  This file should contain the following line:

  set ( hdfPath /Users/justinpounders/TPS/hdf5/HDF_Group/HDF5/1.8.14 )

Make sure to change the path (i.e., /Users/justinpounders/TPS/...) to point to the root directory of the HDF5 installation on your system.  The HDF5 root directory is the directory containing the "bin", "include", "lib", and "share" subdirectories.

2. Setup a local build directory
--------------------------------

I *strongly* recommend an out-of-source build for RTS.  This means that all of the local build-files are stored in a directory other than rts/src and rts/include.  For these instructions, I will assume that you have created a directory named "build" in the rts directory.  Your directory structure should therefore look like

    rts
     |---build
     |---include
     |---postProcess
     |---src


3. Build the code
-----------------

You are now ready to build the code.  The RTS code uses CMake to coordinate the build process.  To initialize CMake, enter the command "cmake ../src" from the "build" directory.  If this step is successful, then enter the command "make" to compile and link the code.  The executable will be placed in the "build" directory upon success.


 EXAMPLE
=========

Following is a terminal dump of the build process.  This example assumes that the HDF5 library has already been installed and the LocalConfig.cmake file has been created with the appropriate pointer to the HDF5 root directory.  Note that I have added some extraneous white space to improve readability.


    [~]% mkdir codeTest

    [~]% cd codeTest/

    [~/codeTest]% git clone https://jmpounders@bitbucket.org/jmpounders/rts.git
    Cloning into 'rts'...
    remote: Counting objects: 12104, done.
    remote: Compressing objects: 100% (8541/8541), done.
    remote: Total 12104 (delta 4075), reused 11298 (delta 3479)
    Receiving objects: 100% (12104/12104), 12.40 MiB | 853.00 KiB/s, done.
    Resolving deltas: 100% (4075/4075), done.
    Checking connectivity... done.
    Checking out files: 100% (10552/10552), done.

    [~/codeTest]% ls -al
    total 0
    drwxr-xr-x   3 justinpounders  staff   102 Mar 26 13:27 .
    drwxr-xr-x+ 46 justinpounders  staff  1564 Mar 26 13:26 ..
    drwxr-xr-x  22 justinpounders  staff   748 Mar 26 14:12 rts

    [~/codeTest]% cd rts

    [rts]% ls -al
    total 328
    drwxr-xr-x  21 justinpounders  staff     714 Mar 26 13:27 .
    drwxr-xr-x   3 justinpounders  staff     102 Mar 26 13:27 ..
    drwxr-xr-x  13 justinpounders  staff     442 Mar 26 13:27 .git
    -rw-r--r--   1 justinpounders  staff      71 Mar 26 13:27 .gitignore
    -rw-r--r--   1 justinpounders  staff  102087 Mar 26 13:27 Doxyfile
    drwxr-xr-x  29 justinpounders  staff     986 Mar 26 13:27 include
    lrwxr-xr-x   1 justinpounders  staff      18 Mar 26 13:27 input.diffusion -> input.diffusion.2g
    -rw-r--r--   1 justinpounders  staff     594 Mar 26 13:27 input.diffusion.1g
    -rw-r--r--   1 justinpounders  staff     620 Mar 26 13:27 input.diffusion.2g
    -rw-r--r--   1 justinpounders  staff     997 Mar 26 13:27 inputWriter.py
    -rwxr-xr-x   1 justinpounders  staff    2900 Mar 26 13:27 makeInput.py
    -rwxr-xr-x   1 justinpounders  staff     991 Mar 26 13:27 outputReader.py
    -rwxr-xr-x   1 justinpounders  staff     605 Mar 26 13:27 postProc.py
    drwxr-xr-x   4 justinpounders  staff     136 Mar 26 13:27 postProcess
    -rwxr-xr-x   1 justinpounders  staff    4259 Mar 26 13:27 referenceSolutions.py
    drwxr-xr-x  28 justinpounders  staff     952 Mar 26 13:27 src
    -rwxr-xr-x   1 justinpounders  staff    2528 Mar 26 13:27 test1.py
    -rwxr-xr-x   1 justinpounders  staff    2710 Mar 26 13:27 test2.py
    -rwxr-xr-x   1 justinpounders  staff    2906 Mar 26 13:27 test3.py
    -rwxr-xr-x   1 justinpounders  staff    2635 Mar 26 13:27 test4.py
    -rwxr-xr-x   1 justinpounders  staff    4818 Mar 26 13:27 test5.py

    [rts]% ls -ld */
    drwxr-xr-x  29 justinpounders  staff  986 Mar 26 13:27 include/
    -rw-r--r--   1 justinpounders  staff  620 Mar 26 13:27 input.diffusion/
    drwxr-xr-x   4 justinpounders  staff  136 Mar 26 13:27 postProcess/
    drwxr-xr-x  28 justinpounders  staff  952 Mar 26 13:27 src/

    [rts]% mkdir build

    [rts]% cd build

    [build]% cmake ../src
    -- The CXX compiler identification is AppleClang 6.0.0.6000057
    -- Check for working CXX compiler: /usr/bin/c++
    -- Check for working CXX compiler: /usr/bin/c++ -- works
    -- Detecting CXX compiler ABI info
    -- Detecting CXX compiler ABI info - done
    -- Configuring done
    -- Generating done
    -- Build files have been written to: /Users/justinpounders/codeTest/rts/build

    [build]% make
    Scanning dependencies of target Transport
    [  4%] Building CXX object CMakeFiles/Transport.dir/main.cpp.o
    [  8%] Building CXX object CMakeFiles/Transport.dir/dataset.cpp.o
    [ 12%] Building CXX object CMakeFiles/Transport.dir/dofobj.cpp.o
    [ 16%] Building CXX object CMakeFiles/Transport.dir/eigenvalue.cpp.o
    [ 20%] Building CXX object CMakeFiles/Transport.dir/element.cpp.o
    [ 25%] Building CXX object CMakeFiles/Transport.dir/fixedsource.cpp.o
    [ 29%] Building CXX object CMakeFiles/Transport.dir/hdf5interface.cpp.o
    [ 33%] Building CXX object CMakeFiles/Transport.dir/inputparser.cpp.o
    [ 37%] Building CXX object CMakeFiles/Transport.dir/logpolicy.cpp.o
    [ 41%] Building CXX object CMakeFiles/Transport.dir/log.cpp.o
    [ 45%] Building CXX object CMakeFiles/Transport.dir/legendre.cpp.o
    [ 50%] Building CXX object CMakeFiles/Transport.dir/material.cpp.o
    [ 54%] Building CXX object CMakeFiles/Transport.dir/mesh.cpp.o
    [ 58%] Building CXX object CMakeFiles/Transport.dir/node.cpp.o
    [ 62%] Building CXX object CMakeFiles/Transport.dir/output.cpp.o
    [ 66%] Building CXX object CMakeFiles/Transport.dir/outputgenerator.cpp.o
    [ 70%] Building CXX object CMakeFiles/Transport.dir/perfstats.cpp.o
    [ 75%] Building CXX object CMakeFiles/Transport.dir/solutionmanager.cpp.o
    [ 79%] Building CXX object CMakeFiles/Transport.dir/solverbase.cpp.o
    [ 83%] Building CXX object CMakeFiles/Transport.dir/solverdiffusion.cpp.o
    [ 87%] Building CXX object CMakeFiles/Transport.dir/solversweepsi.cpp.o
    [ 91%] Building CXX object CMakeFiles/Transport.dir/timing.cpp.o
    [ 95%] Building CXX object CMakeFiles/Transport.dir/transient.cpp.o
    [100%] Building CXX object CMakeFiles/Transport.dir/transportproblem.cpp.o
    Linking CXX executable Transport
    [100%] Built target Transport

    [build]% ls -al
    total 3856
    drwxr-xr-x   8 justinpounders  staff      272 Mar 26 14:14 .
    drwxr-xr-x  22 justinpounders  staff      748 Mar 26 14:12 ..
    -rw-r--r--@  1 justinpounders  staff     6148 Mar 26 14:13 .DS_Store
    -rw-r--r--   1 justinpounders  staff    11274 Mar 26 14:13 CMakeCache.txt
    drwxr-xr-x  12 justinpounders  staff      408 Mar 26 14:14 CMakeFiles
    -rw-r--r--   1 justinpounders  staff    22008 Mar 26 14:13 Makefile
    -rwxr-xr-x   1 justinpounders  staff  1923968 Mar 26 14:14 Transport
    -rw-r--r--   1 justinpounders  staff     1270 Mar 26 14:13 cmake_install.cmake

---------------------------------------------------------------------
> Author: Justin Pounders
> 
> Last updated: March 26, 2015

