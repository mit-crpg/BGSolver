BGSolver
========
4/15/2013

--------------------------
To install BGSolver v1.03:

1. Navigate MATLAB to the /Installer directory.
2. Execute the following script:

installBGSolver 

----------------------------
To uninstall BGSolver v1.03:

0. This will only work if BGSolver v1.03 has been installed, which means its subdirectories are on MATLAB path.
1. Execute the following script:

uninstallBGSolver

------------------------------
Note on .mex file compilation:

BGSolver v1.03 contains two .c files that can be compiled into .mex. /BG_Processor/evalxdot.c function, in particular, greatly improves the code's performance if it is compiled into a .mex file. Note, that the code is fully functional WITHOUT compiling any of the .c files as well, but ensuring that .mex files are present may improve the code performance.

To compile the .mex files:

0. This will only work with a correctly configured compiler. Compilers are configured through "mex -setup" command. To compile /BG_Processor/evalxdot.c, the compiler must support OpenMP 2.0. Note, that Microsoft Software Development Kit (SDK) 7.1 does NOT support OpenMP 2.0, and therefore will not be able to compile /BG_Processor/evalxdot.c.
1. Ensure that BGSolver v1.03 has been installed. See the instructions above for how to install BGSolver v1.03.
2. Execute the following script:

compileBGSolver

3. If the script does not report any errors or warnings, the compilation was successful.

The code is distributed with .mexw64 (64-bit Windows) executables, so there is no need to recompile the .c files on a 64-bit Windows platform.

To familiarize one's self with how the MATLAB Executable files (.mex) work:
http://www.mathworks.com/help/matlab/create-mex-files.html

To compile a .mex file, a correctly configured compiler is necessary. See:
http://www.mathworks.com/help/matlab/matlab_external/building-mex-files.html

Compiled executables for /Utilities/GetFullPath/GetFullPath.c can also be found at:
http://www.n-simon.de/mex
