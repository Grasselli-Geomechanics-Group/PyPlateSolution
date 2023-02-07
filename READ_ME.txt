Use the script “PYTHONPlateSolution_example.m” for an example of how to use the code. 

The main part is the function 
“PYTHONPlateSolutionFunc.m”
Which also includes comments and notes within the function for additional reference.

Some notes for “PYTHONPlateSolutionFunc”                				January, 2023

This PYTHON code is translated, by Edouard Kravchinsky, from the MATLAB code of Greg McLaskey. 
The MATLAB code was translated, by Greg McLaskey, from a FORTRAN code written by Nelson Hsu at the NBS. There are a few known bugs in this code. 

Notes: Greg McLaskey
1)	The units in the original FORTRAN code were not entirely correct. This was fixed in the MATLAB version. Now the output is in units of meters, as described in the comments of the MATLAB function, and reprinted below
2)	The code does not work if the source is directly opposite the sensor (i.e. xdist = 0 in the MATLAB version, or “source-detector horizontal distance” = 0 in the FORTRAN version). To circumvent this problem, a very small value for xdist can be substituted without any loss of accuracy. 
3)	I believe the code will output non-zero values for all indexes, yet, according to L. R. Johnson 1974 (Geophys J. R. astr. Soc 37, 99-131) the only nonzero indexes should be: 11 31 22 13, 33, 111, 311, 221, 131, 331, 212, 322, 232, 113, 313, 223, 133, 333. All other indexes should be identically zero. To be consistent with the notation of the same paper, the polarity of some of the Green’s functions should be switched. For example, 31 311, 322, 333… etc. 
4)	For buried source problems (i.e. when zdist is less than h) some of the later arriving reflections are incorrect!

Notes: Edouard Kravchinsky
5)	As mentioned in original FORTRAN code / paper, if the sensor and detector are far apart or the maximum time is large, the number of rays arriving at approximately the same time may be too large to compute in a reasonable time and the accumulated error may grow.

def PYTHONPlateSolutionFunc(Cs, Cp, rho, h, xdist, zdist, T, npoint, INDEX, printOption):
    # [output RAYTIME]= PYTHONPlateSolutionFunc(Cs,Cp,rho,h,xdist,zdist,T,npoint,INDEX,printOption)
    # PYTHONPlateSolution -
    # Written by Nelson Hsu from the National Bureau of Standards (NBS)
    # Translated into MATLAB by Gregory McLaskey in 2009
    # Translated into Python from MATLAB version by Edouard Kravchinsky in 2022
    #
    # this code computes the Green's functions for an infinite elastic plate
    #
    # OUTPUTS ----
    # output = the displacements (in meters) that would be produced by:
    # a) for the case of INDEX == 11, 12, 13, 21, 22, 23, 31, 32, 33.
    #     a step force of 1 Newton
    # b) for the case if INDEX = 111, 112, etc.
    #     a step moment of 1 Newton meter
    # RAYTIME = a matrix containing wave arrival information
    #
    # INPUTS -----
    # Cs = Swave speed (m/s)
    # Cp = Pwave speed (m/s)
    # rho = density of the material - only needed for the amplitude scale(kq/m^3)
    # h = plate thickness (m)
    # xdist horizontal distance from source to sensor (m)
    # zdist vertical dist from source to sensor (m)
    # T = the sampling period (s)
    # npoint = the total number of time points to be computed
    # INDEX = the greens function subscript index (i.e. G33 G331 etc)
    # printOption ---
    #   if 0: the program will compute the Green's function and RAYTIME
    #   if 1: the program will also display RAYTIME and print an update every 5
    #           time steps.
    #   if 2: the program will only compute RAYTIME

An example:

[output RAYTIME]=MATLABPlateSolutionFunc(3230,5900,7850,0.05,0.005,0.05,1e-7,300,33,1)
plt.figure(1)
plt.subplot(211)
plt.plot(np.arange(1, npoint + 1) * T, output)
plt.xlabel('time (s)')
plt.ylabel('displacement (m) \n due to 1 N step force')
plt.title('x = {:g} z = {:g}'.format(xdist, zdist))

plt.subplot(212)
plt.plot(np.arange(1, npoint) * T, np.diff(output, axis=0) / T)
plt.xlabel('time (s)')
plt.ylabel('Greens function (m/N/s)')
