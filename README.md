This is a PYTHON program to compute the Green's function of an infinite plate. It was originally written by [Nelson Hsu (1985)](https://nvlpubs.nist.gov/nistpubs/Legacy/IR/nbsir85-3234.pdf) mainly for its application to calibrate acoustic emission systems and sensors.

This PYTHON code is translated, by Edouard Kravchinsky, from the MATLAB code of Greg McLaskey. The MATLAB code was translated, by [Greg McLaskey (2010)](https://doi.org/10.1121/1.3466847), from a FORTRAN code written by [Nelson Hsu (1985)](https://nvlpubs.nist.gov/nistpubs/Legacy/IR/nbsir85-3234.pdf) at the National Bureau of Standards. There are a few known bugs in this code.

Notes: Greg McLaskey

1)	The units in the original FORTRAN code were not entirely correct. This was fixed in the MATLAB version. Now the output is in units of meters, as described in the comments of the MATLAB function, and reprinted below.
2)	The code does not work if the source is directly opposite the sensor (i.e. ```xdist = 0``` in the MATLAB version, or ```source-detector horizontal distance = 0``` in the FORTRAN version). To circumvent this problem, a very small value for ```xdist``` can be substituted without any loss of accuracy. 
3)	I believe the code will output non-zero values for all indexes, yet, according to [L. R. Johnson (1974)](https://doi.org/10.1111/j.1365-246X.1974.tb02446.x) the only nonzero indexes should be: 11 31 22 13, 33, 111, 311, 221, 131, 331, 212, 322, 232, 113, 313, 223, 133, 333. All other indexes should be identically zero. To be consistent with the notation of the same paper, the polarity of some of the Green’s functions should be switched. For example, 31 311, 322, 333… etc. 
4)	For buried source problems (i.e. when ``zdist`` is less than ```h```) some of the later arriving reflections are incorrect!

Notes: Edouard Kravchinsky

5)	As mentioned in the original FORTRAN code / paper, if the sensor and detector are far apart or the maximum time is large, the number of rays arriving at approximately the same time may be too large to compute in a reasonable time and the accumulated error may grow.

Use the script ```python PYTHONPlateSolution_example.py``` for an example on how to use the code or the Jupyter example ```PYTHONPlateSolution_example.ipynb```. All of the functions are in ```function_modules.py``` with the main function to call being ```function_modules.PYTHONPlateSolutionFunc```.

```python 
def PYTHONPlateSolutionFunc(Cs, Cp, rho, h, xdist, zdist, T, npoint, INDEX, printOption):
    # [output RAYTIME]= PYTHONPlateSolutionFunc(Cs,Cp,rho,h,xdist,zdist,T,npoint,INDEX,printOption)
    # PYTHONPlateSolution -
    # Written by Nelson Hsu from the National Bureau of Standards (NBS)
    # Translated into MATLAB by Gregory McLaskey in 2009
    # Translated into Python by Edouard Kravchinsky in 2022
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
```

An example:

```python
import numpy as np
import matplotlib.pyplot as plt
import function_modules

def PYTHONSolution_example():

    Cs = 3230  # Swave speed (m/s)
    Cp = 5900  # Pwave speed (m/s)
    h = 0.05  # plate thickness (m)
    xdist = 0.005  # horizontal distance from source to sensor (m)
    zdist = h  # vertical dist from source to sensor (should be 0 or h) (m)
    t_tot = 30*10**-6  # total runtime (s)
    T = 1e-7  # the sampling period (s)
    npoint = int(t_tot / T)  # the total number of time points to be computed
    INDEX = 33  # the greens function subscript index (i.e. G33 G331 etc)
    rho = 7850  # density of medium (kg/m3)
    printOption = 1

    [output, RAYTIME] = function_modules.PYTHONPlateSolutionFunc(Cs, Cp, rho, h, xdist, zdist, T, npoint, INDEX,
                                                                 printOption)

    print('Done')

    plt.figure(1)

    plt.subplot(211)
    plt.plot(np.arange(1, npoint + 1) * T, output)
    plt.xlabel('time (s)')
    plt.ylabel('displacement (m) \n due to 1 N step force')
    plt.title('x = {:g} z = {:g}'.format(xdist, zdist))

    plt.subplot(212)
    plt.plot(np.arange(1, npoint) * T, np.diff(output, axis=0) / T)
    plt.xlabel('time (s)')
    plt.ylabel('Greens function (m/(N*s))')

    plt.show()

if __name__ == '__main__':
    PYTHONSolution_example()
```
