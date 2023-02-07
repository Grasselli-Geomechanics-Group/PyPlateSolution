import numpy as np
import matplotlib.pyplot as plt
import function_modules

def PYTHONSolution_example():

    ### Ray Theory Inputs ###

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
    plt.ylabel('Greens function (m/N/s)')

    plt.show()

if __name__ == '__main__':
    PYTHONSolution_example()