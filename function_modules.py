import numpy as np
import math
import cmath
import time


def COEFF(N1, N2, JS, J2):
    kkkk = N2 + KDELTA(2, J2) * KDELTA(1, JS) - KDELTA(1, J2) * KDELTA(2, JS)
    if N1 <= kkkk:
        KU = int(N1)
    else:
        KU = int(kkkk)

    KL = int(KDELTA(2, J2) + KDELTA(1, J2) * KDELTA(1, JS))
    KD = (KDELTA(0, N1) * KDELTA(2, J2) + KDELTA(0, N2) * KDELTA(1, J2)) * KDELTA(JS, J2)

    if (KU + 2 - KL) != np.floor(KU + 2 - KL):
        print(N1, N2, KU, KL)

    KSUM = np.zeros([1, KU + 2 - KL + 1])
    # print(KL,KU,KU+2-KL+1)
    # print(KSUM)
    for NM in range(KL, KU + 1):
        KSUM[0, NM + 1] = IBINO(N1 - KDELTA(2, J2), NM - KDELTA(2, J2)) * IBINO(N2 - KDELTA(1, J2), NM - KDELTA(1, JS))
    KSUM = KSUM[:, 1:]

    if np.size(KSUM) == 0:
        KSUM = np.zeros([1, 1])
    # a = np.shape(KSUM)
    # for i in range(0,len(a)): # Not sure MATLAB length(a)==2???
    #    if a[i] == 0:
    #        KSUM = np.zeros([1,1])
    #        break

    return KD, KU, KL, KSUM


def COMMON(DS, Y, BLK0, BLK1):
    # note that DS = PHI(U)
    AA = BLK0['alpha'] ** 2
    DN1 = BLK1['AN1']
    DN2 = BLK1['AN2']
    YAA = np.copy(AA)
    YA = BLK0['alpha']
    DXD = BLK0['XD']

    ETA1 = DS / DXD
    ETA2SQ = 1 - ETA1 ** 2
    YETA1 = ETA1.real
    YETA2S = ETA2SQ.real
    YETA2 = cmath.sqrt(YETA2S)

    YD = DS.real
    YY = Y ** 2  # this comes from the YROOT
    YP = cmath.sqrt(YAA * YY - 1)
    YQ = cmath.sqrt(YY - 1)
    YPSQP1 = YAA * YY
    YQSQP1 = YY
    YQSQM1 = YQ ** 2 - 1
    YPSQM1 = YP ** 2 - 1

    YSQRP1 = Y * YA
    YSQRQ1 = Y
    YDELTA = YQSQM1 * YQSQM1 + 4 * YP * YQ
    YM = YQSQM1 * YQSQM1 - 4 * YP * YQ
    YR0 = YM / YDELTA
    YRP = -(4 * YP * YQSQM1) / (YA * YDELTA)
    YRM = (YR0 ** 2 - 1) / YRP
    YDN1 = DN1.real
    YDN2 = DN2.real
    YDPHI = 1 / (YD + YDN1 / YP + YDN2 / YQ)

    # re-save the variables into BLK2
    BLK2 = {}
    BLK2.update({'YD': YD, 'Y': Y, 'YP': YP, 'YQ': YQ, 'YR0': YR0, 'YRP': YRP, 'YRM': YRM, 'YDELTA': YDELTA,
                 'YPSQP1': YPSQP1, 'YQSQP1': YQSQP1, 'YPSQM1': YPSQM1, 'YQSQM1': YQSQM1, 'YSQRP1': YSQRP1,
                 'YSQRQ1': YSQRQ1, 'YY': YY, 'YA': YA, 'YAA': YAA, 'YETA1': YETA1, 'YETA2': YETA2, 'YDPHI': YDPHI})

    return BLK2


def DCRTNI(XST, EPS, IEND, T, R, alpha, C1, C2):
    # compute the root of the equation PHI based on an initial guess (which is
    # damn good)
    #
    # OUTPUTS: X = root
    #          F = function value at root (should be zero)
    #          DERF = function derivative at the root
    #          IER is an error function = 1 if convergence never met
    #                                   = 2 if slope of zero encountered
    #                                   = 0 if its all good
    #
    #
    # XST is also known as YOLD and is the initial guess
    # EPS = 1e-12 # error upper bound
    # IEND = 100  # limit on iterations #
    #
    # we are trying to evaluate an integral from 0 to -DXD. 15 points are taken
    # along this interval... DXD = XD = xdist/h. It could be something very
    # very small or something on the order of 10 R is the variable PHI(U) from the GL15T function which is the
    # dependent variable of the integral. U is the point along the integration
    # path. For XD very small (sensor below the source) PHI(U) = some very
    # small negative number.
    #
    # here T = TIMEI

    # most of the inputs come from the PHIEQ0 program
    # C1 C2, these are number of P and S trips, A = alpha, Z = the dependent
    # variable of PHIO, which is solved for in DCRTNI. T = the current
    # time normalized by h/cs
    # R is the variable U from GL15T function which is the
    # dependent variable of the integral.

    IER = 1  # this means no convergence
    X = XST  # start at the initial value
    TOL = XST
    F = PHIO(TOL, T, R, C1, C2, alpha)
    DERF = DPHIO(TOL, T, C1, C2, alpha)
    DTOLF = 100 * EPS

    for i in range(1, IEND + 1):
        if abs(F) != 0:  # if the function is not yet zero
            if abs(DERF) == 0:
                print('the slope is zero')
                IER = 2
                break
            DX = F / DERF
            X = X - DX
            TOL = X
            F = PHIO(TOL, T, R, C1, C2, alpha)
            DERF = DPHIO(TOL, T, C1, C2, alpha)

            # test on satisfactory accuracy
            DTOL = EPS
            if abs(X) > 1:  # if X is larger than 1, then increase the tolerance?
                DTOL = DTOL * abs(X)
            if abs(DX) - DTOL <= 0:  # if the dependent variable Z only moved a tiny amount
                if (abs(F) - DTOLF) <= 0:  # if the function is basically zero here
                    IER = 0  # (debug) it seems that whenever this function has to iterate, it doesn't ever find the root
                    break  # then the root has been found
                    # otherwise keep iterating!
        else:  # if the function is already zero
            IER = 0
            break

    return X, F, DERF, IER


def DEQ0(T, A, C1, C2):
    #     COMPUTE THE ROOT BY SOLVING QUARTIC EQ.

    #       YF(Y)=-DCMPLX(DT,0.D0)*Y-DCMPLX(DK1,0.D0)*
    #      1       CDSQRT(Y*Y*DCMPLX(DA*DA,0.D0)-DCMPLX(1.D0,0.D0))
    #      2      -DCMPLX(DK2,0.D0)*CDSQRT(Y*Y-DCMPLX(1.D0,0.D0))

    DT = T
    DA = A
    DK1 = C1
    DK2 = C2

    DT2 = DT * DT
    DK12 = DK1 * DK1
    DK22 = DK2 * DK2
    DAA = DA * DA

    B = [((DT2 - DK12 * DAA - DK22) ** 2 - 4 * DK12 * DK22 * DAA),
         2 * ((DT2 - DK12 * DAA - DK22) * (DK12 + DK22) + 2 * DK12 * DK22 * (DAA + 1)),
         (DK12 + DK22) ** 2 - 4 * DK12 * DK22]

    rootsAns = np.roots(B)

    NOROOT = max(rootsAns.shape)

    # if NOROOT==1
    #    rootsAns = [rootsAns 0]
    # elseif NOROOT == 0
    #     rootsAns = [0 0]
    # 

    Y1 = cmath.sqrt(rootsAns[0])
    Y2 = cmath.sqrt(rootsAns[1])
    if np.imag(Y1) > 0:
        Y1 = -Y1

    if np.imag(Y2) > 0:
        Y2 = -Y2

    if abs(PHIO(Y1, T, 0, C1, C2, A)) < abs(PHIO(Y2, T, 0, C1, C2, A)):
        Y = Y1
    else:
        Y = Y2

    if np.imag(Y) >= 0:
        Y = np.conj(Y)

    return Y, NOROOT


def DGLQ1(DINIT, DEND, EPS, NINT, RST, W, NMAX, TIMEI, BLK0, BLK1, itloop):
    # [R E FMIN FMAX KF IFLAG] =
    # DGLQ1(DINIT,DEND,EPS,NINT,RST,W,NMAX,TIMEI,BLK0,BLK1)
    #
    # notes: in the original FORTRAN code it seemed as if the criterion for subdividing the
    # integral into more than two parts was never reached. In other words, the
    # tolerance EB was so high that E>EB was never true (see E and EB definitions below)
    #
    # this version was updated on Oct 20, 2009 to incorporate surface waves.
    # This adds the parameter code. O = not a surface wave,
    #                               1 = DQR1
    #                               2 = DQR2
    #                               3 = DQR3

    A = DINIT
    B = DEND
    SIGN = -1
    oldYROOT = 0

    NMAX = 10
    IFLAG = float('NaN')

    # maybe he calls UFLOW and EPMACH
    EPMACH = 16 ** (-13)  # program seems to be very insensitive to these variables.
    UFLOW = 16e-64

    goto20 = 0
    goto15 = 0

    if A != B:  # integrate over no distance
        if not RST:  # this should be false on the initial call ALWAYS FALSE IN MATLAB??
            KF = 0
            if not (EPS < 0 or NMAX <= 1 or NINT < 0):  # no improper use of inputs
                if NINT == 1:  # it should be equal to 1
                    W = np.zeros([NMAX, 6], dtype=np.complex128)
                    W[0, :] = [A, A + (B - A) / 2, 0, 0, A, B]
                    W[1, :] = [A + (B - A) / 2, B, 0, 0, A, B]
                    NINT = 2
                else:
                    print('found my way to line 21 of DGLQ1')
                    if W[0, 0] != A and W[NINT, 0] != B:
                        IFLAG = 8
                        goto20 = 1
                    else:
                        W[0, 4] = A
                        for i in range(1, NINT + 1):
                            W[i - 1, 1] = W[i, 0]
                            W[i - 1, 4] = W[i - 1, 0]
                            W[i - 1, 5] = W[i - 1, 1]
            else:
                IFLAG = 7
                print('improper use of inputs')
                goto20 = 1

        else:  # RST  is true
            if IFLAG < 3:
                goto15 = 1
            else:
                goto20 = 1
    else:  # in this case A == B and integrate over no distance
        R = 0
        E = 0
        NINT = 0
        IFLAG = 0
        KF = 1
        [FMIN, YROOT] = FFUNC(PHI(A, B, A), BLK0, BLK1, oldYROOT, itloop)
        FMAX = np.copy(FMIN)
        goto20 = 1

    if not goto20 and not goto15:
        DEBUG = 0
        IFLAG = 0
        IROFF = 0
        RABS = 0

        for i in range(1,
                       NINT + 1):  # # of intervals the integral was divided up into (should be 2 intervals initially)
            [Rtemp, AEtemp, RAB, RAV, FMN, FMX, YROOT] = GL15T(W[i - 1, 0], W[i - 1, 1], W[i - 1, 4], W[i - 1, 5],
                                                               TIMEI, BLK0, BLK1, oldYROOT, itloop, i)
            oldYROOT = np.copy(YROOT)
            W[i - 1, 3] = Rtemp
            W[i - 1, 2] = AEtemp
            KF = KF + 15
            if i == 1:
                R = np.copy(Rtemp)
                E = np.copy(AEtemp)
                RABS = RABS + RAB
                FMIN = np.copy(FMN)
                FMAX = np.copy(FMX)
            else:
                R = R + Rtemp  # at this point it seems like R and E are the sum of all intervals
                E = E + AEtemp
                RABS = RABS + RAB
                FMIN = np.copy(min(np.array([FMIN, FMN]), key=np.abs))
                FMAX = np.copy(max(np.array([FMAX, FMX]), key=np.abs))

        for i in range(NINT + 1, NMAX + 1):  # NMAX is an input parameter
            W[i - 1, 2] = 0

    # print('completed the full integral (DGLQ1)')
    # this is the 15 marker
    if not goto20:
        iterations2 = 0
        while iterations2 < 100:
            iterations2 = iterations2 + 1
            # MAIN SUBPROGRAM LOOP
            if not (100 * EPMACH * RABS >= abs(R) and E < EPS):
                EB = max(100 * UFLOW, max(EPS, 50 * EPMACH) * abs(R))
                if E > EB and BLK0['HighSpeed'] == 0:  # go past here if the function does not pass the error spec
                    if NINT < NMAX - 1:
                        NINT = NINT + 1
                        C = NINT
                    else:  # if the max number of intervals has been reached
                        for i in range(1, NMAX + 1):
                            if W[i - 1, 2] == 0:
                                goto20 = 0
                                IFLAG = 0
                                break

                            goto20 = 1
                            IFLAG = 5

                    LOC = np.argmax(
                        np.abs(W[:, 2]))  # find the max of the errors of all of the intervals in the integral
                    val = max(W[:, 2], key=np.abs)
                    XM = W[LOC, 0] + (
                                W[LOC, 1] - W[LOC, 0]) / 2  # pick a new point dividing the current interval in half
                    if max(abs(W[LOC, 0]), abs(W[LOC, 1])) > (1 + 100 * EPMACH) * (abs(XM) + 1000 * UFLOW):
                        [TR1, TE1, RAB, RAV, FMINL, FMAXL, YROOT] = GL15T(W[LOC, 0], XM, W[LOC, 4], W[LOC, 5], TIMEI,
                                                                          BLK0, BLK1, oldYROOT, itloop, iterations2)
                        oldYROOT = np.copy(YROOT)
                        KF = KF + 15
                        if TE1 < (EB * (XM - W[LOC, 0]) / (B - A)):
                            TE1 = TE1 * SIGN

                        [TR2, TE2, RAB, RAV, FMINR, FMAXR, YROOT] = GL15T(XM, W[LOC, 1], W[LOC, 4], W[LOC, 5], TIMEI,
                                                                          BLK0, BLK1, oldYROOT, itloop, iterations2)
                        oldYROOT = np.copy(YROOT)
                        KF = KF + 15
                        # FMIN = max(FMIN, FMINL, FMINR) # !!! it seems like we never use the FMIN or MAX functions anyway
                        # FMAX = max(FMAX, FMAXL, FMAXR)
                        if TE2 < (EB * (W[LOC, 1] - XM) / (B - A)):
                            TE2 = TE2 * SIGN

                        TE = abs(W[LOC, 2])
                        TR = W[LOC, 3]
                        W[C - 1, :] = [XM, W[LOC, 1], TE2, TR2, W[LOC, 4], W[LOC, 5]]  # this is the new row
                        W[[LOC], [1, 2, 3]] = [XM, TE1, TR1]  # this modifies the the current row to change bounds
                        E = E - TE + abs(TE1) + abs(TE2)
                        R = R - TR + TR1 + TR2
                        if abs(abs(TE1) + abs(
                                TE2) - TE) < .001 * TE:  # if the sum of the errors of the new intervals is less than 1/1000 of the old error
                            IROFF = IROFF + 1
                            if IROFF > 10:
                                IFLAG = 4
                                goto20 = 1
                                break
                    else:
                        if EB > W[LOC, 2]:
                            W[LOC, 2] = 0
                        else:
                            IFLAG = 6
                            goto20 = 1
                            break
                else:  # not E>EB
                    goto20 = 1

            else:  # the thing before
                goto20 = 1

            if goto20 == 1:
                break

        #  the while loop

    # this is the 20
    if IFLAG < 4:
        IFLAG = 3
        T = EPS * abs(R)
        if E < EPS and E < T:
            IFLAG = 0
        elif T < E < EPS:  # E < EPS and E > T
            IFLAG = 1
        elif EPS < E < T:  # E > EPS and E < T
            IFLAG = 2
        else:
            IFLAG = 3

    # debug(W,BLK0,BLK1)

    return R, E, FMIN, FMAX, KF, IFLAG


def DINTG(DINIT, DEND, TIMEI, BLK0, BLK1, itloop):
    EPS = 1e-6
    RST = 0
    NMAX = 50
    NINT = 1
    W = np.array([[DINIT, 0, 0, 0], [DEND, 0, 0, 0]])
    [R, E, FMIN, FMAX, KF, IFLAG] = DGLQ1(DINIT, DEND, EPS, NINT, RST, W, NMAX, TIMEI, BLK0, BLK1, itloop)

    # IFLAG

    return R


def DPHIO(Z, T, C1, C2, A):
    denominator1 = cmath.sqrt(A ** 2 * Z ** 2 - 1)
    denominator2 = cmath.sqrt(Z ** 2 - 1)

    F = -T - C1 * A ** 2 * Z / denominator1 - C2 * Z / denominator2

    return F


def DQR(DL, BLK0, BLK1, oldYROOT, itloop):
    # this is the equivalent of FFUNC, but for surface waves
    code = BLK1['code']

    M = BLK1['M']
    L = BLK1['L']
    K = BLK1['K']

    XD = BLK0['XD']
    DT = BLK0['T']
    alpha = BLK0['alpha']
    # in the future DS and DL and R will be the same thing
    DV = XD / DT

    if code == 1:
        Y = np.real(DV) + cmath.sqrt(-1) * DL
    elif code == 2:
        Y = np.real(DL) + cmath.sqrt(-1) * -0.5
    elif code == 3:
        Y = np.real(1.2 / alpha) + cmath.sqrt(-1) * DL
    else:
        print('Error:code must be 1,2,or 3')

    BLK2 = SFCOM(DL, Y, BLK0, BLK1)
    YQD = QD(BLK0, BLK1, BLK2)
    YUP = QU(BLK2)

    YFEE = np.copy([BLK2['YAA'] / BLK2['YP'] / BLK2['YD'], 1 / BLK2['YQ'] / BLK2['YD'], 1 / BLK2['YQ'] / BLK2['YD']])

    YQR = 0

    if K == 0:
        YQK = [-1, -1, -1]
    elif K == 1:
        YQK = BLK2['YETA1'] / BLK2['Y'] * [1, 1, 1]
    elif K == 2:
        YQK = BLK2['YETA2'] / BLK2['Y'] * [1, 1, 1]
    elif K == 3:
        YQK = [-BLK2['YP'] / BLK2['Y'], -BLK2['YQ'] / BLK2['Y'], -BLK2['YQ'] / BLK2['Y']]
    else:
        print('Error:K must be 0, 1, 2, or 3')

    for i in range(0, 3):
        YQR = YQR + YQD[M - 1, i] * YUP[L - 1, i] * YQK[i] * YFEE[i]

    if code == 1 or code == 3:
        YQR = YQR * cmath.sqrt(-1) / cmath.sqrt(np.real(DV ** 2) - BLK2['YY'])
    else:
        YQR = YQR / cmath.sqrt(np.real(DV ** 2) - BLK2['YY'])

    YINTG = YQR
    YROOT = Y

    return YINTG, YROOT


def DSFING(DINIT, DEND, TIMEI, IPATH, BLK0, BLK1, itloop):
    EPS = 1e-8
    RST = 0
    NMAX = 50
    NINT = 1
    W = np.array([[DINIT, 0, 0, 0, 0, 0], [DEND, 0, 0, 0, 0, 0]])

    if IPATH == 1:
        BLK1['code'] = 1
        [R, E, FMIN, FMAX, KF, IFLAG] = DGLQ1(DINIT, DEND, EPS, NINT, RST, W, NMAX, TIMEI, BLK0, BLK1, itloop)
    #     disp('made it Through IPATH = 1 ____________')
    #     R
    elif IPATH == 2:
        BLK1['code'] = 2
        [R, E, FMIN, FMAX, KF, IFLAG] = DGLQ1(DINIT, DEND, EPS, NINT, RST, W, NMAX, TIMEI, BLK0, BLK1, itloop)
    #     disp('made it Through IPATH = 2 ____________')
    #     R
    elif IPATH == 3:
        BLK1['code'] = 3
        [R, E, FMIN, FMAX, KF, IFLAG] = DGLQ1(DINIT, DEND, EPS, NINT, RST, W, NMAX, TIMEI, BLK0, BLK1, itloop)
    #     disp('made it Through IPATH = 3 ____________')
    #     R
    else:
        print('Error:IPATH is not 1 2 or 3')

    return R


def FFUNC(DS, BLK0, BLK1, oldYROOT, itloop):
    # need to add TIMEI as an input to this function and also where GL15T calls
    # this function

    M = BLK1['M']
    L = BLK1['L']

    # in the future DS and R will be the same thing

    [YROOT, DERR, NROOT] = PHIEQ0(BLK1['AN1'], BLK1['AN2'], BLK0['alpha'], BLK0['T'], BLK0['XD'], DS, oldYROOT, itloop)
    Y = np.copy(YROOT)

    if NROOT == 0:
        YINTG = 0
    else:
        BLK2 = COMMON(DS, Y, BLK0, BLK1)
        YRSTAR = RSTAR(BLK0, BLK1, BLK2, itloop)

        YQD = QD(BLK0, BLK1, BLK2)
        YQS = QSU(BLK0, BLK1, BLK2)
        YINTG = 0
        for i in range(0, 2):
            for j in range(0, 2):
                YINTG = YINTG + YQD[M - 1, i] * YRSTAR[i, j] * YQS[L - 1, j]

        if BLK1['N1'] == 0 and BLK1['N2'] > 0:
            YINTG = YINTG + YQD[M - 1, 2] * YQS[L - 1, 2]

        YDPHI = BLK2['YDPHI']
        YINTG = YINTG * YDPHI / cmath.sqrt(BLK0['XD'] ** 2 - DS ** 2)

    return YINTG, YROOT


def firstOrderFirstDerivative(sig, T):
    vsig = np.zeros(np.size(sig))
    for j in range(2, len(sig) + 1):
        vsig[j - 1] = (sig[j - 1] - sig[j - 2]) / T

    return vsig


def GL15T(A, B, XL, XR, TIMEI, BLK0, BLK1, oldYROOT, itloop, itloop2):
    # computes the integral of G(X) where G(X) = PHI(X)*d(PHI(x))/dx
    #
    # this version was updated on Oct 20, 2009 to incorporate surface waves.
    # This adds the parameter code. O = not a surface wave,
    #                               1 = DQR1
    #                               2 = DQR2
    #                               3 = DQR3

    code = BLK1['code']

    # constants for the Gauss-Kronrod integration
    XGK = [0.9914553711208126,
           0.9491079123427585,
           0.8648644233597691,
           0.7415311855993944,
           0.5860872354676911,
           0.4058451513773972,
           0.2077849550078985,
           0.0]

    WGK = [0.2293532201052922e-1,
           0.6309209262997855e-1,
           0.1047900103222502,
           0.1406532597155259,
           0.1690047266392679,
           0.1903505780647854,
           0.2044329400752989,
           0.2094821410847278]

    WG = [0.1294849661688697,
          0.2797053914892767,
          0.3818300505051189,
          0.4179591836734694]

    FV1 = np.zeros(7, dtype=np.complex128)
    FV2 = np.copy(FV1)

    # PHI(U) = XR - (XR-XL)*U**2*(2*U+3) #!!! I think that this is wrong or maybe (XR-XL)=1               the problem may lie  in the fact that there is a discrepancy between what is PHI(U) and what is supposedly the derivative of this function.
    # PHIP(U) = -6*U*(U+1)

    HLGTH = 0.5 * (B - A)
    CENTR = A + HLGTH
    DHLGTH = abs(HLGTH)

    U = (CENTR - XR) / (XR - XL)
    if code == 0:
        [FMIN, YROOT] = FFUNC(PHI(XL, XR, U), BLK0, BLK1, oldYROOT, itloop)
    else:
        [FMIN, YROOT] = DQR(PHI(XL, XR, U), BLK0, BLK1, oldYROOT, itloop)

    FMAX = np.copy(FMIN)
    FC = FMIN * PHIP(U)
    RESG = FC * WG[3]
    RESK = FC * WGK[7]
    RA = abs(RESK)
    for j in range(1, 3 + 1):  # run through the seven points needed for Gauss Quadrature estimate
        JTW = j * 2
        ABSC = HLGTH * XGK[JTW - 1]
        U = (CENTR - ABSC - XR) / (XR - XL)
        if code == 0:
            [FVAL1, YROOT] = FFUNC(PHI(XL, XR, U), BLK0, BLK1, oldYROOT, itloop)
        else:
            [FVAL1, YROOT] = DQR(PHI(XL, XR, U), BLK0, BLK1, oldYROOT, itloop)

        FMAX = max(np.array([FMAX, FVAL1]), key=np.abs)
        FMIN = min(np.array([FMAX, FVAL1]), key=np.abs)
        FVAL1 = FVAL1 * PHIP(U)

        U = (CENTR + ABSC - XR) / (XR - XL)
        if code == 0:
            [FVAL2, YROOT] = FFUNC(PHI(XL, XR, U), BLK0, BLK1, oldYROOT, itloop)
        else:
            [FVAL2, YROOT] = DQR(PHI(XL, XR, U), BLK0, BLK1, oldYROOT, itloop)

        FMAX = max(np.array([FMAX, FVAL2]), key=np.abs)
        FMIN = min(np.array([FMAX, FVAL2]), key=np.abs)
        FVAL2 = FVAL2 * PHIP(U)
        FV1[JTW - 1] = FVAL1
        FV1[JTW - 1] = FVAL2
        FSUM = FVAL1 + FVAL2
        RESG = RESG + WG[j - 1] * FSUM
        RESK = RESK + WGK[JTW - 1] * FSUM
        RA = RA + WGK[JTW - 1] * (abs(FVAL1) + abs(FVAL2))

    for j in range(1, 4 + 1):  # run through the additional 8 points needed for Kronrod estimate
        JTWM1 = j * 2 - 1
        ABSC = HLGTH * XGK[JTWM1 - 1]
        U = (CENTR - ABSC - XR) / (XR - XL)
        if code == 0:
            [FVAL1, YROOT] = FFUNC(PHI(XL, XR, U), BLK0, BLK1, oldYROOT, itloop)
        else:
            [FVAL1, YROOT] = DQR(PHI(XL, XR, U), BLK0, BLK1, oldYROOT, itloop)

        FMAX = max(np.array([FMAX, FVAL1]), key=np.abs)
        FMIN = min(np.array([FMAX, FVAL1]), key=np.abs)
        FVAL1 = FVAL1 * PHIP(U)

        U = (CENTR + ABSC - XR) / (XR - XL)
        if code == 0:
            [FVAL2, YROOT] = FFUNC(PHI(XL, XR, U), BLK0, BLK1, oldYROOT, itloop)
        else:
            [FVAL2, YROOT] = DQR(PHI(XL, XR, U), BLK0, BLK1, oldYROOT, itloop)

        FMAX = max(np.array([FMAX, FVAL2]), key=np.abs)
        FMIN = min(np.array([FMAX, FVAL2]), key=np.abs)
        FVAL2 = FVAL2 * PHIP(U)
        FV1[JTWM1 - 1] = FVAL1
        FV1[JTWM1 - 1] = FVAL2
        FSUM = FVAL1 + FVAL2
        RESK = RESK + WGK[JTWM1 - 1] * FSUM
        RA = RA + WGK[JTWM1 - 1] * (abs(FVAL1) + abs(FVAL2))

    RESKH = RESK * 0.5
    RASC = WGK[7] * abs(FC - RESKH)
    for j in range(1, 7 + 1):
        RASC = RASC + WGK[j - 1] * abs(FV1[j - 1] - RESKH) + abs(FV2[j - 1] - RESKH)

    # undo normalization
    R = RESK * HLGTH
    RA = RA * DHLGTH
    RASC = RASC * DHLGTH
    AE = abs((RESK - RESG) * HLGTH)
    # !!! didn't include the next part because of EPMACH and UFLOW

    # if RASC~=0 && AE~=0
    #     AE = RASC*min(1,(200*AE/RASC)**1.5)
    # 
    # if RA>UFLOW/(50*EPMACH)
    #     AE = max(EPMACH*50*RA,AE)
    # 

    return R, AE, RA, RASC, FMIN, FMAX, YROOT


def GREENSUB(alpha, INDEX, XD, ZD, NRAY, TIMEI, RAYTIME, KASE, BLK0, BLK1, itloop):
    T = BLK0['T']
    NMAX = NRAY  # !!! not sure if NRAY is correct or not
    DSURF = 0

    if T > RAYTIME[NRAY - 1, 2]:  # if the current time is larger than the arrival of the last ray
        print('Error:not enough rays entered')

    DT = T
    DSURF = 0
    if KASE == 1:  # top surface
        DSURF = SFWAVE(BLK0, BLK1, TIMEI, itloop) / 2 / np.pi

    DISPL = 0

    for N in range(1, NMAX + 1):  # for the number of RAYS
        if T > RAYTIME[N - 1, 2]:
            J2P = 1
            J2S = 2
            AN1 = RAYTIME[N - 1, 0]
            N1 = int(np.real(AN1))  # !!! not sure if this is the proper translation of INT
            if AN1 > N1:  # if the source is buried...
                N1 = N1 + 1
                J2S = 1

            AN2 = RAYTIME[N - 1, 1]
            N2 = int(np.real(AN2))  # !!! not sure if this is the proper translation of INT
            if AN2 > np.real(N2):  # if the source is buried...
                N2 = N2 + 1
                J2P = 2

            if AN2 < 0 and AN2 < N2:
                N2 = N2 - 1
                J2P = 2

            if N1 != 0 or N2 != 0:  # if it is not a Rayleigh wave
                if N1 >= 0 or N2 > 0:  # if it is not a Surface wave
                    if N2 < 0:  # if it is a head wave
                        N2 = -N2
                        AN2 = -AN2
                        DK2 = AN2
                        DTSARR = math.sqrt(XD ** 2 + DK2 ** 2)
                        if T < DTSARR:
                            DHEAD = DT / alpha - DK2 * math.sqrt(1 - alpha ** 2) / alpha
                            BLK1.update({'KASE': KASE})
                            BLK1.update({'code': 0})  # this means that it is not a surface wave
                            BLK1.update({'AN1': AN1})
                            BLK1.update({'AN2': AN2})
                            BLK1.update({'N1': N1})
                            BLK1.update({'N2': N2})
                            DRAYI = DINTG(-XD, -DHEAD, TIMEI, BLK0, BLK1, itloop)
                        else:
                            DRAYI = 0

                    else:  # it is a regular ray
                        BLK1.update({'KASE': KASE})
                        BLK1.update({'code': 0})  # this means that it is not a surface wave
                        BLK1.update({'AN1': AN1})
                        BLK1.update({'AN2': AN2})
                        BLK1.update({'N1': N1})
                        BLK1.update({'N2': N2})
                        DRAYI = DINTG(0, -XD, TIMEI, BLK0, BLK1, itloop)

                    DISPL = DISPL + DRAYI
        else:
            Y = 0

    # done with all the rays
    DISPL = DISPL / 2 / np.pi
    DISPL = DISPL + DSURF  # !!! probably not necessary

    return DISPL


def IBINO(M, N):
    if M == N:
        x = 1
    elif M < 0 or N < 0:
        x = 0
    else:
        K = M - N
        if K < 0:
            x = 0
        else:
            if N == 0:
                x = 1
            else:
                IN = 1
                IP = M
                IQ = 1
                for i in range(1, K + 1):
                    IN = IN * IP / IQ
                    IP = IP - 1
                    IQ = IQ + 1

                x = IN

    return x


def KDELTA(i, j):
    if i == j:
        x = 1
    else:
        x = 0
    return x


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
    #
    #
    # _________________________________________________________________________
    #     the following is a sketch of the subroutine architecture
    # PYTHONPlateSolutionFunc
    #     GREENFCT
    #         INIT
    #             RAYRT
    #         TIMEARRI
    #             RAYRT
    #             TARRV
    #         GREENSUB
    #             SFWAVE (only if surface waves needed)
    #                 DSFING
    #                     DGLQ1
    #                         GL15T - Gauss Kronrod integration
    #                             DQR
    #                                 SFCOM # assigns the common block 2 variables
    #                                 QD
    #                                 QU
    #             DINTG (body waves and head waves)
    #                 DGLQ1
    #                     GL15T - Gauss Kronrod integration
    #                         F
    #                             PHIEQ0 - (Classify and Initial Guess) breaks down the ray into three cases (P only, S only, and P and S) and gives an initial estimate of the root  (YOLD or XST)
    #                                 DCRTNI - implements the minimization by Newtons' method
    #                                     PHIO (function) (done - but not sure about input vars)
    #                                     DPHIO (function) (done - but not sure about input vars)
    #                                 DEQ0
    #                             COMMON # assigns the common block 2 variables
    #                             RSTAR
    #                                 COEFF
    #                                     IBINO
    #                                     KDELTA
    #                                 KDELTA
    #                             QD
    #                             QSU
    # _________________________________________________________________________

    # _________________________________________________________________________
    # some variable definitions
    startTimeProc = time.time()
    alpha = Cs / Cp
    ZD = 0.5 - (zdist / h)  # (0.5 for the top surface, -0.5 for the bottom surface)

    if xdist < 0:
        print('Assuming xdist is positive abs(xdist)')
        xdist = abs(xdist)
    elif xdist == 0:
        print('xdist/=0, setting xdist=1e-16')
        xdist = 1e-12

    XD = xdist / h  # the number of plate thicknesses horizontally away from the source the sensor is
    TDELTA = T / (h / Cs)  # the sampling period normalized to the time take from S wave to traverse the plate thickness
    TMAX = TDELTA * npoint  # total time
    NRAY = 501  # this makes it so there is a maximum of 501 rays
    HighSpeed = 0  # this is a variable that was used occasionally during the debugging process.
    # if HighSpeed = 1, it will skip multiple iterations of the integration
    # subroutine regardless of errors.

    # _________________________________________________________________________
    # define three cases for source sensor location/orientation
    if ZD == 0.5:  # top surface (source on the same side as sensor)
        KASE = 1
    elif ZD == -0.5:  # bottom surface (source opposite sensor)
        KASE = 2
    elif -0.5 < ZD < 0.5:  # source inside the plate somewhere or ZD < 0.5 and ZD > -0.5
        KASE = 3
    else:
        print('Error: Bad source/sensor location: outside the plate')

    # _________________________________________________________________________
    # convert INDEX into three digits --- M is the first digit, then L then K
    if INDEX < 100:
        K = int(0)
        M = int(np.floor(INDEX / 10))
        L = int(INDEX - 10 * M)
    elif 100 < INDEX < 334:  # INDEX > 100 and INDEX < 334
        M = int(np.floor(INDEX / 100))
        L = int(np.floor((INDEX - M * 100) / 10))
        K = int(INDEX - M * 100 - L * 10)
    else:
        print('Error: invalid Greens function index')

    # _________________________________________________________________________
    # create the matrix RAYTIME which contains wave arrival information
    [RCN, TA, CN] = TIMEARRI(alpha, XD, ZD, NRAY)  # this creates a list of all of the rays
    # now sort the result CN
    TAsorted = np.sort(TA)
    TAindexes = np.argsort(TA)
    CN2 = np.zeros(CN.shape, dtype=np.complex128)
    for i in range(len(TA)):
        CN2[i, :] = CN[TAindexes[i], :]

    CN = np.copy(CN2)
    for i in range(1, NRAY + 1):  # Only include those arriving in the timeframe specified
        if CN[i - 1, 2] > TMAX:
            NRAY = i  # find the last ray
            break

    RAYTIME = np.array(CN[0:NRAY, :], dtype=np.complex128)  # cut off all non-necessary rays
    # cut out any repeating rays
    count = 1
    NRAYTIME = np.zeros(RAYTIME.shape, dtype=np.complex128)
    NRAYTIME[0, :] = RAYTIME[0, :]
    for i in range(2, NRAY + 1):
        if RAYTIME[i - 1, 0] != RAYTIME[i - 2, 0] or RAYTIME[i - 1, 1] != RAYTIME[
                i - 2, 1]:  # if this is not a repeat wave
            NRAYTIME[count, :] = RAYTIME[i - 1, :]
            count = count + 1

    NRAY = count
    RAYTIME = NRAYTIME[0:NRAY, :]
    if printOption > 0:
        print('Raytime array\n', [RAYTIME[:, 0:1], RAYTIME[:, 2] * h / Cs])

    # _________________________________________________________________________
    # compute the Green's function
    if printOption < 2:
        # put needed variables into common block BLK0 (a structure is used for
        # MATLAB version
        BLK0 = {}
        BLK1 = {}
        BLK0.update({'alpha': alpha, 'XD': XD, 'T': T, 'HighSpeed': HighSpeed})
        BLK1.update({'M': M, 'L': L, 'K': K})

        # the main part of the program
        DISPL = np.zeros((npoint, 1), dtype=np.complex128)
        for i in range(0, npoint):
            TIMEI = (i + 1) * TDELTA
            BLK0['T'] = (i + 1) * TDELTA
            # this is the main subroutine to compute the Green's function
            currentDisp = GREENSUB(alpha, INDEX, XD, ZD, NRAY, TIMEI, RAYTIME, KASE, BLK0, BLK1, i)
            DISPL[i] = currentDisp
            if i % 1 == 0 and printOption == 1:
                print(i)

        # _________________________________________________________________________
        # dimensionalize the output
        if INDEX < 100:  # if this is the plain Green's function (i.e. Gij)
            DISPL = DISPL / h / Cs ** 2 / rho / np.pi  # re-normalize to plate thickness (m)
        else:  # if it is the first spatial derivative of the Green's function
            DISPL = np.diff(DISPL) / T / h / Cs ** 3 / rho / np.pi

        if printOption == 1:
            totalTime = round(time.time() - startTimeProc)
            print('Total Time = ' + str(totalTime) + ' seconds')

        output = DISPL.real

    else:
        output = 0

    RAYTIME[:, 2] = RAYTIME[:, 2] * h / Cs  # add dimensions to the arrival times

    return output, RAYTIME


def PHI(A, B, U):
    return B - (B - A) * U ** 2 * (2 * U + 3)


def PHIEQ0(DK1, DK2, alpha, T, XD, R, oldYROOT, itloop):
    # this function makes an initial guess of the root of the function PHI,
    # then calls DCRTNI to make sure it is the root or iterate
    #
    # OUTPUTS: YROOT = the root... the Z that satisfies the equation PHIO
    #          DERR = the abs of the function PHIO at this root (should be zero)
    #          NROOT = 1 if its all good, = 0 if convergence was not achieved
    #
    #
    # we are trying to evaluate an integral from 0 to -DXD. 15 points are taken
    # along this interval... DXD = XD = xdist/h. It could be something very
    # very small or something on the order of 10 R is the variable PHI(U) from the GL15T function which is the
    # dependent variable of the integral. U is the point along the integration
    # path. For XD very small (sensor below the source) PHI(U) = some very
    # small negative number.

    # DK1 = AN1, DK2 = AN2, DA = alpha, DT = T, DXD = XD)
    # XD = xdist/h - this in probably only needed for the initial guess
    # here T = TIMEI
    # okay, here are all the inputs: DK1 = # of P trips, DK2 = # of S trips,
    # alpha is the constant cs/cp, T is the current time normalized by h/cs, XD
    # is xdist/h, and R is the variable U from the GL15T function which is the
    # dependent variable of the integral.

    YROOT = oldYROOT
    YOLD = YROOT

    IER = 0
    DYREAL = 0
    YPHI = 0
    A = alpha
    C1 = np.real(DK1)
    C2 = np.real(DK2)
    EPS = 1e-12  # error upper bound
    IEND = 100  # limit on iterations #
    NROOT = 1

    if np.real(DK2) == 0:  # case 1  =  P-Wave only
        DSTSQ = (T / A) ** 2 - C1 ** 2
        DSTAR = cmath.sqrt(DSTSQ)
        if R < 0:  # should be zero or negative for normal waves
            DSTAR = -DSTAR  # therefore DSTAR is turned negative for after the wave arrives.

        if R == 0:
            DYIMAG = -DK1 / (DSTAR * alpha)
            YROOT = 0 + DYIMAG * cmath.sqrt(-1)
        elif R == DSTAR:
            YROOT = T / (alpha ** 2 * DSTAR)
        else:
            DYREAL = R * T / (alpha ** 2 * DSTSQ)
            DSMD = DSTSQ - R ** 2
            # !!! if DSMD < 0 , write into 6
            DYIMAG = -(DK1 / (DSTSQ * alpha)) * cmath.sqrt(DSMD)
            YROOT = DYREAL + cmath.sqrt(-1) * DYIMAG

        YOLD = YROOT
        [YROOT, YPHI, YDPHI, IER] = DCRTNI(YOLD, EPS, IEND, T, R, alpha, C1, C2)

    elif np.real(DK1) == 0:  # case 2: S-Wave only
        DTHARR = alpha * abs(XD) + DK2 * cmath.sqrt(1 - alpha ** 2)
        DTSARR = cmath.sqrt(XD ** 2 + DK2 ** 2)
        DSTSQ = T ** 2 - DK2 ** 2
        DSTAR = cmath.sqrt(DSTSQ)
        if T > np.real(DTSARR):  # ??? if the wave arrival has already passed?
            if R < 0:
                DSTAR = -DSTAR

            if R == 0:
                DYIMAG = -DK2 / abs(DSTAR)
                YROOT = 0 + DYIMAG * cmath.sqrt(-1)
            elif R == DSTAR:
                YROOT = T / DSTAR
            else:
                DYREAL = R * T / DSTSQ
                DYIMAG = -DK2 / DSTSQ * cmath.sqrt(DSTSQ - R ** 2)
                YROOT = DYREAL + cmath.sqrt(-1) * DYIMAG
        elif np.real(DTHARR) < T < np.real(DTSARR):  # T < np.real(DTSARR) and T > np.real(DTHARR)
            DHEAD = T / alpha - DK2 * cmath.sqrt(1 - alpha ** 2) / alpha
            if R < 0:
                DHEAD = -DHEAD
                DSTAR = -DSTAR

            if R == 0:
                DYIMAG = -DK2 / DSTAR
                YROOT = 0 + DYIMAG * cmath.sqrt(-1)
            elif R == DSTAR:
                YROOT = T / DSTAR
            elif abs(R) > abs(DSTAR):
                if R < 0:
                    YROOT = (R * T - DK2 * cmath.sqrt(R ** 2 - DSTSQ)) / DSTSQ
                elif R > 0:
                    YROOT = (R * T + DK2 * cmath.sqrt(R ** 2 - DSTSQ)) / DSTSQ
        else:  # wave has not arrived yet?? should have...
            print('Error:wave has not arrived yet')

        YOLD = YROOT
        [YROOT, YPHI, YDPHI, IER] = DCRTNI(YOLD, EPS, IEND, T, R, alpha, C1, C2)
    elif DK1 != 0 and DK2 != 0:  # case3 C2 not 0 and C1 not 0
        if R == 0 or YOLD == 0:  # if this is the first YROOT to be computed or R = 0 --- need some sort of a guess
            [YROOT, NOROOT] = DEQ0(T, A, C1, C2)  # compute the root using a quartic equation

        YOLD = YROOT  # it seems that YROOT hasn't been defined yet maybe the previous if statement is always true
        [YROOT, YPHI, YDPHI, IER] = DCRTNI(YOLD, EPS, IEND, T, R, alpha, C1, C2)

    else:
        print('Error:should not be here PHIEQ0')

    if IER != 0:  # if there was a problem and the root was not found
        NROOT = 0  # flag NROOT

    DERR = abs(YPHI)  # this is the abs of the function at where the root was (should be zero)

    return YROOT, DERR, NROOT


def PHIO(Z, T, R, C1, C2, A):
    # most of the inputs come from the PHIEQ0 program
    # C1 C2, these are number of P and S trips, A = alpha, Z = the dependent
    # variable of PHIO, which is solved for in DCRTNI. T = the current
    # time normalized by h/cs
    # R is the variable U from GL15T function which is the
    # dependent variable of the integral.
    F = R - T * Z - C1 * cmath.sqrt(A ** 2 * Z ** 2 - 1) - C2 * cmath.sqrt(Z ** 2 - 1)

    return F


def PHIP(U):
    return -6 * U * (U + 1)


def QD(BLK0, BLK1, BLK2):
    M = BLK1['M']
    YP = BLK2['YP']
    YQ = BLK2['YQ']
    YDELTA = BLK2['YDELTA']
    YQSQP1 = BLK2['YQSQP1']
    YQSQM1 = BLK2['YQSQM1']
    YSQRP1 = BLK2['YSQRP1']
    YSQRQ1 = BLK2['YSQRQ1']
    YETA1 = BLK2['YETA1']
    YETA2 = BLK2['YETA2']

    YQD = np.zeros([3, 3], dtype=np.complex128)

    if M == 1:
        YQD[0, 0] = YETA1 * 4 * YP * YQ * YQSQP1 / (YSQRP1 * YDELTA)
        YQD[0, 1] = -2 * YETA1 * YQ * YQSQM1 * YSQRQ1 / YDELTA
        YQD[0, 2] = -2 * YETA2
    elif M == 2:
        YQD[1, 0] = YETA2 * 4 * YP * YQ * YQSQP1 / (YSQRP1 * YDELTA)
        YQD[1, 1] = -2 * YETA2 * YQ * YQSQM1 * YSQRQ1 / YDELTA
        YQD[1, 2] = 2 * YETA1
    elif M == 3:
        YQD[2, 0] = -2 * YP * YQSQM1 * YQSQP1 / (YSQRP1 * YDELTA)
        YQD[2, 1] = -4 * YP * YQ * YSQRQ1 / YDELTA
        YQD[2, 2] = 0

    return YQD


def QSU(BLK0, BLK1, BLK2):
    K = BLK1['K']
    KASE = BLK1['KASE']
    N1 = BLK1['N1']
    N2 = BLK1['N2']

    Y = BLK2['Y']
    YP = BLK2['YP']
    YQ = BLK2['YQ']
    YR0 = BLK2['YR0']
    YRP = BLK2['YRP']
    YRM = BLK2['YRM']
    YSQRP1 = BLK2['YSQRP1']
    YSQRQ1 = BLK2['YSQRQ1']
    YAA = BLK2['YAA']
    YETA1 = BLK2['YETA1']
    YETA2 = BLK2['YETA2']

    YQS = np.zeros([3, 3], dtype=np.complex128)
    YRRP = np.zeros([3, 3], dtype=np.complex128)
    YUP = np.zeros([3, 3], dtype=np.complex128)

    YUP[0, 0] = YETA1 / YSQRP1
    YUP[0, 1] = -YETA1 * YQ / YSQRQ1
    YUP[0, 2] = -YETA2
    YUP[1, 0] = YETA2 / YSQRP1
    YUP[1, 1] = -YETA2 * YQ / YSQRQ1
    YUP[1, 2] = YETA1
    YUP[2, 0] = -YP / YSQRP1
    YUP[2, 1] = -1 / YSQRQ1
    YUP[2, 2] = 0

    YUM = np.copy(YUP)
    YUM[0, 1] = -YUP[0, 1]
    YUM[1, 1] = -YUP[1, 1]
    YUM[2, 0] = -YUP[2, 0]

    YRRP[0, 0] = YR0
    YRRP[0, 1] = YRM
    YRRP[0, 2] = 0
    YRRP[1, 0] = YRP
    YRRP[1, 1] = YR0
    YRRP[1, 2] = 0
    YRRP[2, 0] = 0
    YRRP[2, 1] = 0
    YRRP[2, 2] = -1

    YRRM = np.copy(YRRP)
    YRRM[0, 1] = -YRRP[0, 1]
    YRRM[1, 0] = -YRRP[1, 0]

    if K == 0:
        YP1 = YAA / YP
        YP2 = 1 / YQ
    elif K == 1:
        YP1 = (YETA1 * YAA) / (Y * YP)
        YP2 = YETA1 / (Y * YQ)
    elif K == 2:
        YP1 = (YETA2 * YAA) / (Y * YP)
        YP2 = YETA2 / (Y * YQ)
    elif K == 3:
        YP1 = -YAA / Y
        YP2 = -1 / Y
    else:
        print('Error:SOMETHING IS WRONG, NO if == FITS K OR INDEX')

    if KASE == 1:  # source at the top surface
        for i in range(0, 3):
            YQS[i, 0] = YRRP[0, 0] * YUP[i, 0] * YP1 + YRRP[0, 1] * YUP[i, 1] * YP2 - YUM[i, 0] * YP1
            YQS[i, 1] = YRRP[1, 0] * YUP[i, 0] * YP1 + YRRP[1, 1] * YUP[i, 1] * YP2 - YUM[i, 1] * YP2
            YQS[i, 2] = (YRRP[2, 2] * YUP[i, 2] - YUM[i, 2]) * YP2

    elif KASE == 2:  # source at the bottom surface
        for i in range(0, 3):
            YQS[i, 0] = YRRM[0, 0] * YUM[i, 0] * YP1 + YRRM[0, 1] * YUM[i, 1] * YP2 - YUP[i, 0] * YP1
            YQS[i, 1] = YRRM[1, 0] * YUM[i, 0] * YP1 + YRRM[1, 1] * YUM[i, 1] * YP2 - YUP[i, 1] * YP2
            YQS[i, 2] = (YRRM[2, 2] * YUM[i, 2] - YUP[i, 2]) * YP2

    elif KASE == 3:  # buried source
        if np.remainder(N1 + N2, 2) == 1:  # !!! this assumes that N1 and N2 are integer values
            # N1+N2 ODD, FIRST RAY FROM THE SOURCE IS UP
            for i in range(0, 3):
                YQS[i, 0] = -YUP[i, 0] * YP1
                YQS[i, 1] = -YUP[i, 1] * YP2
                YQS[i, 2] = -YUP[i, 2] * YP2

        else:  # N1+N2 even, FIRST RAY FROM THE SOURCE IS down
            for i in range(0, 3):
                YQS[i, 0] = -YUM[i, 0] * YP1
                YQS[i, 1] = -YUM[i, 1] * YP2
                YQS[i, 2] = -YUM[i, 2] * YP2

    else:
        print('Error:SOMETHING IS WRONG, NO if == FITS KASE')

    return YQS


def QU(BLK2):
    YP = BLK2['YP']
    YQ = BLK2['YQ']
    YSQRP1 = BLK2['YSQRP1']
    YSQRQ1 = BLK2['YSQRQ1']
    YETA1 = BLK2['YETA1']
    YETA2 = BLK2['YETA2']

    YUP = np.zeros([3, 3], dtype=np.complex128)

    YUP[0, 0] = YETA1 / YSQRP1
    YUP[0, 1] = -YETA1 * YQ / YSQRQ1
    YUP[0, 2] = -YETA2
    YUP[1, 0] = YETA2 / YSQRP1
    YUP[1, 1] = -YETA2 * YQ / YSQRQ1
    YUP[1, 2] = YETA1
    YUP[2, 0] = -YP / YSQRP1
    YUP[2, 1] = -1 / YSQRQ1
    YUP[2, 2] = 0

    return YUP


def RAYRT(alpha):
    if alpha <= 0 or alpha >= 1:
        raise Exception('bad value of alpha, the ratio of shear wave speed to long wave speed')

    AA = alpha ** 2
    A = [1, -8, 8 * (3 - 2 * AA), -16 * (1 - AA)]  # these are the coefficients to the Rayleigh equation
    Ans = np.roots(A)  # this is a MATLAB function
    if np.real(Ans[0]) < np.real(Ans[1]) and np.real(Ans[0]) < np.real(Ans[2]):
        YR = Ans[0]
        YI1 = Ans[1]
        YI2 = Ans[2]
    elif np.real(Ans[1]) < np.real(Ans[0]) and np.real(Ans[1]) < np.real(Ans[2]):
        YR = Ans[1]
        YI1 = Ans[0]
        YI2 = Ans[2]
    elif np.real(Ans[2]) < np.real(Ans[0]) and np.real(Ans[2]) < np.real(Ans[1]):
        YR = Ans[2]
        YI1 = Ans[1]
        YI2 = Ans[0]
    else:
        raise Exception('in function RAYRT the if statement')

    return YR, YI1, YI2


def RSTAR(BLK0, BLK1, BLK2, itloop):
    N1 = BLK1['N1']
    N2 = BLK1['N2']
    YR0 = BLK2['YR0']
    YRP = BLK2['YRP']

    YRSTAR = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 1]], dtype=np.complex128)

    for JS in range(1, 2 + 1):
        for J2 in range(1, 2 + 1):
            if J2 == 1:
                N1 = N1 - 1

            if J2 == 2:
                N2 = N2 - 1

            ISIGN = (-1) ** (J2 * (N1 + N2) + N2)

            [KD, KU, KL, KSUM] = COEFF(N1, N2, JS, J2)

            YSUM = 0

            for NM in range(KL, KU + 1):
                YSUM = YSUM + np.real(KSUM[0, NM]) * ((YR0 ** 2 - 1) / YR0 ** 2) ** NM

            KD2 = KDELTA(1, JS) - KDELTA(1, J2)
            YSUM = ((YR0 / YRP) ** KD2) * YSUM
            RSIGN = np.real(ISIGN)
            DRSIGN = RSIGN
            DKD = np.real(KD)
            YRSTAR[JS - 1, J2 - 1] = np.real(DRSIGN) * YR0 ** (N1 + N2) * (np.real(DKD) + YSUM)

            if J2 == 1:
                N1 = N1 + 1

            if J2 == 2:
                N2 = N2 + 1

    return YRSTAR


def SFCOM(DL, Y, BLK0, BLK1):
    AA = (BLK0['alpha'] ** 2)
    YAA = AA
    YA = (BLK0['alpha'])
    DXD = (BLK0['XD'])
    DT = (BLK0['T'])
    DV = (DXD / DT)

    YD = -Y * np.real(DT)
    YETA1 = -Y / np.real(DV)
    YETA2 = cmath.sqrt(1 - YETA1 ** 2)
    YY = Y ** 2
    YP = cmath.sqrt(YAA * YY - 1)
    YQ = cmath.sqrt(YY - 1)
    YPSQP1 = YAA * YY
    YQSQP1 = YY
    YQSQM1 = YQ * YQ - 1
    YPSQM1 = YP * YP - 1
    YSQRP1 = Y * YA
    YSQRQ1 = Y
    YDELTA = YQSQM1 * YQSQM1 + 4 * YP * YQ
    YM = (YQSQM1 * YQSQM1 - 4 * YP * YQ)
    YR0 = YM / YDELTA
    YRP = -(4 * YP * YQSQM1) / (YA * YDELTA)
    YRM = (YR0 * YR0 - 1) / YRP

    BLK2 = {}
    BLK2['YD'] = YD
    BLK2['Y'] = Y
    BLK2['YP'] = YP
    BLK2['YQ'] = YQ
    BLK2['YR0'] = YR0
    BLK2['YRP'] = YRP
    BLK2['YRM'] = YRM
    BLK2['YDELTA'] = YDELTA
    BLK2['YPSQP1'] = YPSQP1
    BLK2['YQSQP1'] = YQSQP1
    BLK2['YPSQM1'] = YPSQM1
    BLK2['YQSQM1'] = YQSQM1
    BLK2['YSQRP1'] = YSQRP1
    BLK2['YSQRQ1'] = YSQRQ1
    BLK2['YY'] = YY
    BLK2['YA'] = YA
    BLK2['YAA'] = YAA
    BLK2['YETA1'] = YETA1
    BLK2['YETA2'] = YETA2
    # BLK2['YDPHI = YDPHI

    return BLK2


def SFWAVE(BLK0, BLK1, TIMEI, itloop):
    N1 = 0
    N2 = 0

    XR = BLK0['XD']
    XD = BLK0['XD']
    alpha = BLK0['alpha']
    DT = BLK0['T']
    TTIME = DT
    DPF = 1.2 / alpha

    BLK1['N1'] = 'hi'  # this is a little test to see if this variable is needed
    BLK1['N2'] = 'hi'

    if TTIME <= alpha * XR:
        DISPL = 0
    else:
        DV = XD / DT
        DISPL = 0
        DISPL = DISPL + DSFING(0, -0.5, TIMEI, 1, BLK0, BLK1, itloop)
        DISPL = DISPL + DSFING(DV, DPF, TIMEI, 2, BLK0, BLK1, itloop)
        DISPL = DISPL + DSFING(-0.5, 0, TIMEI, 3, BLK0, BLK1, itloop)

    return DISPL


def TARRV(alpha, C1, C2, XD, ITER):
    NUM = 0
    DT1 = np.pi / 2  # these are the starting bounds on the ray angle
    DT2 = 0  #

    DA = alpha
    DAA = alpha ** 2
    DN1 = C1
    DN2 = C2
    DR = XD

    if C1 == 0:  # there is no P component to the ray path
        if C2 < 0:  # head wave
            TARR = DR * alpha - DN2 * math.sqrt(1 - alpha ** 2)
            ERR = 0
        else:
            TARR = math.sqrt(DR ** 2 + DN2 ** 2)
            ERR = 0

    elif C2 == 0:  # there is no S component to the ray path
        TARR = math.sqrt(DR ** 2 + DN1 ** 2) * alpha
        ERR = 0
    else:  # both P and S wave paths are present, we must iterate!
        while NUM < ITER:
            NUM = NUM + 1
            DT = (DT1 + DT2) / 2
            DCS2 = math.sqrt(1 - (alpha ** 2 * np.sin(DT) ** 2))
            DRT = C1 * np.tan(DT) + C2 * alpha * np.sin(DT) / DCS2
            # out = [DT DRT C1 C2]
            if abs(DRT - DR) <= 5e-11:
                break

            if DRT <= DR:
                DT2 = DT
            else:
                DT1 = DT

        ITER = NUM
        TARR = alpha * C1 / np.cos(DT) + C2 / DCS2
        ERR = abs(XD - DRT)

    return TARR, ERR


def TIMEARRI(alpha, XD, ZD, NRAY):
    # FORTRAN subroutine TIMEARRI
    # !!! need to find the size of CN first
    CN = np.zeros([1, 3], dtype=np.complex128)

    ITER = 100
    i = 0
    LM = 1
    Z = ZD + 0.5
    Z1 = [0, Z, -Z, 0]
    Z2 = [Z, 0, 0, -Z]
    TA = np.array([])

    if Z == 1:
        # print('sensor is on the top surface of the plate')
        Z2[0] = 0

    OREVEN = Z * 0.5

    for j in range(1, 51 + 1):  # this is the number of C1 rays calculated
        for k in range(1, 22 + 1):  # this is the number of C2 rays calculated
            JJ = j - 1
            KK = k - 1
            if k == 1 and j == 1:  # if it is the first time through the loops
                if Z == 1:  # if the sensor is on the same side as the source (top of plate)
                    # !!! these next few lines could be simplified
                    CN[i, :] = np.array([-1, 0, XD * alpha])
                    TA = np.append(TA, CN[i, 2])
                    i = i + 1
                    CN = np.vstack((CN, [-1, -1, XD]))
                    TA = np.append(TA, CN[i, 2])
                    [YR, YI1, YI2] = RAYRT(alpha)
                    i = i + 1
                    CN = np.vstack((CN, [0, 0, XD / cmath.sqrt(YR)]))
                    TA = np.append(TA, CN[i, 2])

            elif np.remainder(JJ + KK + Z, 2) == 0:
                # !!! go back to the beginning of the k loop
                continue
            else:  # if it is not the first time through the k j loop
                for LL in range(1, 4 + 1):  # !!! increment in steps of LM
                    C1 = JJ + Z1[LL - 1]
                    C2 = KK + Z2[LL - 1]
                    if C1 < 0 or C2 < 0:
                        # !!! goto beginning of LL loop
                        # do nothing
                        if Z == 0 or Z == 1:  # (top or bottom surface of the plate)
                            LM = 4
                            break

                    else:
                        i = i + 1
                        [TARR, ERR] = TARRV(alpha, C1, C2, XD, ITER)
                        if ERR > 1e-6:
                            print('the error in finding the ray takeoff angle by iteration is too great')
                            # error('the error in finding the ray takeoff angle by iteration is too great')

                        if i == 1:
                            TA = np.append(TA, TARR)
                            CN[i - 1, :] = np.array([C1, C2, TARR])
                        else:
                            TA = np.append(TA, TARR)
                            CN = np.vstack((CN, [C1, C2, TARR]))
                        if C1 <= 0:  # if there is no P-wave component
                            if XD > alpha * C2 / math.sqrt(1 - alpha ** 2):
                                C2 = -C2
                                [TARR, ERR] = TARRV(alpha, C1, C2, XD, ITER)
                                i = i + 1
                                CN = np.vstack((CN, [C1, C2, TARR]))
                                TA = np.append(TA, TARR)

                        if Z == 0 or Z == 1:  # (top or bottom surface of the plate)
                            LM = 4
                            break

                    if i > NRAY - 1:
                        # !!! goto out of all the loops
                        break

                #  the LL loop
                if i > NRAY - 1:
                    # !!! goto out of all the loops
                    break

        #  the k loop
        if i > NRAY - 1:
            # !!! goto out of all the loops
            break

    #  the j loop
    RCN = 0  # why was this not assigned before?

    return RCN, TA, CN