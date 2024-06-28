import numpy as np


def legendre(lat, nMax):

    """
    Computes  all of the Schmidt-semi normalized associated Legendre
    functions up to degree nMax. If nMax <= 16, function MAG_PcupLow is used.
    Otherwise MAG_PcupHigh is called.
    INPUT  Latitude (degrees)
    nMax    integer    (Maximum degree of spherical harmonic secular model)
        LegendreFunction Pointer to data structure with the following elements
    double *Pcup  (  pointer to store Legendre Function  )
    double *dPcup ( pointer to store  Derivative of Legendre function )

    LegendreFunction  Calculated Legendre variables in the data structure
    """

    sin_phi = np.sin(np.radians(lat))

    #    if nMax <= 16 or (1-abs(sin_phi)) < 10**-10:
    [P, dP] = MAG_PcupLow(sin_phi, nMax)
    # hoping that overflow doesn't occur
    #    else:
    #        raise ValueError("Implement Pcup High for this", nMax)

    return [P, dP]


def MAG_PcupLow(x, nMax):
    """ This function evaluates all of the Schmidt-semi normalized associated Legendre
        functions up to degree nMax.

        Calling Parameters:
                INPUT
                        nMax:	 Maximum spherical harmonic degree to compute.
                        x:		cos(colatitude) or sin(latitude).

                OUTPUT
                        Pcup:	A vector of all associated Legendgre polynomials evaluated at
                                        x up to nMax.
                   dPcup: Derivative of Pcup(x) with respect to latitude

        Notes: Overflow may occur if nMax > 20 , especially for high-latitudes.
        Use MAG_PcupHigh for large nMax.

    Written by Manoj Nair, June, 2009 . Manoj.C.Nair@Noaa.Gov.

    Note: In geomagnetism, the derivatives of ALF are usually found with
    respect to the colatitudes. Here the derivatives are found with respect
    to the latitude. The difference is a sign reversal for the derivative of
    the Associated Legendre Functions.
    """

    NumTerms = int((nMax + 1) * (nMax + 2) / 2)
    try:
        Pcup = np.zeros((NumTerms, x.size))
        dPcup = np.zeros((NumTerms, x.size))
    except AttributeError:
        Pcup = np.zeroes((NumTerms, 1))
        dPcup = np.zeroes((NumTerms, 1))
    Pcup[0] = 1.0
    dPcup[0] = 0.0
    # sin (geocentric latitude) - sin_phi
    z = np.sqrt((1.0 - x) * (1.0 + x))
    schmidtQuasiNorm = np.zeros((NumTerms, x.size))

    # schmidtQuasiNorm = (double *) malloc((NumTerms + 1) * sizeof ( double))

    # First, Compute the Gauss-normalized associated Legendre functions
    for n in range(1, nMax+1):
        for m in range(n+1):
            index = int(n * (n + 1) / 2 + m)
            if n == m:
                index1 = int((n - 1) * n / 2 + m - 1)
                Pcup[index] = z * Pcup[index1]
                dPcup[index] = z * dPcup[index1] + x * Pcup[index1]
            elif n == 1 and m == 0:
                index1 = int((n - 1) * n / 2 + m)
                Pcup[index] = x * Pcup[index1]
                dPcup[index] = x * dPcup[index1] - z * Pcup[index1]
            elif n > 1 and n != m:
                index1 = int((n - 2) * (n - 1) / 2 + m)
                index2 = int((n - 1) * n / 2 + m)
                if m > n - 2:
                    Pcup[index] = x * Pcup[index2]
                    dPcup[index] = x * dPcup[index2] - z * Pcup[index2]
                else:
                    k = (
                        np.float64(((n - 1) * (n - 1)) - (m * m))
                        / np.float64((2 * n - 1) * (2 * n - 3))
                    )
                    Pcup[index] = x * Pcup[index2] - k * Pcup[index1]
                    dPcup[index] = (
                        x * dPcup[index2]
                        - z * Pcup[index2]
                        - k * dPcup[index1]
                    )

    # Compute the ratio between the the Schmidt quasi-normalized associated
    # Legendre functions and the Gauss-normalized version.

    schmidtQuasiNorm[0] = 1.0
    for n in range(1, nMax+1):
        index = int(n * (n + 1) / 2)
        index1 = int((n - 1) * n / 2)
        # for m = 0
        schmidtQuasiNorm[index] = (
            schmidtQuasiNorm[index1]
            * np.float64((2 * n - 1))
            / np.float64(n)
        )
        for m in range(1, n+1):
            index = int((n * (n + 1) / 2 + m))
            index1 = int((n * (n + 1) / 2 + m - 1))
            schmidtQuasiNorm[index] = (
                schmidtQuasiNorm[index1]
                * np.sqrt(
                    np.float64((n - m + 1) * (2 if m == 1 else 1))
                    / np.float64(n + m)
                )
            )

    # Converts the  Gauss-normalized associated Legendre
    # functions to the Schmidt quasi-normalized version using pre-computed
    # relation stored in the variable schmidtQuasiNorm

    for n in range(1, nMax+1):
        for m in range(n+1):
            index = int(n * (n + 1) / 2 + m)
            Pcup[index] = Pcup[index] * schmidtQuasiNorm[index]
            dPcup[index] = -dPcup[index] * schmidtQuasiNorm[index]
            # The sign is changed since the new WMM routines use derivative
            # with respect to latitude instead of co-latitude

    return [Pcup, dPcup]
    # MAG Pcup low
