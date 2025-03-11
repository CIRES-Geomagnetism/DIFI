# Autogenerated with SMOP 0.29
import datetime
import numpy as np
from typing import Union

from DIFI import jd2000_dt
from DIFI import getmut
from DIFI import gg2gm_2010
from DIFI import design_SHA_Sq_i_Re_v2
from DIFI import design_SHA_Sq_e_Re_v2


def forward_Sq_d_Re(r: Union[float, list], theta: Union[float, list], phi: Union[float, list], t: Union[float, list], f107: Union[float, list], s: dict) -> tuple[np.ndarray, np.ndarray]:
    # [B_1,B_2] = forward_Sq_d_Re(r,theta,phi,t,f107,s)

    # Calculate the primary (1) and induced (2) Sq magnetic field from a model
    # in DIPOLE coordinates (parameterization of Sabaka et al., 2002).
    # The output values B_1 and B_2 are in SPHERICAL Coordinates

    # Inputs:
    #     r(:)                    radius [km]
    #     theta(:), phi(:)        geographic co-latitude, longitude [deg]
    #     t(:)                    time [MJD2000]
    #     f107(:)                 F10.7 index
    #     s is a structure with the following elements
    #         m_e_d_Re(:,:)           primary Sq model (REAL coefficients)
    #         m_i_d_Re(:,:)           secondary Sq model (REAL coefficients)
    #         nmax, mmax              max degree and order in dipole coord.
    #         p_vec(:)                diurnal wavenumbers
    #         s_vec(:)                seasonal wavenumbers
    #         theta_NGP               co-latitude of N. Geomagnetic Pole [deg]
    #         phi_NGP                 longitude of N. Geomagnetic Pole [deg]
    #         h                       altitude of Sq currents [km]
    #         N                       Woldf ratio of the F10.7 dependence

    # Dependendies :
    # - "gg2gm", "jd2000", "getmut"
    # - "design_SHA_Sq_i_Re" and "design_SHA_Sq_e_Re"

    # A. Chulliat, 2016-09-22
    # (from an earlier version dated 2011-04-23, with inputs from N. Olsen)

    # Translated to Python by A. Woods 2018-1-1

    # calculate radii in units of reference radius
    a = 6371.2
    rho = r / a
    rho_Sq = (a + s['h']) / a
    # convert to matrix if input parameter is scalar
    max_size = max(
        np.array(
            (
                [np.size(rho)],
                [np.size(theta)],
                [np.size(phi)],
                [np.size(t)],
            )
        )
    )
    if np.isscalar(rho) or len(rho) == 1:
        rho = np.multiply(rho, np.ones(max_size))

    if np.isscalar(theta) or len(theta) == 1:
        theta = np.multiply(theta, np.ones(max_size))

    if np.isscalar(phi) or len(phi) == 1:
        phi = np.multiply(phi, np.ones(max_size))

    if np.isscalar(t) or len(t) == 1:
        t = np.multiply(t, np.ones(max_size))

    if np.isscalar(f107) or len(f107) == 1:
        f107 = np.multiply(f107, np.ones(max_size))

    if (
        np.size(t) != np.size(rho)
        or np.size(t) != np.size(theta)
        or np.size(t) != np.size(phi)
        or np.size(t) != np.size(f107)
    ):
        print(
            np.size(t),"np.size(t)",
            np.size(rho),"np.size(rho)",
            np.size(theta),"np.size(theta)",
            np.size(phi),"np.size(phi)",
            np.size(f107),"np.size(f107)"
        )
        raise Exception(r"Variables must be of equal size (or scalars)")

    w_s = 2*np.pi
    w_p = 2*np.pi / 24
    N_data = np.size(theta, 0)

    # calculate time in year (season) and MUT
    t_1 = [
        datetime.timedelta(int(t_iter), (t_iter - int(t_iter)) * 3600)
        + datetime.datetime(2000, 1, 1)
        for t_iter in t
    ]
    year = [t_iter.year for t_iter in t_1]
    ndays = np.zeros(np.size(year))
    for i in range(len(year)):
        ndays[i] = (
            datetime.datetime(year[i]+1, 1, 1)
            - datetime.datetime(year[i], 1, 1)
        ).days
    # list comprehension must be converted to ndarray
    t_season = np.array(
        [
            (t[i] - jd2000_dt.jd2000_dt(year[i], 1, 1, 0)) / ndays[i]
            for i in range(len(ndays))
        ]
    )
    # not inlcuded in HDGM version
    t_season = t_season.flatten()
    t_mut = getmut.getmut(t, s['theta_NGP'], s['phi_NGP'])

    # calculate dipolar coordinates + matrix R
    theta_d, phi_d, rotmat = gg2gm_2010.gg2gm_2010(theta, phi, get_R=True)

    arr_internal = np.array(
        design_SHA_Sq_i_Re_v2.design_SHA_Sq_i_Re_v2(
            rho,
            theta_d,
            phi_d,
            s['nmax'],
            s['mmax'],
        )
    )

    s_vec = np.array(s['s_vec'])
    p_vec = np.array(s['p_vec'])

    beta = (
        w_s * s_vec[:, None, None] * t_season[None, None, :]
        + w_p * p_vec[None, :, None] * t_mut[None, None, :]
    )
    beta = beta.reshape(
        -1, N_data
    )
    time_arr = np.array(
        [
            np.cos(beta),
            np.sin(beta),
        ]
    )
    time_arr = time_arr.transpose(1, 0, 2).reshape(-1, N_data)

    # CASE #1: above Sq currents

    if (min(rho) > rho_Sq):
        B_1_tmp = np.einsum(
            'ij, ikl, jk ->lk',
            s['m_e_d_Re'],
            arr_internal,
            time_arr,
        )

        B_2_tmp = np.einsum(
            'ij, ikl, jk ->lk',
            s['m_i_d_Re'],
            arr_internal,
            time_arr,
        )
    # CASE #2: below Sq currents
    elif (max(rho) < rho_Sq):
        arr_external = np.array(
            design_SHA_Sq_e_Re_v2.design_SHA_Sq_e_Re_v2(
                rho,
                theta_d,
                phi_d,
                s['nmax'],
                s['mmax'],
            )
        )

        B_1_tmp = np.einsum(
            'ij, lik, jk ->lk',
            s['m_e_d_Re'],
            arr_external,
            time_arr,
        )

        B_2_tmp = np.einsum(
            'ij, lik, jk ->lk',
            s['m_i_d_Re'],
            arr_internal,
            time_arr,
        )

    # CASE #3: error
    else:
        raise Exception(r"data in both regions (below and above Sq currents)")

    # Rotate into geomagnetic frame
    # fix minus sign
    B_1_tmp[[1, 2]] = np.einsum(
        'kij, jk -> ik',
        rotmat.transpose(0, 2, 1),
        B_1_tmp[[1, 2]],
    )

    B_2_tmp[[1, 2]] = np.einsum(
        'kij, jk -> ik',
        rotmat.transpose(0, 2, 1),
        B_2_tmp[[1, 2]],
    )

    # correct for F10.7 dependence
    w = (1 + s['N']*f107)
    B_1 = B_1_tmp * w
    B_2 = B_2_tmp * w

    return B_1, B_2
