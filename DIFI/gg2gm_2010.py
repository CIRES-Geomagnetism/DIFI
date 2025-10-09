import math
import numpy as np



def gg2gm_2010(
    theta_gg: np.ndarray=None,
    phi_gg: np.ndarray=None,
    theta_gm=None,
    phi_gm=None,
    B_theta=None,
    B_phi=None,
    get_R: bool=False, 
    theta_b = 9.92,
    phi_b = 287.78
):
    # transformation between geographic (geocentric as opposed to geodetic)
    # and geomagnetic (dipole) co-ordinates
    #          and components

    # theta_gm                = gg2gm(theta_gg, phi_gg)
    # [theta_gm, phi_gm]      = gg2gm(theta_gg, phi_gg)
    # [theta_out, phi_out]    = gg2gm(theta_in, phi_in, i_trans)
    #                         i_trans = +1: gg -> gm
    #                         i_trans = -1: gm -> gg

    # theta(:) co-latitude [deg]
    # phi(:)   longitude [deg]

    # R(:,2,2) is the matrix that rotates a horizontal vector from GEO to MAG:
    #          [B_theta_gm; B_phi_gm] = R*[B_theta_geo; B_phi_geo]

    # March 2003, Nils Olsen, DSRI

    # September 2004: Output of matrix R added

    # February 2011: coordinates of geomagnetic North pole updated to 2010.0
    # (Arnaud Chulliat, IPGP)
    # December 2016: translated to python
    # (Adam Woods, NCEI)
    # October 2025: changed pole location to be argument | Default IGRF 2010 pole
    # (Collin Kadlecek, NCEI)

    rad = math.pi / 180

    if (
        (theta_gg is not None and phi_gg is not None)
        and (theta_gm is None or phi_gm is None)
    ):
        i_trans = 1
        theta_in = theta_gg
        phi_in = phi_gg
    elif (
        (theta_gg is None and phi_gg is None)
        and (theta_gm is not None or phi_gm is not None)
    ):
        i_trans = -1
        theta_in = theta_gm
        phi_in = phi_gm
    else:
        print(
            "Need either values for phi_gg and theta_gg xor phi_gm and "
            "theta_gm, not neither or both"
        )

    s_p_b = np.sin(phi_b*rad)
    c_p_b = np.cos(phi_b*rad)
    c_t_b = np.cos(theta_b*rad)
    s_t_b = np.sin(theta_b*rad)
    A = np.array(
        [
            [c_t_b*c_p_b, c_t_b*s_p_b,  -s_t_b],
            [-s_p_b, c_p_b, 0],
            [s_t_b*c_p_b, s_t_b*s_p_b, c_t_b],
        ]
    )

    if i_trans == - 1:
        A = A.transpose()

    c_t = np.cos(theta_in*rad)
    s_t = np.sin(theta_in*rad)
    c_p = np.cos(phi_in*rad)
    s_p = np.sin(phi_in*rad)
    z = c_t
    x = s_t*c_p
    y = s_t*s_p
    x_out = A[0, 0] * x + A[0, 1] * y + A[0, 2] * z
    y_out = A[1, 0] * x + A[1, 1] * y + A[1, 2] * z
    z_out = A[2, 0] * x + A[2, 1] * y + A[2, 2] * z
    theta_out = 90 - np.arctan2(z_out, np.sqrt(x_out ** 2 + y_out ** 2)) / rad
    phi_out = (np.arctan2(y_out, x_out) / rad) % 360

    # Still need to pythonify the below...
    if B_theta is not None:
        raise NotImplementedError(
            "The rotation of the Magnetic field vector is not implemented"
        )
        BE = B_theta*c_t
        Bx = BE*c_p - B_phi*s_p
        By = BE*s_p + B_phi*c_p
        Bz = -B_theta*s_t
        Bx_gm = A[0, 0] * Bx + A[0, 1] * By + A[0, 2] * Bz
        By_gm = A[1, 0] * Bx + A[1, 1] * By + A[1, 2] * Bz
        Bz_gm = A[2, 0] * Bx + A[2, 1] * By + A[2, 2] * Bz
        c_t = np.cos(theta_gm*rad)
        s_t = np.sin(theta_gm*rad)
        c_p = np.cos(phi_gm*rad)
        s_p = np.sin(phi_gm*rad)
        BE = Bx_gm*c_p + By_gm*s_p
        varargout[1] = np.multiply(BE, c_t) - np.multiply(Bz_gm, s_t)
        varargout[2] = np.multiply(By_gm, c_p) - np.multiply(Bx_gm, s_p)

    # Still need to pythonify the above
    if get_R:
        if i_trans != 1:
            raise ValueError(r"Calculation of R only for i_trans = 1")
        theta_gm = theta_out
        phi_gm = phi_out
        c_t = np.cos(theta_gg*rad)
        s_t = np.sin(theta_gg*rad)
        c_p = np.cos(phi_gg*rad)
        s_p = np.sin(phi_gg*rad)
        R_1 = np.zeros((np.size(c_p, 0), 3, 3))
        R_1[:, 0, :] = np.array([s_t*c_p, c_t*c_p, -s_p]).transpose()
        R_1[:, 1, :] = np.array([s_t*s_p, c_t*s_p, c_p]).transpose()
        R_1[:, 2, :] = np.array(
            [c_t, -s_t, np.zeros(np.size(c_p))]
        ).transpose()
        c_t = np.cos(theta_gm*rad)
        s_t = np.sin(theta_gm*rad)
        c_p = np.cos(phi_gm*rad)
        s_p = np.sin(phi_gm*rad)
        R_2 = np.zeros((np.size(c_p, 0), 3, 3))
        R_2[:, 0, :] = np.array([s_t*c_p, s_t*s_p, c_t]).transpose()
        R_2[:, 1, :] = np.array([c_t*c_p, c_t*s_p, -s_t]).transpose()
        R_2[:, 2, :] = np.array(
            [-s_p, c_p, np.zeros(np.size(c_p))]
        ).transpose()
        s_p_b = s_p_b*np.ones(np.size(c_p, 0))
        c_p_b = c_p_b*np.ones(np.size(c_p, 0))
        s_t_b = s_t_b*np.ones(np.size(c_p, 0))
        c_t_b = c_t_b*np.ones(np.size(c_p, 0))
        # While they look like extra parentheticals these are necessary as
        # zeros takes a tuple as a parameter
        R_mag_geo = np.zeros((np.size(c_p, 0), 3, 3))
        R_mag_geo[:, 0, :] = np.array(
            [c_t_b*c_p_b, c_t_b*s_p_b, -s_t_b]
            ).transpose()
        R_mag_geo[:, 1, :] = np.array(
            [-s_p_b, c_p_b, np.zeros(np.size(c_p, 0))]
            ).transpose()
        R_mag_geo[:, 2, :] = np.array(
            [s_t_b*c_p_b, s_t_b*s_p_b, c_t_b]
            ).transpose()
        R_tmp = np.einsum('ijk,ikl->ijl', R_mag_geo, R_1)
        R_3 = np.einsum('ijk,ikl->ijl', R_2, R_tmp)
        sheet_1, sheet_2, sheet_3 = np.dsplit(R_3, 3)
        sm_R3 = np.dstack((sheet_2, sheet_3))
        column_1, column_2, column_3 = np.hsplit(sm_R3, 3)
        sm_R3 = np.hstack((column_2, column_3))
        # peeling off top sheet and column
        return [theta_out, phi_out, sm_R3]
    else:
        return [theta_out, phi_out]
