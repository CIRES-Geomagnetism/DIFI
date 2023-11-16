# Autogenerated with SMOP 0.29
import numpy as np
#from scipy import special
import math
import legendre
# 

    
def design_SHA_Sq_i_Re_v2(rho,theta,phi,t,t_ut,nmax,mmax,p_vec,s_vec):

    # [A_r, A_theta, A_phi] = design_SHA_Sq_i_Re_v2(rho,theta,phi,t,t_ut,...
#    nmax,mmax,p_vec,s_vec);
    
    # Calculate design matrices A_i that connects the vector of REAL,
# Schmidt-normalized, spherical harmonic expansion coefficients,
#        x = g_{nsp}^{m} and h_{nsp}^{m}
# and the magnetic component B_r, B_theta, B_phi:
#        B_r     = x'*A_r
#        B_theta = x'*A_theta
#        B_phi   = x'*A_phi
# for the region above (primary/secondary) Sq currents.
    
    # Inputs:   rho(:)                  radius [units of reference radius]
#           theta(:), phi(:)        co-latitude, longitude [deg]
#           t(:), t_ut(:)           time [yr]; UT [h]
#           nmax, mmax              maximum degree and order
#           p_vec(:)                diurnal wavenumbers
#           s_vec(:)                seasonal wavenumbers
    
    # (Optimized version)
    
    # A. Chulliat, 2016-09-22
# (from an earlier version dated 2011-04-21, with inputs from N. Olsen)
    
    rad=np.pi / 180
    w_s=2*np.pi
    w_p=2*np.pi / 24
    N_data=np.size(theta,0)
    if (np.size(phi,0) != N_data):
        raise ValueError(r"design_SHA_Sq_i: theta and phi have different dimensions")
    if np.isscalar(theta) is not True and theta.ndim > 1:
        dim2=np.size(theta,1)
    else:
        dim2 = 1
    if dim2 > 1:
        theta=theta.transpose()
        phi=phi.transpose()
    
    cos_theta=np.cos(theta*rad)
    sin_theta=np.sin(theta*rad)
    # convert to row vector if input parameter is scalar
    if np.isscalar(rho) is True:
        rho=rho*np.ones(N_data)
    
    if np.isscalar(t) is True:
        t=t*np.ones(N_data)

    if np.isscalar(t_ut) is True:
        t_ut=t_ut*np.ones(N_data)
    
    # number of parameters
    
    N_nm=mmax*(mmax + 2) + (nmax - mmax)*(2*mmax + 1)
    N_sp=np.size(p_vec)*np.size(s_vec)
    N_coeff=2*N_nm*N_sp
    # calculate sub-matrices for s=0, p=0
    #This looks strange, but np.zeros takes a tuple as input, thus the extra paren is necessary
    A_r_0=np.zeros((N_nm,N_data)) 
    A_theta_0=np.zeros((N_nm,N_data))
    A_phi_0=np.zeros((N_nm,N_data))
    k=0 #Cycles through all the coefficients for degree n, order m model.
    Pnm = np.zeros((mmax+2,nmax+1,np.size(cos_theta,0)))
    dPnm = np.zeros((mmax+2,nmax+1,np.size(cos_theta,0)))

    [Pnm, dPnm] = legendre.legendre(90-theta, nmax)


    #Cycle through all of n and m
    for n in xrange(1, nmax+1):
        rn1=rho**(-(n+2))
        for m in xrange(min(n,mmax)+1): #The maximum m is less than the maximum n in the default DIFI
            index = n * (n + 1) / 2 + m
            if m==0:
                #no h terms for g10, g20, g30 etc...
                A_r_0[k] = (n+1)*rn1*Pnm[index]
                A_theta_0[k] = rn1*dPnm[index]
                A_phi_0[k] = -rn1*0
                k+=1
            else:
                alpha = m*phi*rad
                cos_phi = np.cos(alpha)
                sin_phi = np.sin(alpha)
                #g terms, as in g11, g21, g22 etc...
                A_r_0[k] = (n+1)*rn1*Pnm[index]*cos_phi
                A_theta_0[k] = rn1*dPnm[index]*cos_phi
                A_phi_0[k] = -rn1*Pnm[index]/sin_theta*(-m*sin_phi)
                k+=1
                #h terms as in h11, h21, h22 etc...
                A_r_0[k] = (n+1)*rn1*Pnm[index]*sin_phi
                A_theta_0[k] = rn1*dPnm[index]*sin_phi
                A_phi_0[k] = -rn1*Pnm[index]/sin_theta*m*cos_phi;
                k+=1


    for s in np.squeeze(s_vec): #s_vec and p_vec are two dimensional, the first dimension is 1
        for p in np.squeeze(p_vec):
            if s==0 and p==0:
                A_r = np.vstack((A_r, A_r_0))
                A_theta = np.vstack((A_theta, A_theta_0))
                A_phi = np.vstack((A_phi, A_phi_0))

                A_r = np.vstack((A_r,0*A_r_0))
                A_theta = np.vstack((A_theta, 0*A_theta_0))
                A_phi = np.vstack((A_phi,0*A_phi_0))
            else:
                beta = w_s*s*t + w_p*p*t_ut
                cos_beta = np.cos(beta)
                sin_beta = np.sin(beta)
                try:
                    A_r = np.vstack((A_r, A_r_0*cos_beta))
                    A_theta = np.vstack((A_theta, A_theta_0*cos_beta))
                    A_phi = np.vstack((A_phi, A_phi_0*cos_beta))
                except UnboundLocalError:
                    A_r = A_r_0*cos_beta
                    A_theta = A_theta_0*cos_beta
                    A_phi = A_phi_0*cos_beta
                    

                A_r = np.vstack((A_r, A_r_0*sin_beta))
                A_theta = np.vstack((A_theta, A_theta_0*sin_beta))
                A_phi = np.vstack((A_phi, A_phi_0*sin_beta))

    return A_r, A_theta, A_phi
