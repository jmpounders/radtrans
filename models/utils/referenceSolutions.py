#!/usr/bin/python

import numpy as np
import quadrature as quad
import aux

import matplotlib.pyplot as plt

def oneGroupHomogSlabDiffusion( x, w, D, sigma, q ) :
    """
    Return the solution of
        -D phi'' + sigma * phi = q
    in homogeneous slab geometry where 
       sigma = constant absorption cross section
       D = constant diffusion coefficient
       w = slab width
       q = constant external souce
       x = possition array
    """
    L = np.sqrt( D/sigma )
    A = -q/sigma/(np.exp(w/(2*L)) + np.exp(-w/(2*L)))
    return A*(np.exp(x/L) + np.exp(-x/L)) + q/sigma


def twoGroupHomogSlabDiffusionNoUpScatter( x, w, D, sigma_a, sigma_s, q ) :
    """
    Return the solution of
        -D_1 phi_1'' + sigma_a_1 * phi_1 = q_1
        -D_2 phi_2'' + sigma_a_2 * phi_2 = sigma_s_12 * phi_1 + q_2
    in homogeneous slab geometry where 
       sigma_a = constant absorption cross section array
       sigma_s = constant scattering matrix
       D = constant diffusion coefficient array
       w = slab width
       q = constant external souce array
       x = possition array
    """
    L1 = np.sqrt( D[0]/sigma_a[0] )
    L2 = np.sqrt( D[1]/sigma_a[1] )
    C1 = q[0]/sigma_a[0]
    C2 = sigma_s[1,0]/sigma_a[1] * C1 + q[1]/sigma_a[1]
    A1 = -q[0]/sigma_a[0]/(np.exp(w/(2*L1)) + np.exp(-w/(2*L1)))
    B2 = sigma_s[1,0]*A1 / ( -D[1]/(L1**2) + sigma_a[1] )
    A2 = (-C2 - B2*(np.exp(w/(2*L1)) + np.exp(-w/(2*L1))))/(np.exp(w/(2*L2)) + np.exp(-w/(2*L2)))

    flux1 = A1*(np.exp(x/L1) + np.exp(-x/L1)) + C1
    flux2 = A2*(np.exp(x/L2) + np.exp(-x/L2)) + B2*(np.exp(x/L1) + np.exp(-x/L1)) + C2
    return np.column_stack((flux1, flux2))


def multiGroupHomogSlabDiffusion( x, w, D, sigma_a, sigma_s, nu_sigma_f, chi, q ) :
    """
    Return the solution of the multigroup homogeneous diffusion equation with fission and scattering
    present where
       x = possition array
       w = slab width
       D = constant diffusion coefficient array
       sigma_a = constant absorption cross section array
       sigma_s = constant scattering matrix
       nu_sigma_f = constant production cross section
       chi = fission spectrum
       q = constant external souce array
    """
    numGroups = D.shape[0]

    for g in range(numGroups) :
        sigma_s[g,g] = 0

    # Get the eigenvalues and eigenfunctions for the homogenous part of the solution
    F = np.outer(chi, nu_sigma_f)
    AmSmF = np.diag(sigma_a) - sigma_s - F
    M = np.dot( np.diag(1./D), AmSmF )
    d, V = np.linalg.eig( M )

    # Get the constant (particular) solution
    L = 1/np.sqrt(d)
    C = np.linalg.solve(AmSmF, q)

    # Solve for the eigenfunction expansion coefficients
    W = np.zeros( V.shape )
    for i in range(numGroups) :
        for j in range(numGroups) :
            W[i,j] = V[i,j] * (np.exp( w/(2*L[j]) ) + np.exp( -w/(2*L[j]) ))
    a = np.linalg.solve( W, -C )

    # Form solution
    flux =  np.zeros( (numGroups, x.shape[0]) )
    for i in range(x.shape[0]) :
        flux[:,i] = C

    for k in range(numGroups) :
        for i in range(numGroups) :
            flux[k,:] = flux[k,:] + V[k,i]*a[i]*(np.exp(x/L[i]) + np.exp(-x/L[i]))

    return flux.T

def oneGroupEigenvalueSlabDiffusion( x, w, D, sigma_a, nu_sigma_f, norm=1.0 ) :
    Bn = np.pi/w
    k = nu_sigma_f/(D*Bn**2 + sigma_a)
    flux = (Bn*norm/2.0)*np.cos( Bn * x )
    return (flux, k)

def oneGroupUniformTransientSlabDiffusion( x, t, w, D, sigma_a, nu_sigma_f, beta, lmbda, v, kcrit, norm=1.0 ) :
    Bn = np.pi/w
    nu_sigma_f = nu_sigma_f/kcrit

    alpha = -D*Bn**2 - sigma_a + (1-beta)*nu_sigma_f
    k1 = (v*alpha - lmbda + np.sqrt((lmbda+v*alpha)**2 + 4*lmbda*v*beta*nu_sigma_f))/2.0
    k2 = (v*alpha - lmbda - np.sqrt((lmbda+v*alpha)**2 + 4*lmbda*v*beta*nu_sigma_f))/2.0

    phi0 = Bn*norm/2.0
    A1 = -phi0*k2/(lmbda*(k2+lmbda)) / ( k1/(lmbda*(k1+lmbda)) - k2/(lmbda*(k2+lmbda)) )
    A2 = phi0 - A1

    flux0 = np.cos( Bn * x )
    T = A1*np.exp(k1*t) + A2*np.exp(k2*t)
    flux = np.outer(flux0,T)
    P = 2/Bn*T

    return (flux, P)

def oneGroupHomogAbsorbingSlabTransport(p, slabWidth, q, sigma, N) :
    """
    Returns the scalar flux for a 2D square, purely absorbing transport problem with
    a discrete quadrature set in angle.
    Input:
      p = (x,y) location
      slabWidth = the width of the problem in each direction
      q = isotropic uniform source strength (per solid angle)
      sigma = total macroscopic cross section
      N = quadrature order
    """

    (omega_x,omega_y,omega_z,wgt) = quad.getQuadrature(N)
    
    bdries = ((slabWidth, slabWidth),
              (      0.0, slabWidth),
              (      0.0,       0.0),
              (slabWidth,       0.0),
              (slabWidth, slabWidth))

    scalarFlux = 0.0
    angularFlux = np.zeros(len(wgt))
    for b in range(4) :
        b0 = bdries[b]
        b1 = bdries[b+1]
        phi0 = aux.getAngleFromVector( b0[0]-p[0], b0[1]-p[1] )
        phi1 = aux.getAngleFromVector( b1[0]-p[0], b1[1]-p[1] )
        if phi0 > phi1 : phi1 = phi1 + 2.0*np.pi
        d01 = np.sqrt( (p[0]-b0[0])**2 + (p[1]-b0[1])**2 )
        d12 = np.sqrt( (b1[0]-b0[0])**2 + (b1[1]-b0[1])**2 )
        d20 = np.sqrt( (p[0]-b1[0])**2 + (p[1]-b1[1])**2 )
        theta0 = phi1 - phi0
        theta1 = np.arcsin(min(1.0,d20*np.sin(theta0)/d12))
        if b==4 :
            print theta0*180/np.pi, theta1*180/np.pi, d12, d20
        if d01 < 1e-8 or d20 < 1e-8:
            print p
            continue
        #plt.plot( [p[0],b0[0],b1[0],p[0]], [p[1],b0[1],b1[1],p[1]],'-o' )
        
        for n in range(len(wgt)) :
            phi = aux.getAngleFromVector( omega_x[n], omega_y[n] )
            phi = np.mod( phi+np.pi, 2.0*np.pi )
            if (phi>=phi0 and phi<phi1) or (phi1>=phi0-2.0*np.pi and phi<phi1-2.0*np.pi) :
                theta = phi - phi0
                S = d01*np.sin(theta1)/np.sin(theta1+theta)/np.sqrt(1.0-omega_z[n]**2)
                #plt.plot( [p[0],p[0]+S*np.cos(phi)*np.sqrt(1.0-omega_z[n]**2)], [p[1],p[1]+S*np.sin(phi)*np.sqrt(1.0-omega_z[n]**2)])
                angularFlux[n] = (1.0 - np.exp(-sigma*S))*q/sigma
                scalarFlux += wgt[n]*angularFlux[n]

    #plt.show()
    #return np.concatenate([np.array([scalarFlux/(4.0*np.pi)]), angularFlux])
    return scalarFlux/(4.0*np.pi)
    
