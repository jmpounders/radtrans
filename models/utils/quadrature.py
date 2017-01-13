import numpy as np
import numpy.polynomial.polynomial as poly

def PnCoeffs(n) :
    if n==0 :
        return np.array((1))

    pcnm1 = np.array((1))
    pcn = np.array((0, 1))
    for i in range(2,n+1) :
        pcnp1 = poly.polysub((2.0*i-1.0)/i * poly.polymulx(pcn), (i-1.0)/i*pcnm1)
        pcnm1 = pcn
        pcn = pcnp1

    return pcn

def getQuadrature(N) :
    mu = poly.polyroots( PnCoeffs(N) )
    mu = mu[ mu>0.0 ]

    Ntheta = 4*N
    dtheta = 2.0*np.pi/Ntheta
    theta = np.arange(dtheta/2.0, 2*np.pi-dtheta/2.0, dtheta)

    omega_x = np.zeros( N*4.0 )
    omega_y = omega_x.copy()
    omega_z = omega_x.copy()
    wgt = np.ones( N*4.0 ) * 4.0*np.pi/(N*4.0)

    for (i,thetai) in enumerate(theta) :
        omega_x[i] = np.cos(thetai)
        omega_y[i] = np.sin(thetai)
        omega_z[i] = 0.0
            
    return (omega_x,omega_y,omega_z,wgt)

def getQuadrature3DHalf(N) :
    mu = poly.polyroots( PnCoeffs(N) )
    mu = mu[ mu>0.0 ]

    Ntheta = 4*N
    dtheta = 2.0*np.pi/Ntheta
    theta = [dtheta/2.0+dtheta*i for i in range(4*N)]
    theta = np.array(theta)

    omega_x = np.zeros( (N*N*2) )
    omega_y = omega_x.copy()
    omega_z = omega_x.copy()
    wgt = np.ones( (N*N*2) ) * 4.0*np.pi/(2.0*N**2)

    for n in range(N/2) :
        for i in range(N*4) :
            omega_x[ N*4*n + i ] = np.sqrt(1.0-mu[n]**2)*np.cos(theta[i])
            omega_y[ N*4*n + i ] = np.sqrt(1.0-mu[n]**2)*np.sin(theta[i])
            omega_z[ N*4*n + i ] = mu[n]
            
    return (omega_x,omega_y,omega_z,wgt)
