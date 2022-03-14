# File produced automatically by PNCodeGen.ipynb
from scipy.integrate import solve_ivp
import numpy as np
from numpy import dot, cross, log, sqrt, pi
from numpy import euler_gamma as EulerGamma
from numba import jit, njit, float64, boolean
from numba.typed import List
from scipy.interpolate import InterpolatedUnivariateSpline as Spline
from scipy.special import zeta
import quaternionic

qmul = njit(quaternionic.algebra.multiply)
qexp=njit(quaternionic.algebra.exp)
qconj=njit(quaternionic.algebra.conj)
qinverse=njit(quaternionic.algebra.reciprocal)

@njit
def mul(A,B):
    C=np.empty(4)
    qmul(A,B,C)
    return C
    
@njit
def exp(A):
    B=np.empty(4)
    qexp(A,B)
    return B
    
@njit
def conjugate(A):
    B=np.empty(4)
    qconj(A,B)
    return B
    
@njit
def inverse(A):
    B=np.empty(4)
    qinverse(A,B)
    return B

@njit
def FrameFromAngularVelocity_2D_Integrand(rfrak_x, rfrak_y, Omega):
    rfrakMag = np.sqrt(rfrak_x*rfrak_x+rfrak_y*rfrak_y)
    rfrakDot_x = Omega[0]/2.0
    rfrakDot_y = Omega[1]/2.0
    if np.abs(np.sin(rfrakMag)) > 1e-12 and np.abs(np.cos(rfrakMag)) > 1e-12:
        omega_v = (Omega[0]*(-rfrak_y/rfrakMag)+Omega[1]*(rfrak_x/rfrakMag))*np.tan(rfrakMag)-Omega[2]
        Omega[0] += -omega_v*np.sin(2*rfrakMag)*(-rfrak_y/rfrakMag)
        Omega[1] += -omega_v*np.sin(2*rfrakMag)*(rfrak_x/rfrakMag)
        Omega[2] +=  omega_v*np.cos(2*rfrakMag)
        dotTerm = (rfrak_x*Omega[0]+rfrak_y*Omega[1])/(rfrakMag*rfrakMag)
        cotTerm = rfrakMag/(2*np.tan(rfrakMag))
        rfrakDot_x = (Omega[0] - rfrak_x*dotTerm)*cotTerm + rfrak_x*dotTerm/2. - 0.5*Omega[2]*rfrak_y
        rfrakDot_y = (Omega[1] - rfrak_y*dotTerm)*cotTerm + rfrak_y*dotTerm/2. + 0.5*Omega[2]*rfrak_x
    return rfrakDot_x, rfrakDot_y

@njit
def FrameFromAngularVelocityIntegrand(rfrak, Omega):
    rfrakMag = np.sqrt(rfrak[0] * rfrak[0] + rfrak[1] * rfrak[1] + rfrak[2] * rfrak[2])
    OmegaMag = np.sqrt(Omega[0] * Omega[0] + Omega[1] * Omega[1] + Omega[2] * Omega[2])
    # If the matrix is really close to the identity, return
    if rfrakMag < 1e-12*OmegaMag:
        return np.array([Omega[0] / 2.0, Omega[1] / 2.0, Omega[2] / 2.0])
    # If the matrix is really close to singular, it's equivalent to the identity, so return
    if np.abs(np.sin(rfrakMag)) < 1e-12:
        return np.array([Omega[0] / 2.0, Omega[1] / 2.0, Omega[2] / 2.0])
    OmegaOver2 = np.array([Omega[0] / 2.0, Omega[1] / 2.0, Omega[2] / 2.0])
    rfrakHat = np.array([rfrak[0] / rfrakMag, rfrak[1] / rfrakMag, rfrak[2] / rfrakMag])
    return ((OmegaOver2 - rfrakHat * np.dot(rfrakHat, OmegaOver2)) * (rfrakMag / np.tan(rfrakMag))\
        + rfrakHat * np.dot(rfrakHat, OmegaOver2) + np.cross(OmegaOver2, rfrak))  
