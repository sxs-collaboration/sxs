# File produced automatically by PNCodeGen.ipynb
from scipy.integrate import solve_ivp
import numpy as np
from numpy import dot, cross, log, sqrt, pi
from numpy import euler_gamma as EulerGamma
from numba import jit, njit, float64, boolean
from numba.experimental import jitclass
from scipy.interpolate import InterpolatedUnivariateSpline as Spline
from scipy.special import zeta
import quaternionic

qmul = njit(quaternionic.algebra.multiply)
qexp=njit(quaternionic.algebra.exp)
qconj=njit(quaternionic.algebra.conj)
qinverse=njit(quaternionic.algebra.reciprocal)

@njit(cache=True)
def mul(A,B):
    C=np.empty(4)
    qmul(A,B,C)
    return C
    
@njit(cache=True)
def exp(A):
    B=np.empty(4)
    qexp(A,B)
    return B
    
@njit(cache=True)
def conjugate(A):
    B=np.empty(4)
    qconj(A,B)
    return B
    
@njit(cache=True)
def inverse(A):
    B=np.empty(4)
    qinverse(A,B)
    return B
    
@njit(cache=True)
def normalized(A):
    return A/np.linalg.norm(A)

@njit(cache=True)
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

@njit(cache=True)
def FrameFromAngularVelocityIntegrand(R, Omega):
    return 0.5*mul(np.append(0.0,Omega),R)  

ConsSpec=[('wHat', float64[:]),('xHat', float64[:]),('yHat', float64[:]),('zHat', float64[:]),('M1', float64[:]),('M2', float64[:]),('e_i', float64[:]),('xi_i', float64[:]),('v_i', float64[:]),('S_chi1', float64[:]),('S_chi2', float64[:]),('M', float64[:]),('delta', float64[:]),('nu', float64[:]),('chi1chi1', float64[:]),('chi1chi2', float64[:]),('chi2chi2', float64[:]),('Fcal_0', float64[:]),('Fcal_2', float64[:]),('Fcal_3', float64[:]),('Fcal_4', float64[:]),('Fcal_5', float64[:]),('Fcal_6', float64[:]),('Fcal_lnv_6', float64[:]),('Fcal_7', float64[:]),('Fcal_8', float64[:]),('Fcal_lnv_8', float64[:]),('Fecc_0', float64[:]),('Fecc_2', float64[:]),('E_0', float64[:]),('E_2', float64[:]),('E_4', float64[:]),('E_6', float64[:]),('E_8', float64[:]),('E_lnv_8', float64[:]),('Eecc_2', float64[:]),('EvolveSpin1',boolean),('EvolveSpin2',boolean)]
@jitclass(ConsSpec)
class Cons:
    def __init__(self,wHat,xHat,yHat,zHat,M1,M2,e_i,xi_i,v_i,S_chi1,S_chi2,M,delta,nu,chi1chi1,chi1chi2,chi2chi2,Fcal_0,Fcal_2,Fcal_3,Fcal_4,Fcal_5,Fcal_6,Fcal_lnv_6,Fcal_7,Fcal_8,Fcal_lnv_8,Fecc_0,Fecc_2,E_0,E_2,E_4,E_6,E_8,E_lnv_8,Eecc_2,EvolveSpin1,EvolveSpin2):
        self.wHat=wHat
        self.xHat=xHat
        self.yHat=yHat
        self.zHat=zHat
        self.M1=M1
        self.M2=M2
        self.e_i=e_i
        self.xi_i=xi_i
        self.v_i=v_i
        self.S_chi1=S_chi1
        self.S_chi2=S_chi2
        self.M=M
        self.delta=delta
        self.nu=nu
        self.chi1chi1=chi1chi1
        self.chi1chi2=chi1chi2
        self.chi2chi2=chi2chi2
        self.Fcal_0=Fcal_0
        self.Fcal_2=Fcal_2
        self.Fcal_3=Fcal_3
        self.Fcal_4=Fcal_4
        self.Fcal_5=Fcal_5
        self.Fcal_6=Fcal_6
        self.Fcal_lnv_6=Fcal_lnv_6
        self.Fcal_7=Fcal_7
        self.Fcal_8=Fcal_8
        self.Fcal_lnv_8=Fcal_lnv_8
        self.Fecc_0=Fecc_0
        self.Fecc_2=Fecc_2
        self.E_0=E_0
        self.E_2=E_2
        self.E_4=E_4
        self.E_6=E_6
        self.E_8=E_8
        self.E_lnv_8=E_lnv_8
        self.Eecc_2=Eecc_2
        self.EvolveSpin1=EvolveSpin1
        self.EvolveSpin2=EvolveSpin2

VarsSpec=[('v', float64[:]),('rfrak_chi1', float64[:]),('rfrak_chi2', float64[:]),('Ecc', float64[:]),('xi', float64[:]),('R', float64[:]),('nHat', float64[:]),('lambdaHat', float64[:]),('ellHat', float64[:]),('R_S1', float64[:]),('R_S2', float64[:]),('chiVec1', float64[:]),('chiVec2', float64[:]),('chi1_n', float64[:]),('chi1_lambda', float64[:]),('chi1_ell', float64[:]),('chi2_n', float64[:]),('chi2_lambda', float64[:]),('chi2_ell', float64[:]),('S_ell', float64[:]),('S_n', float64[:]),('S_lambda', float64[:]),('Sigma_ell', float64[:]),('Sigma_n', float64[:]),('Sigma_lambda', float64[:]),('chi_s_ell', float64[:]),('chi_a_ell', float64[:]),('logv', float64[:]),('Fcal_coeff', float64[:]),('Fcal_SQ_4', float64[:]),('Fcal_SO_3', float64[:]),('Fcal_SO_5', float64[:]),('Fcal_SO_6', float64[:]),('Fcal_SO_7', float64[:]),('Fcal_SO_8', float64[:]),('Fecc_coeff', float64[:]),('E_SQ_4', float64[:]),('E_SO_3', float64[:]),('E_SO_5', float64[:]),('E_SO_7', float64[:])]
@jitclass(VarsSpec)
class Vars:
    def __init__(self,v,rfrak_chi1,rfrak_chi2,Ecc,xi,R,nHat,lambdaHat,ellHat,R_S1,R_S2,chiVec1,chiVec2,chi1_n,chi1_lambda,chi1_ell,chi2_n,chi2_lambda,chi2_ell,S_ell,S_n,S_lambda,Sigma_ell,Sigma_n,Sigma_lambda,chi_s_ell,chi_a_ell,logv,Fcal_coeff,Fcal_SQ_4,Fcal_SO_3,Fcal_SO_5,Fcal_SO_6,Fcal_SO_7,Fcal_SO_8,Fecc_coeff,E_SQ_4,E_SO_3,E_SO_5,E_SO_7):
        self.v=v
        self.rfrak_chi1=rfrak_chi1
        self.rfrak_chi2=rfrak_chi2
        self.Ecc=Ecc
        self.xi=xi
        self.R=R
        self.nHat=nHat
        self.lambdaHat=lambdaHat
        self.ellHat=ellHat
        self.R_S1=R_S1
        self.R_S2=R_S2
        self.chiVec1=chiVec1
        self.chiVec2=chiVec2
        self.chi1_n=chi1_n
        self.chi1_lambda=chi1_lambda
        self.chi1_ell=chi1_ell
        self.chi2_n=chi2_n
        self.chi2_lambda=chi2_lambda
        self.chi2_ell=chi2_ell
        self.S_ell=S_ell
        self.S_n=S_n
        self.S_lambda=S_lambda
        self.Sigma_ell=Sigma_ell
        self.Sigma_n=Sigma_n
        self.Sigma_lambda=Sigma_lambda
        self.chi_s_ell=chi_s_ell
        self.chi_a_ell=chi_a_ell
        self.logv=logv
        self.Fcal_coeff=Fcal_coeff
        self.Fcal_SQ_4=Fcal_SQ_4
        self.Fcal_SO_3=Fcal_SO_3
        self.Fcal_SO_5=Fcal_SO_5
        self.Fcal_SO_6=Fcal_SO_6
        self.Fcal_SO_7=Fcal_SO_7
        self.Fcal_SO_8=Fcal_SO_8
        self.Fecc_coeff=Fecc_coeff
        self.E_SQ_4=E_SQ_4
        self.E_SO_3=E_SO_3
        self.E_SO_5=E_SO_5
        self.E_SO_7=E_SO_7

@njit(cache=True)
def Initialization(Cons, wHat_i, xHat_i, yHat_i, zHat_i, M1_i, M2_i, e_i_i, xi_i_i, v_i_i, v_i, S_chi1_i, S_chi2_i): 
    Cons.wHat=wHat_i
    Cons.xHat=xHat_i
    Cons.yHat=yHat_i
    Cons.zHat=zHat_i
    Cons.M1=np.array([M1_i])
    Cons.M2=np.array([M2_i])
    Cons.e_i=np.array([e_i_i])
    Cons.xi_i=np.array([xi_i_i])
    Cons.v_i=np.array([v_i_i])
    Cons.S_chi1=S_chi1_i
    Cons.S_chi2=S_chi2_i
    rfrak_chi1=np.array([0.0,0.0])
    rfrak_chi2=np.array([0.0,0.0])
    Cons.M=Cons.M1 + Cons.M2
    Cons.delta=(Cons.M1 - Cons.M2)/Cons.M
    Cons.nu=Cons.M1*Cons.M2/Cons.M**2
    R_S1=exp(rfrak_chi1[0]*Cons.xHat + rfrak_chi1[1]*Cons.yHat)
    R_S2=exp(rfrak_chi2[0]*Cons.xHat + rfrak_chi2[1]*Cons.yHat)
    chiVec1=mul(mul(mul(Cons.S_chi1,R_S1),Cons.zHat),mul(conjugate(R_S1),conjugate(Cons.S_chi1)))
    chiVec2=mul(mul(mul(Cons.S_chi2,R_S2),Cons.zHat),mul(conjugate(R_S2),conjugate(Cons.S_chi2)))
    Cons.chi1chi1=np.array([dot(chiVec1[1:],chiVec1[1:])])
    Cons.chi1chi2=np.array([dot(chiVec1[1:],chiVec2[1:])])
    Cons.chi2chi2=np.array([dot(chiVec2[1:],chiVec2[1:])])
    Cons.Fcal_0=np.array([1.0])
    Cons.Fcal_2=-35*Cons.nu/12 - 1247/336
    Cons.Fcal_3=np.array([4*pi])
    Cons.Fcal_4=65*Cons.nu**2/18 + 9271*Cons.nu/504 - 44711/9072
    Cons.Fcal_5=pi*(-583*Cons.nu/24 - 8191/672)
    Cons.Fcal_6=-775*Cons.nu**3/324 - 94403*Cons.nu**2/3024 + Cons.nu*(-134543/7776 + 41*pi**2/48) - 1712*log(4)/105 - 1712*EulerGamma/105 + 16*pi**2/3 + 6643739519/69854400
    Cons.Fcal_lnv_6=np.array([-1712/105])
    Cons.Fcal_7=pi*(193385*Cons.nu**2/3024 + 214745*Cons.nu/1728 - 16285/504)
    Cons.Fcal_8=np.array([-1369*pi**2/126 - 323105549467/3178375200 - 47385*log(3)/1568 + 232597*EulerGamma/4410 + 39931*log(2)/294])
    Cons.Fcal_lnv_8=np.array([232597/4410])
    Cons.Fecc_0=Cons.v_i**2*(2833/1008 - 197*Cons.nu/36) + 1.0
    Cons.Fecc_2=-6355*Cons.nu/5652 - 67387/22608
    Cons.E_0=np.array([1.0])
    Cons.E_2=-Cons.nu/12 - 3/4
    Cons.E_4=-Cons.nu**2/24 + 19*Cons.nu/8 - 27/8
    Cons.E_6=-35*Cons.nu**3/5184 - 155*Cons.nu**2/96 + Cons.nu*(34445/576 - 205*pi**2/96) - 675/64
    Cons.E_8=77*Cons.nu**4/31104 + 301*Cons.nu**3/1728 + Cons.nu**2*(-498449/3456 + 3157*pi**2/576) + Cons.nu*(-123671/5760 + 896*EulerGamma/15 + 9037*pi**2/1536 + 1792*log(2)/15) - 3969/128
    Cons.E_lnv_8=896*Cons.nu/15
    Cons.Eecc_2=np.array([-2.0])
    Cons.EvolveSpin1=np.linalg.norm(mul(Cons.S_chi1,conjugate(Cons.S_chi1)))>1e-8
    Cons.EvolveSpin2=np.linalg.norm(mul(Cons.S_chi2,conjugate(Cons.S_chi2)))>1e-8

@njit(cache=True)
def Recalculate_0(Cons,Vars,y):
    Vars.v = np.array([y[0]])
    Vars.rfrak_chi1 = np.array([y[1],y[2]])
    Vars.rfrak_chi2 = np.array([y[3],y[4]])
    Vars.Ecc = Cons.e_i*(Cons.v_i/Vars.v)**(19/6)
    Vars.xi = y[9]+Cons.xi_i
    Vars.R = y[5]*Cons.wHat + y[6]*Cons.xHat + y[7]*Cons.yHat + y[8]*Cons.zHat
    Vars.nHat = mul(mul(Vars.R,Cons.xHat),inverse(Vars.R))
    Vars.ellHat = mul(mul(Vars.R,Cons.zHat),inverse(Vars.R))
    Vars.R_S1 = exp(Vars.rfrak_chi1[0]*Cons.xHat + Vars.rfrak_chi1[1]*Cons.yHat)
    Vars.R_S2 = exp(Vars.rfrak_chi2[0]*Cons.xHat + Vars.rfrak_chi2[1]*Cons.yHat)
    Vars.chiVec1 = mul(mul(mul(Cons.S_chi1,Vars.R_S1),Cons.zHat),mul(conjugate(Vars.R_S1),conjugate(Cons.S_chi1)))
    Vars.chiVec2 = mul(mul(mul(Cons.S_chi2,Vars.R_S2),Cons.zHat),mul(conjugate(Vars.R_S2),conjugate(Cons.S_chi2)))
    Vars.chi1_n = np.array([dot(Vars.chiVec1[1:],Vars.nHat[1:])])
    Vars.chi1_ell = np.array([dot(Vars.chiVec1[1:],Vars.ellHat[1:])])
    Vars.chi2_n = np.array([dot(Vars.chiVec2[1:],Vars.nHat[1:])])
    Vars.chi2_ell = np.array([dot(Vars.chiVec2[1:],Vars.ellHat[1:])])
    Vars.S_ell = Cons.M1**2*Vars.chi1_ell + Cons.M2**2*Vars.chi2_ell
    Vars.S_n = Cons.M1**2*Vars.chi1_n + Cons.M2**2*Vars.chi2_n
    Vars.Sigma_ell = Cons.M*(-Cons.M1*Vars.chi1_ell + Cons.M2*Vars.chi2_ell)
    Vars.Sigma_n = Cons.M*(-Cons.M1*Vars.chi1_n + Cons.M2*Vars.chi2_n)
    Vars.chi_s_ell = Vars.chi1_ell/2 + Vars.chi2_ell/2
    Vars.chi_a_ell = Vars.chi1_ell/2 - Vars.chi2_ell/2
    Vars.Fcal_coeff = 32*Cons.nu**2*Vars.v**10/5
    Vars.Fecc_coeff = 157*Cons.e_i**2.0*Cons.v_i**6.33333333333333/(24*Vars.v**6.33333333333333)

@njit
def OmegaVec_chiVec_1_0(Cons,Vars):
    Omega1_coeff = Vars.v**5/Cons.M
    return Omega1_coeff*Vars.ellHat*(-0.75*Cons.delta + 0.5*Cons.nu + 0.75)

@njit
def OmegaVec_chiVec_2_0(Cons,Vars):
    Omega2_coeff = Vars.v**5/Cons.M
    return Omega2_coeff*Vars.ellHat*(0.75*Cons.delta + 0.5*Cons.nu + 0.75)

@njit
def OmegaVec_0(Cons,Vars):
    gamma_PN_0 = 1.00000000000000
    a_ell_0 = 7.0*Vars.S_n + 3.0*Vars.Sigma_n*Cons.delta
    return Vars.ellHat*Vars.v**3*np.sqrt(1.0 - Vars.Ecc**2.0)*(-Vars.Ecc*np.cos(Vars.xi) + 1.0)**(-2.0)/Cons.M + a_ell_0*gamma_PN_0*Vars.nHat*Vars.v**6/Cons.M**3


@njit(cache=True)
def TaylorT1_0(Cons,Vars):
    Flux = Cons.Fcal_0*Vars.Fcal_coeff + Vars.Fcal_coeff*Cons.Fecc_0*Vars.Fecc_coeff
    dEdV = -Cons.E_0*Cons.M*Cons.nu*Vars.v
    Absorption = 0
    dvdt_T1 = (-Absorption - Flux)/dEdV
    dydt=np.zeros(10)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_0(Cons,Vars)[1:])
    dydt[0] = dvdt_T1[0] 
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_0(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_0(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    dydt[9] = np.linalg.norm(OmegaVec_0(Cons,Vars)[1:])
    return dydt

@njit(cache=True)
def TaylorT4_0(Cons,Vars):
    dvdt_T4 = -2.0*Vars.Fcal_coeff*(-Cons.Fcal_0*Vars.Fcal_coeff - Vars.Fcal_coeff*Cons.Fecc_0*Vars.Fecc_coeff)/Vars.Fcal_coeff/(Cons.nu*Vars.v*2*Cons.E_0*Cons.M)
    dydt=np.zeros(10)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_0(Cons,Vars)[1:])
    dydt[0] = dvdt_T4[0]
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_0(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_0(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    dydt[9] = np.linalg.norm(OmegaVec_0(Cons,Vars)[1:])
    return dydt      

@njit(cache=True)
def TaylorT5_0(Cons,Vars):
    dtdv = -0.5*Cons.nu*Vars.v*2*Cons.E_0*Cons.M/(Vars.Fcal_coeff*(-Cons.Fcal_0*Vars.Fcal_coeff - Vars.Fcal_coeff*Cons.Fecc_0*Vars.Fecc_coeff)/Vars.Fcal_coeff)
    dvdt_T5 = 1.0/dtdv
    dydt=np.zeros(10)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_0(Cons,Vars)[1:])
    dydt[0] = dvdt_T5[0]
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_0(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_0(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    dydt[9] = np.linalg.norm(OmegaVec_0(Cons,Vars)[1:])
    return dydt      

@njit(cache=True)
def Recalculate_0p50(Cons,Vars,y):
    Vars.v = np.array([y[0]])
    Vars.rfrak_chi1 = np.array([y[1],y[2]])
    Vars.rfrak_chi2 = np.array([y[3],y[4]])
    Vars.Ecc = Cons.e_i*(Cons.v_i/Vars.v)**(19/6)
    Vars.xi = y[9]+Cons.xi_i
    Vars.R = y[5]*Cons.wHat + y[6]*Cons.xHat + y[7]*Cons.yHat + y[8]*Cons.zHat
    Vars.nHat = mul(mul(Vars.R,Cons.xHat),inverse(Vars.R))
    Vars.ellHat = mul(mul(Vars.R,Cons.zHat),inverse(Vars.R))
    Vars.R_S1 = exp(Vars.rfrak_chi1[0]*Cons.xHat + Vars.rfrak_chi1[1]*Cons.yHat)
    Vars.R_S2 = exp(Vars.rfrak_chi2[0]*Cons.xHat + Vars.rfrak_chi2[1]*Cons.yHat)
    Vars.chiVec1 = mul(mul(mul(Cons.S_chi1,Vars.R_S1),Cons.zHat),mul(conjugate(Vars.R_S1),conjugate(Cons.S_chi1)))
    Vars.chiVec2 = mul(mul(mul(Cons.S_chi2,Vars.R_S2),Cons.zHat),mul(conjugate(Vars.R_S2),conjugate(Cons.S_chi2)))
    Vars.chi1_n = np.array([dot(Vars.chiVec1[1:],Vars.nHat[1:])])
    Vars.chi1_ell = np.array([dot(Vars.chiVec1[1:],Vars.ellHat[1:])])
    Vars.chi2_n = np.array([dot(Vars.chiVec2[1:],Vars.nHat[1:])])
    Vars.chi2_ell = np.array([dot(Vars.chiVec2[1:],Vars.ellHat[1:])])
    Vars.S_ell = Cons.M1**2*Vars.chi1_ell + Cons.M2**2*Vars.chi2_ell
    Vars.S_n = Cons.M1**2*Vars.chi1_n + Cons.M2**2*Vars.chi2_n
    Vars.Sigma_ell = Cons.M*(-Cons.M1*Vars.chi1_ell + Cons.M2*Vars.chi2_ell)
    Vars.Sigma_n = Cons.M*(-Cons.M1*Vars.chi1_n + Cons.M2*Vars.chi2_n)
    Vars.chi_s_ell = Vars.chi1_ell/2 + Vars.chi2_ell/2
    Vars.chi_a_ell = Vars.chi1_ell/2 - Vars.chi2_ell/2
    Vars.Fcal_coeff = 32*Cons.nu**2*Vars.v**10/5
    Vars.Fecc_coeff = 157*Cons.e_i**2.0*Cons.v_i**6.33333333333333/(24*Vars.v**6.33333333333333)

@njit
def OmegaVec_chiVec_1_0p50(Cons,Vars):
    Omega1_coeff = Vars.v**5/Cons.M
    return Omega1_coeff*(Vars.ellHat*(-0.75*Cons.delta + 0.5*Cons.nu + 0.75) + Vars.nHat*Vars.v*(3.0*Vars.chi1_n*Cons.nu + 3.0*Cons.M2**2*Vars.chi2_n/Cons.M**2) - Cons.M2**2*Vars.chiVec2*Vars.v/Cons.M**2)

@njit
def OmegaVec_chiVec_2_0p50(Cons,Vars):
    Omega2_coeff = Vars.v**5/Cons.M
    return Omega2_coeff*(Vars.ellHat*(0.75*Cons.delta + 0.5*Cons.nu + 0.75) + Vars.nHat*Vars.v*(3.0*Vars.chi2_n*Cons.nu + 3.0*Cons.M1**2*Vars.chi1_n/Cons.M**2) - Cons.M1**2*Vars.chiVec1*Vars.v/Cons.M**2)

@njit
def OmegaVec_0p50(Cons,Vars):
    gamma_PN_0 = 1.00000000000000
    a_ell_0 = 7.0*Vars.S_n + 3.0*Vars.Sigma_n*Cons.delta
    return Vars.ellHat*Vars.v**3*np.sqrt(1.0 - Vars.Ecc**2.0)*(-Vars.Ecc*np.cos(Vars.xi) + 1.0)**(-2.0)/Cons.M + a_ell_0*gamma_PN_0*Vars.nHat*Vars.v**6/Cons.M**3


@njit(cache=True)
def TaylorT1_0p50(Cons,Vars):
    Flux = Cons.Fcal_0*Vars.Fcal_coeff + Vars.Fcal_coeff*Cons.Fecc_0*Vars.Fecc_coeff
    dEdV = -Cons.E_0*Cons.M*Cons.nu*Vars.v
    Absorption = 0
    dvdt_T1 = (-Absorption - Flux)/dEdV
    dydt=np.zeros(10)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_0p50(Cons,Vars)[1:])
    dydt[0] = dvdt_T1[0] 
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_0p50(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_0p50(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    dydt[9] = np.linalg.norm(OmegaVec_0p50(Cons,Vars)[1:])
    return dydt

@njit(cache=True)
def TaylorT4_0p50(Cons,Vars):
    dvdt_T4 = -2.0*Vars.Fcal_coeff*((-Cons.Fcal_0*Vars.Fcal_coeff - Vars.Fcal_coeff*Cons.Fecc_0*Vars.Fecc_coeff)/Vars.Fcal_coeff + 0*Vars.v - 0*(-Cons.Fcal_0*Vars.Fcal_coeff - Vars.Fcal_coeff*Cons.Fecc_0*Vars.Fecc_coeff)/Vars.Fcal_coeff*Vars.v/2*Cons.E_0*Cons.M)/(Cons.nu*Vars.v*2*Cons.E_0*Cons.M)
    dydt=np.zeros(10)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_0p50(Cons,Vars)[1:])
    dydt[0] = dvdt_T4[0]
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_0p50(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_0p50(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    dydt[9] = np.linalg.norm(OmegaVec_0p50(Cons,Vars)[1:])
    return dydt      

@njit(cache=True)
def TaylorT5_0p50(Cons,Vars):
    dtdv = -0.5*Cons.nu*Vars.v*(2*Cons.E_0*Cons.M + 0*Vars.v - 0*2*Cons.E_0*Cons.M*Vars.v/(-Cons.Fcal_0*Vars.Fcal_coeff - Vars.Fcal_coeff*Cons.Fecc_0*Vars.Fecc_coeff)/Vars.Fcal_coeff)/(Vars.Fcal_coeff*(-Cons.Fcal_0*Vars.Fcal_coeff - Vars.Fcal_coeff*Cons.Fecc_0*Vars.Fecc_coeff)/Vars.Fcal_coeff)
    dvdt_T5 = 1.0/dtdv
    dydt=np.zeros(10)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_0p50(Cons,Vars)[1:])
    dydt[0] = dvdt_T5[0]
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_0p50(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_0p50(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    dydt[9] = np.linalg.norm(OmegaVec_0p50(Cons,Vars)[1:])
    return dydt      

@njit(cache=True)
def Recalculate_1p0(Cons,Vars,y):
    Vars.v = np.array([y[0]])
    Vars.rfrak_chi1 = np.array([y[1],y[2]])
    Vars.rfrak_chi2 = np.array([y[3],y[4]])
    Vars.Ecc = Cons.e_i*(Cons.v_i/Vars.v)**(19/6)
    Vars.xi = y[9]+Cons.xi_i
    Vars.R = y[5]*Cons.wHat + y[6]*Cons.xHat + y[7]*Cons.yHat + y[8]*Cons.zHat
    Vars.nHat = mul(mul(Vars.R,Cons.xHat),inverse(Vars.R))
    Vars.ellHat = mul(mul(Vars.R,Cons.zHat),inverse(Vars.R))
    Vars.R_S1 = exp(Vars.rfrak_chi1[0]*Cons.xHat + Vars.rfrak_chi1[1]*Cons.yHat)
    Vars.R_S2 = exp(Vars.rfrak_chi2[0]*Cons.xHat + Vars.rfrak_chi2[1]*Cons.yHat)
    Vars.chiVec1 = mul(mul(mul(Cons.S_chi1,Vars.R_S1),Cons.zHat),mul(conjugate(Vars.R_S1),conjugate(Cons.S_chi1)))
    Vars.chiVec2 = mul(mul(mul(Cons.S_chi2,Vars.R_S2),Cons.zHat),mul(conjugate(Vars.R_S2),conjugate(Cons.S_chi2)))
    Vars.chi1_n = np.array([dot(Vars.chiVec1[1:],Vars.nHat[1:])])
    Vars.chi1_ell = np.array([dot(Vars.chiVec1[1:],Vars.ellHat[1:])])
    Vars.chi2_n = np.array([dot(Vars.chiVec2[1:],Vars.nHat[1:])])
    Vars.chi2_ell = np.array([dot(Vars.chiVec2[1:],Vars.ellHat[1:])])
    Vars.S_ell = Cons.M1**2*Vars.chi1_ell + Cons.M2**2*Vars.chi2_ell
    Vars.S_n = Cons.M1**2*Vars.chi1_n + Cons.M2**2*Vars.chi2_n
    Vars.Sigma_ell = Cons.M*(-Cons.M1*Vars.chi1_ell + Cons.M2*Vars.chi2_ell)
    Vars.Sigma_n = Cons.M*(-Cons.M1*Vars.chi1_n + Cons.M2*Vars.chi2_n)
    Vars.chi_s_ell = Vars.chi1_ell/2 + Vars.chi2_ell/2
    Vars.chi_a_ell = Vars.chi1_ell/2 - Vars.chi2_ell/2
    Vars.Fcal_coeff = 32*Cons.nu**2*Vars.v**10/5
    Vars.Fecc_coeff = 157*Cons.e_i**2.0*Cons.v_i**6.33333333333333/(24*Vars.v**6.33333333333333)

@njit
def OmegaVec_chiVec_1_1p0(Cons,Vars):
    Omega1_coeff = Vars.v**5/Cons.M
    return Omega1_coeff*(Vars.ellHat*(-0.75*Cons.delta + 0.5*Cons.nu + Vars.v**2*(Cons.delta*(0.625*Cons.nu - 0.5625) + Cons.nu*(1.25 - 0.0416666666666667*Cons.nu) + 0.5625) + 0.75) + Vars.nHat*Vars.v*(3.0*Vars.chi1_n*Cons.nu + 3.0*Cons.M2**2*Vars.chi2_n/Cons.M**2) - Cons.M2**2*Vars.chiVec2*Vars.v/Cons.M**2)

@njit
def OmegaVec_chiVec_2_1p0(Cons,Vars):
    Omega2_coeff = Vars.v**5/Cons.M
    return Omega2_coeff*(Vars.ellHat*(0.75*Cons.delta + 0.5*Cons.nu + Vars.v**2*(Cons.delta*(0.5625 - 0.625*Cons.nu) + Cons.nu*(1.25 - 0.0416666666666667*Cons.nu) + 0.5625) + 0.75) + Vars.nHat*Vars.v*(3.0*Vars.chi2_n*Cons.nu + 3.0*Cons.M1**2*Vars.chi1_n/Cons.M**2) - Cons.M1**2*Vars.chiVec1*Vars.v/Cons.M**2)

@njit
def OmegaVec_1p0(Cons,Vars):
    gamma_PN_0 = 1.00000000000000
    a_ell_0 = 7.0*Vars.S_n + 3.0*Vars.Sigma_n*Cons.delta
    a_ell_2 = Vars.S_n*(-9.66666666666667*Cons.nu - 10.0) + Vars.Sigma_n*Cons.delta*(-4.5*Cons.nu - 6.0)
    gamma_PN_2 = 1.0 - 0.333333333333333*Cons.nu
    return Vars.ellHat*Vars.v**3*np.sqrt(1.0 - Vars.Ecc**2.0)*(-Vars.Ecc*np.cos(Vars.xi) + 1.0)**(-2.0)/Cons.M + Vars.nHat*Vars.v**6*(a_ell_0 + a_ell_2*Vars.v**2)*(gamma_PN_0 + gamma_PN_2*Vars.v**2)/Cons.M**3


@njit(cache=True)
def TaylorT1_1p0(Cons,Vars):
    Flux = Vars.Fcal_coeff*Vars.Fecc_coeff*(Cons.Fecc_0 + Cons.Fecc_2*Vars.v**2) + Vars.Fcal_coeff*(Cons.Fcal_0 + Cons.Fcal_2*Vars.v**2)
    dEdV = 1.16666666666667*Cons.Eecc_2*Cons.M*Cons.e_i**2.0*Cons.nu*Vars.v**(-3.33333333333333)*Cons.v_i**6.33333333333333 - Cons.M*Cons.nu*Vars.v*(Cons.E_0 + 2.0*Cons.E_2*Vars.v**2)
    Absorption = 0
    dvdt_T1 = (-Absorption - Flux)/dEdV
    dydt=np.zeros(10)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_1p0(Cons,Vars)[1:])
    dydt[0] = dvdt_T1[0] 
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_1p0(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_1p0(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    dydt[9] = np.linalg.norm(OmegaVec_1p0(Cons,Vars)[1:])
    return dydt

@njit(cache=True)
def TaylorT4_1p0(Cons,Vars):
    dvdt_T4 = -2.0*Vars.Fcal_coeff*(-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff + 0*Vars.v + -Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**2 + (0*(--Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v - 0*Vars.v**2) - 4*Cons.E_2*Cons.M*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**2 + 0**2*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**2/2*Cons.E_0*Cons.M)/2*Cons.E_0*Cons.M)/(Cons.nu*Vars.v*2*Cons.E_0*Cons.M)
    dydt=np.zeros(10)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_1p0(Cons,Vars)[1:])
    dydt[0] = dvdt_T4[0]
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_1p0(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_1p0(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    dydt[9] = np.linalg.norm(OmegaVec_1p0(Cons,Vars)[1:])
    return dydt      

@njit(cache=True)
def TaylorT5_1p0(Cons,Vars):
    dtdv = -0.5*Cons.nu*Vars.v*(2*Cons.E_0*Cons.M + 0*Vars.v + 4*Cons.E_2*Cons.M*Vars.v**2 + (0*(-2*Cons.E_0*Cons.M*Vars.v - 0*Vars.v**2) - -Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*2*Cons.E_0*Cons.M*Vars.v**2 + 0**2*2*Cons.E_0*Cons.M*Vars.v**2/-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff)/-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff)/(Vars.Fcal_coeff*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff)
    dvdt_T5 = 1.0/dtdv
    dydt=np.zeros(10)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_1p0(Cons,Vars)[1:])
    dydt[0] = dvdt_T5[0]
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_1p0(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_1p0(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    dydt[9] = np.linalg.norm(OmegaVec_1p0(Cons,Vars)[1:])
    return dydt      

@njit(cache=True)
def Recalculate_1p5(Cons,Vars,y):
    Vars.v = np.array([y[0]])
    Vars.rfrak_chi1 = np.array([y[1],y[2]])
    Vars.rfrak_chi2 = np.array([y[3],y[4]])
    Vars.Ecc = Cons.e_i*(Cons.v_i/Vars.v)**(19/6)
    Vars.xi = y[9]+Cons.xi_i
    Vars.R = y[5]*Cons.wHat + y[6]*Cons.xHat + y[7]*Cons.yHat + y[8]*Cons.zHat
    Vars.nHat = mul(mul(Vars.R,Cons.xHat),inverse(Vars.R))
    Vars.lambdaHat = mul(mul(Vars.R,Cons.yHat),inverse(Vars.R))
    Vars.ellHat = mul(mul(Vars.R,Cons.zHat),inverse(Vars.R))
    Vars.R_S1 = exp(Vars.rfrak_chi1[0]*Cons.xHat + Vars.rfrak_chi1[1]*Cons.yHat)
    Vars.R_S2 = exp(Vars.rfrak_chi2[0]*Cons.xHat + Vars.rfrak_chi2[1]*Cons.yHat)
    Vars.chiVec1 = mul(mul(mul(Cons.S_chi1,Vars.R_S1),Cons.zHat),mul(conjugate(Vars.R_S1),conjugate(Cons.S_chi1)))
    Vars.chiVec2 = mul(mul(mul(Cons.S_chi2,Vars.R_S2),Cons.zHat),mul(conjugate(Vars.R_S2),conjugate(Cons.S_chi2)))
    Vars.chi1_n = np.array([dot(Vars.chiVec1[1:],Vars.nHat[1:])])
    Vars.chi1_lambda = np.array([dot(Vars.chiVec1[1:],Vars.lambdaHat[1:])])
    Vars.chi1_ell = np.array([dot(Vars.chiVec1[1:],Vars.ellHat[1:])])
    Vars.chi2_n = np.array([dot(Vars.chiVec2[1:],Vars.nHat[1:])])
    Vars.chi2_lambda = np.array([dot(Vars.chiVec2[1:],Vars.lambdaHat[1:])])
    Vars.chi2_ell = np.array([dot(Vars.chiVec2[1:],Vars.ellHat[1:])])
    Vars.S_ell = Cons.M1**2*Vars.chi1_ell + Cons.M2**2*Vars.chi2_ell
    Vars.S_n = Cons.M1**2*Vars.chi1_n + Cons.M2**2*Vars.chi2_n
    Vars.S_lambda = Cons.M1**2*Vars.chi1_lambda + Cons.M2**2*Vars.chi2_lambda
    Vars.Sigma_ell = Cons.M*(-Cons.M1*Vars.chi1_ell + Cons.M2*Vars.chi2_ell)
    Vars.Sigma_n = Cons.M*(-Cons.M1*Vars.chi1_n + Cons.M2*Vars.chi2_n)
    Vars.Sigma_lambda = Cons.M*(-Cons.M1*Vars.chi1_lambda + Cons.M2*Vars.chi2_lambda)
    Vars.chi_s_ell = Vars.chi1_ell/2 + Vars.chi2_ell/2
    Vars.chi_a_ell = Vars.chi1_ell/2 - Vars.chi2_ell/2
    Vars.Fcal_coeff = 32*Cons.nu**2*Vars.v**10/5
    Vars.Fcal_SO_3 = (-4*Vars.S_ell - 5*Vars.Sigma_ell*Cons.delta/4)/Cons.M**2
    Vars.Fecc_coeff = 157*Cons.e_i**2.0*Cons.v_i**6.33333333333333/(24*Vars.v**6.33333333333333)
    Vars.E_SO_3 = (14*Vars.S_ell/3 + 2*Vars.Sigma_ell*Cons.delta)/Cons.M**2

@njit
def OmegaVec_chiVec_1_1p5(Cons,Vars):
    Omega1_coeff = Vars.v**5/Cons.M
    return Omega1_coeff*(Vars.ellHat*(-0.75*Cons.delta + 0.5*Cons.nu + Vars.v**2*(Cons.delta*(0.625*Cons.nu - 0.5625) + Cons.nu*(1.25 - 0.0416666666666667*Cons.nu) + 0.5625) + 0.75) + Vars.nHat*Vars.v*(3.0*Vars.chi1_n*Cons.nu + 3.0*Cons.M2**2*Vars.chi2_n/Cons.M**2) - Cons.M2**2*Vars.chiVec2*Vars.v/Cons.M**2)

@njit
def OmegaVec_chiVec_2_1p5(Cons,Vars):
    Omega2_coeff = Vars.v**5/Cons.M
    return Omega2_coeff*(Vars.ellHat*(0.75*Cons.delta + 0.5*Cons.nu + Vars.v**2*(Cons.delta*(0.5625 - 0.625*Cons.nu) + Cons.nu*(1.25 - 0.0416666666666667*Cons.nu) + 0.5625) + 0.75) + Vars.nHat*Vars.v*(3.0*Vars.chi2_n*Cons.nu + 3.0*Cons.M1**2*Vars.chi1_n/Cons.M**2) - Cons.M1**2*Vars.chiVec1*Vars.v/Cons.M**2)

@njit
def OmegaVec_1p5(Cons,Vars):
    gamma_PN_0 = 1.00000000000000
    a_ell_0 = 7.0*Vars.S_n + 3.0*Vars.Sigma_n*Cons.delta
    gamma_PN_3 = (1.66666666666667*Vars.S_ell + Vars.Sigma_ell*Cons.delta)/Cons.M**2
    a_ell_2 = Vars.S_n*(-9.66666666666667*Cons.nu - 10.0) + Vars.Sigma_n*Cons.delta*(-4.5*Cons.nu - 6.0)
    gamma_PN_2 = 1.0 - 0.333333333333333*Cons.nu
    return Vars.ellHat*Vars.v**3*np.sqrt(1.0 - Vars.Ecc**2.0)*(-Vars.Ecc*np.cos(Vars.xi) + 1.0)**(-2.0)/Cons.M + Vars.nHat*Vars.v**6*(a_ell_0 + a_ell_2*Vars.v**2)*(gamma_PN_0 + Vars.v**2*(gamma_PN_2 + gamma_PN_3*Vars.v))/Cons.M**3


@njit(cache=True)
def TaylorT1_1p5(Cons,Vars):
    Flux = Vars.Fcal_coeff*Vars.Fecc_coeff*(Cons.Fecc_0 + Cons.Fecc_2*Vars.v**2) + Vars.Fcal_coeff*(Cons.Fcal_0 + Vars.v**2*(Cons.Fcal_2 + Vars.v*(Cons.Fcal_3 + Vars.Fcal_SO_3)))
    dEdV = 1.16666666666667*Cons.Eecc_2*Cons.M*Cons.e_i**2.0*Cons.nu*Vars.v**(-3.33333333333333)*Cons.v_i**6.33333333333333 - 0.5*Cons.M*Cons.nu*Vars.v*(2.0*Cons.E_0 + Vars.v**2*(4.0*Cons.E_2 + 5.0*Vars.E_SO_3*Vars.v))
    Absorption = 0
    dvdt_T1 = (-Absorption - Flux)/dEdV
    dydt=np.zeros(10)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_1p5(Cons,Vars)[1:])
    dydt[0] = dvdt_T1[0] 
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_1p5(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_1p5(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    dydt[9] = np.linalg.norm(OmegaVec_1p5(Cons,Vars)[1:])
    return dydt

@njit(cache=True)
def TaylorT4_1p5(Cons,Vars):
    dvdt_T4 = -2.0*Vars.Fcal_coeff*(1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff + 1.0*0*Vars.v + 1.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**2 + 1.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*Vars.v**3 + (0*(-1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v - 1.0*0*Vars.v**2 - 1.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**3) + 4*Cons.E_2*Cons.M*(-1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**2 - 1.0*0*Vars.v**3) - 1.0*5*Vars.E_SO_3*Cons.M*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**3 + (0*(0*(1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**2 + 1.0*0*Vars.v**3) + 2.0*4*Cons.E_2*Cons.M*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**3) - 1.0*0**3*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**3/2*Cons.E_0*Cons.M)/2*Cons.E_0*Cons.M)/2*Cons.E_0*Cons.M)/(Cons.nu*Vars.v*2*Cons.E_0*Cons.M)
    dydt=np.zeros(10)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_1p5(Cons,Vars)[1:])
    dydt[0] = dvdt_T4[0]
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_1p5(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_1p5(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    dydt[9] = np.linalg.norm(OmegaVec_1p5(Cons,Vars)[1:])
    return dydt      

@njit(cache=True)
def TaylorT5_1p5(Cons,Vars):
    dtdv = -0.5*Cons.nu*Vars.v*(1.0*2*Cons.E_0*Cons.M + 1.0*0*Vars.v + 1.0*4*Cons.E_2*Cons.M*Vars.v**2 + 1.0*5*Vars.E_SO_3*Cons.M*Vars.v**3 + (0*(-1.0*2*Cons.E_0*Cons.M*Vars.v - 1.0*0*Vars.v**2 - 1.0*4*Cons.E_2*Cons.M*Vars.v**3) + -Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*(-1.0*2*Cons.E_0*Cons.M*Vars.v**2 - 1.0*0*Vars.v**3) - 1.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*2*Cons.E_0*Cons.M*Vars.v**3 + (0*(0*(1.0*2*Cons.E_0*Cons.M*Vars.v**2 + 1.0*0*Vars.v**3) + 2.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*2*Cons.E_0*Cons.M*Vars.v**3) - 1.0*0**3*2*Cons.E_0*Cons.M*Vars.v**3/-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff)/-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff)/-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff)/(Vars.Fcal_coeff*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff)
    dvdt_T5 = 1.0/dtdv
    dydt=np.zeros(10)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_1p5(Cons,Vars)[1:])
    dydt[0] = dvdt_T5[0]
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_1p5(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_1p5(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    dydt[9] = np.linalg.norm(OmegaVec_1p5(Cons,Vars)[1:])
    return dydt      

@njit(cache=True)
def Recalculate_2p0(Cons,Vars,y):
    Vars.v = np.array([y[0]])
    Vars.rfrak_chi1 = np.array([y[1],y[2]])
    Vars.rfrak_chi2 = np.array([y[3],y[4]])
    Vars.Ecc = Cons.e_i*(Cons.v_i/Vars.v)**(19/6)
    Vars.xi = y[9]+Cons.xi_i
    Vars.R = y[5]*Cons.wHat + y[6]*Cons.xHat + y[7]*Cons.yHat + y[8]*Cons.zHat
    Vars.nHat = mul(mul(Vars.R,Cons.xHat),inverse(Vars.R))
    Vars.lambdaHat = mul(mul(Vars.R,Cons.yHat),inverse(Vars.R))
    Vars.ellHat = mul(mul(Vars.R,Cons.zHat),inverse(Vars.R))
    Vars.R_S1 = exp(Vars.rfrak_chi1[0]*Cons.xHat + Vars.rfrak_chi1[1]*Cons.yHat)
    Vars.R_S2 = exp(Vars.rfrak_chi2[0]*Cons.xHat + Vars.rfrak_chi2[1]*Cons.yHat)
    Vars.chiVec1 = mul(mul(mul(Cons.S_chi1,Vars.R_S1),Cons.zHat),mul(conjugate(Vars.R_S1),conjugate(Cons.S_chi1)))
    Vars.chiVec2 = mul(mul(mul(Cons.S_chi2,Vars.R_S2),Cons.zHat),mul(conjugate(Vars.R_S2),conjugate(Cons.S_chi2)))
    Vars.chi1_n = np.array([dot(Vars.chiVec1[1:],Vars.nHat[1:])])
    Vars.chi1_lambda = np.array([dot(Vars.chiVec1[1:],Vars.lambdaHat[1:])])
    Vars.chi1_ell = np.array([dot(Vars.chiVec1[1:],Vars.ellHat[1:])])
    Vars.chi2_n = np.array([dot(Vars.chiVec2[1:],Vars.nHat[1:])])
    Vars.chi2_lambda = np.array([dot(Vars.chiVec2[1:],Vars.lambdaHat[1:])])
    Vars.chi2_ell = np.array([dot(Vars.chiVec2[1:],Vars.ellHat[1:])])
    Vars.S_ell = Cons.M1**2*Vars.chi1_ell + Cons.M2**2*Vars.chi2_ell
    Vars.S_n = Cons.M1**2*Vars.chi1_n + Cons.M2**2*Vars.chi2_n
    Vars.S_lambda = Cons.M1**2*Vars.chi1_lambda + Cons.M2**2*Vars.chi2_lambda
    Vars.Sigma_ell = Cons.M*(-Cons.M1*Vars.chi1_ell + Cons.M2*Vars.chi2_ell)
    Vars.Sigma_n = Cons.M*(-Cons.M1*Vars.chi1_n + Cons.M2*Vars.chi2_n)
    Vars.Sigma_lambda = Cons.M*(-Cons.M1*Vars.chi1_lambda + Cons.M2*Vars.chi2_lambda)
    Vars.chi_s_ell = Vars.chi1_ell/2 + Vars.chi2_ell/2
    Vars.chi_a_ell = Vars.chi1_ell/2 - Vars.chi2_ell/2
    Vars.Fcal_coeff = 32*Cons.nu**2*Vars.v**10/5
    Vars.Fcal_SQ_4 = Cons.chi1chi1*(-89*Cons.delta/192 + 89*Cons.nu/96 - 89/192) - 103*Cons.chi1chi2*Cons.nu/48 + Cons.chi2chi2*(89*Cons.delta/192 + 89*Cons.nu/96 - 89/192) + Vars.chi_a_ell*(Vars.chi_a_ell*(287/96 - 12*Cons.nu) + 287*Vars.chi_s_ell*Cons.delta/48) + Vars.chi_s_ell**2*(Cons.nu/24 + 287/96)
    Vars.Fcal_SO_3 = (-4*Vars.S_ell - 5*Vars.Sigma_ell*Cons.delta/4)/Cons.M**2
    Vars.Fecc_coeff = 157*Cons.e_i**2.0*Cons.v_i**6.33333333333333/(24*Vars.v**6.33333333333333)
    Vars.E_SQ_4 = -3*Vars.chi_a_ell**2/2 - 3*Vars.chi_s_ell**2/2 - Cons.delta*(Cons.chi2chi2/2 + 3*Vars.chi_a_ell*Vars.chi_s_ell) + Cons.nu*(Cons.chi1chi2 + 6*Vars.chi_a_ell**2) + (Cons.chi1chi1 + Cons.chi2chi2)*(Cons.delta - 2*Cons.nu + 1)/4
    Vars.E_SO_3 = (14*Vars.S_ell/3 + 2*Vars.Sigma_ell*Cons.delta)/Cons.M**2

@njit
def OmegaVec_chiVec_1_2p0(Cons,Vars):
    Omega1_coeff = Vars.v**5/Cons.M
    return Omega1_coeff*(Vars.ellHat*(-0.75*Cons.delta + 0.5*Cons.nu + Vars.v**2*(Cons.delta*(0.625*Cons.nu - 0.5625) + Cons.nu*(1.25 - 0.0416666666666667*Cons.nu) + Vars.v**2*(Cons.delta*(Cons.nu*(4.875 - 0.15625*Cons.nu) - 0.84375) + Cons.nu*(Cons.nu*(-0.0208333333333333*Cons.nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + Vars.nHat*Vars.v*(3.0*Vars.chi1_n*Cons.nu + 3.0*Cons.M2**2*Vars.chi2_n/Cons.M**2) - Cons.M2**2*Vars.chiVec2*Vars.v/Cons.M**2)

@njit
def OmegaVec_chiVec_2_2p0(Cons,Vars):
    Omega2_coeff = Vars.v**5/Cons.M
    return Omega2_coeff*(Vars.ellHat*(0.75*Cons.delta + 0.5*Cons.nu + Vars.v**2*(Cons.delta*(0.5625 - 0.625*Cons.nu) + Cons.nu*(1.25 - 0.0416666666666667*Cons.nu) + Vars.v**2*(Cons.delta*(Cons.nu*(0.15625*Cons.nu - 4.875) + 0.84375) + Cons.nu*(Cons.nu*(-0.0208333333333333*Cons.nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + Vars.nHat*Vars.v*(3.0*Vars.chi2_n*Cons.nu + 3.0*Cons.M1**2*Vars.chi1_n/Cons.M**2) - Cons.M1**2*Vars.chiVec1*Vars.v/Cons.M**2)

@njit
def OmegaVec_2p0(Cons,Vars):
    gamma_PN_0 = 1.00000000000000
    a_ell_0 = 7.0*Vars.S_n + 3.0*Vars.Sigma_n*Cons.delta
    gamma_PN_3 = (1.66666666666667*Vars.S_ell + Vars.Sigma_ell*Cons.delta)/Cons.M**2
    a_ell_4 = Vars.S_n*(5.77777777777778*Cons.nu**2 + 14.75*Cons.nu + 1.5) + Vars.Sigma_n*Cons.delta*(2.83333333333333*Cons.nu**2 + 9.125*Cons.nu + 1.5)
    a_ell_2 = Vars.S_n*(-9.66666666666667*Cons.nu - 10.0) + Vars.Sigma_n*Cons.delta*(-4.5*Cons.nu - 6.0)
    gamma_PN_4 = 1.0 - 5.41666666666667*Cons.nu
    gamma_PN_2 = 1.0 - 0.333333333333333*Cons.nu
    return Vars.ellHat*Vars.v**3*np.sqrt(1.0 - Vars.Ecc**2.0)*(-Vars.Ecc*np.cos(Vars.xi) + 1.0)**(-2.0)/Cons.M + Vars.nHat*Vars.v**6*(a_ell_0 + Vars.v**2*(a_ell_2 + a_ell_4*Vars.v**2))*(gamma_PN_0 + Vars.v**2*(gamma_PN_2 + Vars.v*(gamma_PN_3 + gamma_PN_4*Vars.v)))/Cons.M**3


@njit(cache=True)
def TaylorT1_2p0(Cons,Vars):
    Flux = Vars.Fcal_coeff*Vars.Fecc_coeff*(Cons.Fecc_0 + Cons.Fecc_2*Vars.v**2) + Vars.Fcal_coeff*(Cons.Fcal_0 + Vars.v**2*(Cons.Fcal_2 + Vars.v*(Cons.Fcal_3 + Vars.Fcal_SO_3 + Vars.v*(Cons.Fcal_4 + Vars.Fcal_SQ_4))))
    dEdV = 1.16666666666667*Cons.Eecc_2*Cons.M*Cons.e_i**2.0*Cons.nu*Vars.v**(-3.33333333333333)*Cons.v_i**6.33333333333333 - 0.5*Cons.M*Cons.nu*Vars.v*(2.0*Cons.E_0 + Vars.v**2*(4.0*Cons.E_2 + Vars.v*(5.0*Vars.E_SO_3 + 6.0*Vars.v*(Cons.E_4 + Vars.E_SQ_4))))
    Absorption = 0
    dvdt_T1 = (-Absorption - Flux)/dEdV
    dydt=np.zeros(10)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_2p0(Cons,Vars)[1:])
    dydt[0] = dvdt_T1[0] 
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_2p0(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_2p0(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    dydt[9] = np.linalg.norm(OmegaVec_2p0(Cons,Vars)[1:])
    return dydt

@njit(cache=True)
def TaylorT4_2p0(Cons,Vars):
    dvdt_T4 = -2.0*Vars.Fcal_coeff*(1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff + 1.0*0*Vars.v + 1.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**2 + 1.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*Vars.v**3 + 1.0*-Cons.Fcal_4 - Vars.Fcal_SQ_4*Vars.v**4 + (0*(-1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v - 1.0*0*Vars.v**2 - 1.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**3 - 1.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*Vars.v**4) + 4*Cons.E_2*Cons.M*(-1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**2 - 1.0*0*Vars.v**3 - 1.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**4) + 5*Vars.E_SO_3*Cons.M*(-1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**3 - 1.0*0*Vars.v**4) - 1.0*6*Cons.E_4*Cons.M + 6*Vars.E_SQ_4*Cons.M*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**4 + (0*(0*(1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**2 + 1.0*0*Vars.v**3 + 1.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**4) + 4*Cons.E_2*Cons.M*(2.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**3 + 2.0*0*Vars.v**4) + 2.0*5*Vars.E_SO_3*Cons.M*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**4) + 1.0*4*Cons.E_2*Cons.M**2*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**4 + (0**2*(0*(-1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**3 - 1.0*0*Vars.v**4) - 3.0*4*Cons.E_2*Cons.M*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**4) + 1.0*0**4*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**4/2*Cons.E_0*Cons.M)/2*Cons.E_0*Cons.M)/2*Cons.E_0*Cons.M)/2*Cons.E_0*Cons.M)/(Cons.nu*Vars.v*2*Cons.E_0*Cons.M)
    dydt=np.zeros(10)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_2p0(Cons,Vars)[1:])
    dydt[0] = dvdt_T4[0]
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_2p0(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_2p0(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    dydt[9] = np.linalg.norm(OmegaVec_2p0(Cons,Vars)[1:])
    return dydt      

@njit(cache=True)
def TaylorT5_2p0(Cons,Vars):
    dtdv = -0.5*Cons.nu*Vars.v*(1.0*2*Cons.E_0*Cons.M + 1.0*0*Vars.v + 1.0*4*Cons.E_2*Cons.M*Vars.v**2 + 1.0*5*Vars.E_SO_3*Cons.M*Vars.v**3 + 1.0*6*Cons.E_4*Cons.M + 6*Vars.E_SQ_4*Cons.M*Vars.v**4 + (0*(-1.0*2*Cons.E_0*Cons.M*Vars.v - 1.0*0*Vars.v**2 - 1.0*4*Cons.E_2*Cons.M*Vars.v**3 - 1.0*5*Vars.E_SO_3*Cons.M*Vars.v**4) + -Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*(-1.0*2*Cons.E_0*Cons.M*Vars.v**2 - 1.0*0*Vars.v**3 - 1.0*4*Cons.E_2*Cons.M*Vars.v**4) + -Cons.Fcal_3 - Vars.Fcal_SO_3*(-1.0*2*Cons.E_0*Cons.M*Vars.v**3 - 1.0*0*Vars.v**4) - 1.0*-Cons.Fcal_4 - Vars.Fcal_SQ_4*2*Cons.E_0*Cons.M*Vars.v**4 + (0*(0*(1.0*2*Cons.E_0*Cons.M*Vars.v**2 + 1.0*0*Vars.v**3 + 1.0*4*Cons.E_2*Cons.M*Vars.v**4) + -Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*(2.0*2*Cons.E_0*Cons.M*Vars.v**3 + 2.0*0*Vars.v**4) + 2.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*2*Cons.E_0*Cons.M*Vars.v**4) + 1.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff**2*2*Cons.E_0*Cons.M*Vars.v**4 + (0**2*(0*(-1.0*2*Cons.E_0*Cons.M*Vars.v**3 - 1.0*0*Vars.v**4) - 3.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*2*Cons.E_0*Cons.M*Vars.v**4) + 1.0*0**4*2*Cons.E_0*Cons.M*Vars.v**4/-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff)/-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff)/-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff)/-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff)/(Vars.Fcal_coeff*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff)
    dvdt_T5 = 1.0/dtdv
    dydt=np.zeros(10)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_2p0(Cons,Vars)[1:])
    dydt[0] = dvdt_T5[0]
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_2p0(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_2p0(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    dydt[9] = np.linalg.norm(OmegaVec_2p0(Cons,Vars)[1:])
    return dydt      

@njit(cache=True)
def Recalculate_2p5(Cons,Vars,y):
    Vars.v = np.array([y[0]])
    Vars.rfrak_chi1 = np.array([y[1],y[2]])
    Vars.rfrak_chi2 = np.array([y[3],y[4]])
    Vars.Ecc = Cons.e_i*(Cons.v_i/Vars.v)**(19/6)
    Vars.xi = y[9]+Cons.xi_i
    Vars.R = y[5]*Cons.wHat + y[6]*Cons.xHat + y[7]*Cons.yHat + y[8]*Cons.zHat
    Vars.nHat = mul(mul(Vars.R,Cons.xHat),inverse(Vars.R))
    Vars.lambdaHat = mul(mul(Vars.R,Cons.yHat),inverse(Vars.R))
    Vars.ellHat = mul(mul(Vars.R,Cons.zHat),inverse(Vars.R))
    Vars.R_S1 = exp(Vars.rfrak_chi1[0]*Cons.xHat + Vars.rfrak_chi1[1]*Cons.yHat)
    Vars.R_S2 = exp(Vars.rfrak_chi2[0]*Cons.xHat + Vars.rfrak_chi2[1]*Cons.yHat)
    Vars.chiVec1 = mul(mul(mul(Cons.S_chi1,Vars.R_S1),Cons.zHat),mul(conjugate(Vars.R_S1),conjugate(Cons.S_chi1)))
    Vars.chiVec2 = mul(mul(mul(Cons.S_chi2,Vars.R_S2),Cons.zHat),mul(conjugate(Vars.R_S2),conjugate(Cons.S_chi2)))
    Vars.chi1_n = np.array([dot(Vars.chiVec1[1:],Vars.nHat[1:])])
    Vars.chi1_lambda = np.array([dot(Vars.chiVec1[1:],Vars.lambdaHat[1:])])
    Vars.chi1_ell = np.array([dot(Vars.chiVec1[1:],Vars.ellHat[1:])])
    Vars.chi2_n = np.array([dot(Vars.chiVec2[1:],Vars.nHat[1:])])
    Vars.chi2_lambda = np.array([dot(Vars.chiVec2[1:],Vars.lambdaHat[1:])])
    Vars.chi2_ell = np.array([dot(Vars.chiVec2[1:],Vars.ellHat[1:])])
    Vars.S_ell = Cons.M1**2*Vars.chi1_ell + Cons.M2**2*Vars.chi2_ell
    Vars.S_n = Cons.M1**2*Vars.chi1_n + Cons.M2**2*Vars.chi2_n
    Vars.S_lambda = Cons.M1**2*Vars.chi1_lambda + Cons.M2**2*Vars.chi2_lambda
    Vars.Sigma_ell = Cons.M*(-Cons.M1*Vars.chi1_ell + Cons.M2*Vars.chi2_ell)
    Vars.Sigma_n = Cons.M*(-Cons.M1*Vars.chi1_n + Cons.M2*Vars.chi2_n)
    Vars.Sigma_lambda = Cons.M*(-Cons.M1*Vars.chi1_lambda + Cons.M2*Vars.chi2_lambda)
    Vars.chi_s_ell = Vars.chi1_ell/2 + Vars.chi2_ell/2
    Vars.chi_a_ell = Vars.chi1_ell/2 - Vars.chi2_ell/2
    Vars.Fcal_coeff = 32*Cons.nu**2*Vars.v**10/5
    Vars.Fcal_SQ_4 = Cons.chi1chi1*(-89*Cons.delta/192 + 89*Cons.nu/96 - 89/192) - 103*Cons.chi1chi2*Cons.nu/48 + Cons.chi2chi2*(89*Cons.delta/192 + 89*Cons.nu/96 - 89/192) + Vars.chi_a_ell*(Vars.chi_a_ell*(287/96 - 12*Cons.nu) + 287*Vars.chi_s_ell*Cons.delta/48) + Vars.chi_s_ell**2*(Cons.nu/24 + 287/96)
    Vars.Fcal_SO_3 = (-4*Vars.S_ell - 5*Vars.Sigma_ell*Cons.delta/4)/Cons.M**2
    Vars.Fcal_SO_5 = (Vars.S_ell*(272*Cons.nu/9 - 9/2) + Vars.Sigma_ell*Cons.delta*(43*Cons.nu/4 - 13/16))/Cons.M**2
    Vars.Fecc_coeff = 157*Cons.e_i**2.0*Cons.v_i**6.33333333333333/(24*Vars.v**6.33333333333333)
    Vars.E_SQ_4 = -3*Vars.chi_a_ell**2/2 - 3*Vars.chi_s_ell**2/2 - Cons.delta*(Cons.chi2chi2/2 + 3*Vars.chi_a_ell*Vars.chi_s_ell) + Cons.nu*(Cons.chi1chi2 + 6*Vars.chi_a_ell**2) + (Cons.chi1chi1 + Cons.chi2chi2)*(Cons.delta - 2*Cons.nu + 1)/4
    Vars.E_SO_3 = (14*Vars.S_ell/3 + 2*Vars.Sigma_ell*Cons.delta)/Cons.M**2
    Vars.E_SO_5 = (Vars.S_ell*(11 - 61*Cons.nu/9) + Vars.Sigma_ell*Cons.delta*(3 - 10*Cons.nu/3))/Cons.M**2

@njit
def OmegaVec_chiVec_1_2p5(Cons,Vars):
    Omega1_coeff = Vars.v**5/Cons.M
    return Omega1_coeff*(Vars.ellHat*(-0.75*Cons.delta + 0.5*Cons.nu + Vars.v**2*(Cons.delta*(0.625*Cons.nu - 0.5625) + Cons.nu*(1.25 - 0.0416666666666667*Cons.nu) + Vars.v**2*(Cons.delta*(Cons.nu*(4.875 - 0.15625*Cons.nu) - 0.84375) + Cons.nu*(Cons.nu*(-0.0208333333333333*Cons.nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + Vars.nHat*Vars.v*(3.0*Vars.chi1_n*Cons.nu + 3.0*Cons.M2**2*Vars.chi2_n/Cons.M**2) - Cons.M2**2*Vars.chiVec2*Vars.v/Cons.M**2)

@njit
def OmegaVec_chiVec_2_2p5(Cons,Vars):
    Omega2_coeff = Vars.v**5/Cons.M
    return Omega2_coeff*(Vars.ellHat*(0.75*Cons.delta + 0.5*Cons.nu + Vars.v**2*(Cons.delta*(0.5625 - 0.625*Cons.nu) + Cons.nu*(1.25 - 0.0416666666666667*Cons.nu) + Vars.v**2*(Cons.delta*(Cons.nu*(0.15625*Cons.nu - 4.875) + 0.84375) + Cons.nu*(Cons.nu*(-0.0208333333333333*Cons.nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + Vars.nHat*Vars.v*(3.0*Vars.chi2_n*Cons.nu + 3.0*Cons.M1**2*Vars.chi1_n/Cons.M**2) - Cons.M1**2*Vars.chiVec1*Vars.v/Cons.M**2)

@njit
def OmegaVec_2p5(Cons,Vars):
    gamma_PN_0 = 1.00000000000000
    a_ell_0 = 7.0*Vars.S_n + 3.0*Vars.Sigma_n*Cons.delta
    gamma_PN_3 = (1.66666666666667*Vars.S_ell + Vars.Sigma_ell*Cons.delta)/Cons.M**2
    a_ell_4 = Vars.S_n*(5.77777777777778*Cons.nu**2 + 14.75*Cons.nu + 1.5) + Vars.Sigma_n*Cons.delta*(2.83333333333333*Cons.nu**2 + 9.125*Cons.nu + 1.5)
    a_ell_2 = Vars.S_n*(-9.66666666666667*Cons.nu - 10.0) + Vars.Sigma_n*Cons.delta*(-4.5*Cons.nu - 6.0)
    gamma_PN_5 = (Vars.S_ell*(0.888888888888889*Cons.nu + 3.33333333333333) + 2.0*Vars.Sigma_ell*Cons.delta)/Cons.M**2
    gamma_PN_4 = 1.0 - 5.41666666666667*Cons.nu
    gamma_PN_2 = 1.0 - 0.333333333333333*Cons.nu
    return Vars.ellHat*Vars.v**3*np.sqrt(1.0 - Vars.Ecc**2.0)*(-Vars.Ecc*np.cos(Vars.xi) + 1.0)**(-2.0)/Cons.M + Vars.nHat*Vars.v**6*(a_ell_0 + Vars.v**2*(a_ell_2 + a_ell_4*Vars.v**2))*(gamma_PN_0 + Vars.v**2*(gamma_PN_2 + Vars.v*(gamma_PN_3 + Vars.v*(gamma_PN_4 + gamma_PN_5*Vars.v))))/Cons.M**3


@njit(cache=True)
def TaylorT1_2p5(Cons,Vars):
    Flux = Vars.Fcal_coeff*Vars.Fecc_coeff*(Cons.Fecc_0 + Cons.Fecc_2*Vars.v**2) + Vars.Fcal_coeff*(Cons.Fcal_0 + Vars.v**2*(Cons.Fcal_2 + Vars.v*(Cons.Fcal_3 + Vars.Fcal_SO_3 + Vars.v*(Cons.Fcal_4 + Vars.Fcal_SQ_4 + Vars.v*(Cons.Fcal_5 + Vars.Fcal_SO_5)))))
    dEdV = 1.16666666666667*Cons.Eecc_2*Cons.M*Cons.e_i**2.0*Cons.nu*Vars.v**(-3.33333333333333)*Cons.v_i**6.33333333333333 - 0.5*Cons.M*Cons.nu*Vars.v*(2.0*Cons.E_0 + Vars.v**2*(4.0*Cons.E_2 + Vars.v*(5.0*Vars.E_SO_3 + Vars.v*(6.0*Cons.E_4 + 7.0*Vars.E_SO_5*Vars.v + 6.0*Vars.E_SQ_4))))
    Absorption = 0
    dvdt_T1 = (-Absorption - Flux)/dEdV
    dydt=np.zeros(10)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_2p5(Cons,Vars)[1:])
    dydt[0] = dvdt_T1[0] 
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_2p5(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_2p5(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    dydt[9] = np.linalg.norm(OmegaVec_2p5(Cons,Vars)[1:])
    return dydt

@njit(cache=True)
def TaylorT4_2p5(Cons,Vars):
    dvdt_T4 = -2.0*Vars.Fcal_coeff*(1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff + 1.0*0*Vars.v + 1.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**2 + 1.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*Vars.v**3 + 1.0*-Cons.Fcal_4 - Vars.Fcal_SQ_4*Vars.v**4 + 1.0*-Cons.Fcal_5 - Vars.Fcal_SO_5*Vars.v**5 + (0*(-1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v - 1.0*0*Vars.v**2 - 1.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**3 - 1.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*Vars.v**4 - 1.0*-Cons.Fcal_4 - Vars.Fcal_SQ_4*Vars.v**5) + 4*Cons.E_2*Cons.M*(-1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**2 - 1.0*0*Vars.v**3 - 1.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**4 - 1.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*Vars.v**5) + 5*Vars.E_SO_3*Cons.M*(-1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**3 - 1.0*0*Vars.v**4 - 1.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**5) + 6*Cons.E_4*Cons.M + 6*Vars.E_SQ_4*Cons.M*(-1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**4 - 1.0*0*Vars.v**5) - 1.0*7*Vars.E_SO_5*Cons.M*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**5 + (0*(0*(1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**2 + 1.0*0*Vars.v**3 + 1.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**4 + 1.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*Vars.v**5) + 4*Cons.E_2*Cons.M*(2.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**3 + 2.0*0*Vars.v**4 + 2.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**5) + 5*Vars.E_SO_3*Cons.M*(2.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**4 + 2.0*0*Vars.v**5) + 2.0*6*Cons.E_4*Cons.M + 6*Vars.E_SQ_4*Cons.M*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**5) + 4*Cons.E_2*Cons.M*(4*Cons.E_2*Cons.M*(1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**4 + 1.0*0*Vars.v**5) + 2.0*5*Vars.E_SO_3*Cons.M*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**5) + (0*(0*(0*(-1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**3 - 1.0*0*Vars.v**4 - 1.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**5) + 4*Cons.E_2*Cons.M*(-3.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**4 - 3.0*0*Vars.v**5) - 3.0*5*Vars.E_SO_3*Cons.M*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**5) - 3.0*4*Cons.E_2*Cons.M**2*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**5) + (0**3*(0*(1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**4 + 1.0*0*Vars.v**5) + 4.0*4*Cons.E_2*Cons.M*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**5) - 1.0*0**5*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**5/2*Cons.E_0*Cons.M)/2*Cons.E_0*Cons.M)/2*Cons.E_0*Cons.M)/2*Cons.E_0*Cons.M)/2*Cons.E_0*Cons.M)/(Cons.nu*Vars.v*2*Cons.E_0*Cons.M)
    dydt=np.zeros(10)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_2p5(Cons,Vars)[1:])
    dydt[0] = dvdt_T4[0]
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_2p5(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_2p5(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    dydt[9] = np.linalg.norm(OmegaVec_2p5(Cons,Vars)[1:])
    return dydt      

@njit(cache=True)
def TaylorT5_2p5(Cons,Vars):
    dtdv = -0.5*Cons.nu*Vars.v*(1.0*2*Cons.E_0*Cons.M + 1.0*0*Vars.v + 1.0*4*Cons.E_2*Cons.M*Vars.v**2 + 1.0*5*Vars.E_SO_3*Cons.M*Vars.v**3 + 1.0*6*Cons.E_4*Cons.M + 6*Vars.E_SQ_4*Cons.M*Vars.v**4 + 1.0*7*Vars.E_SO_5*Cons.M*Vars.v**5 + (0*(-1.0*2*Cons.E_0*Cons.M*Vars.v - 1.0*0*Vars.v**2 - 1.0*4*Cons.E_2*Cons.M*Vars.v**3 - 1.0*5*Vars.E_SO_3*Cons.M*Vars.v**4 - 1.0*6*Cons.E_4*Cons.M + 6*Vars.E_SQ_4*Cons.M*Vars.v**5) + -Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*(-1.0*2*Cons.E_0*Cons.M*Vars.v**2 - 1.0*0*Vars.v**3 - 1.0*4*Cons.E_2*Cons.M*Vars.v**4 - 1.0*5*Vars.E_SO_3*Cons.M*Vars.v**5) + -Cons.Fcal_3 - Vars.Fcal_SO_3*(-1.0*2*Cons.E_0*Cons.M*Vars.v**3 - 1.0*0*Vars.v**4 - 1.0*4*Cons.E_2*Cons.M*Vars.v**5) + -Cons.Fcal_4 - Vars.Fcal_SQ_4*(-1.0*2*Cons.E_0*Cons.M*Vars.v**4 - 1.0*0*Vars.v**5) - 1.0*-Cons.Fcal_5 - Vars.Fcal_SO_5*2*Cons.E_0*Cons.M*Vars.v**5 + (0*(0*(1.0*2*Cons.E_0*Cons.M*Vars.v**2 + 1.0*0*Vars.v**3 + 1.0*4*Cons.E_2*Cons.M*Vars.v**4 + 1.0*5*Vars.E_SO_3*Cons.M*Vars.v**5) + -Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*(2.0*2*Cons.E_0*Cons.M*Vars.v**3 + 2.0*0*Vars.v**4 + 2.0*4*Cons.E_2*Cons.M*Vars.v**5) + -Cons.Fcal_3 - Vars.Fcal_SO_3*(2.0*2*Cons.E_0*Cons.M*Vars.v**4 + 2.0*0*Vars.v**5) + 2.0*-Cons.Fcal_4 - Vars.Fcal_SQ_4*2*Cons.E_0*Cons.M*Vars.v**5) + -Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*(-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*(1.0*2*Cons.E_0*Cons.M*Vars.v**4 + 1.0*0*Vars.v**5) + 2.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*2*Cons.E_0*Cons.M*Vars.v**5) + (0*(0*(0*(-1.0*2*Cons.E_0*Cons.M*Vars.v**3 - 1.0*0*Vars.v**4 - 1.0*4*Cons.E_2*Cons.M*Vars.v**5) + -Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*(-3.0*2*Cons.E_0*Cons.M*Vars.v**4 - 3.0*0*Vars.v**5) - 3.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*2*Cons.E_0*Cons.M*Vars.v**5) - 3.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff**2*2*Cons.E_0*Cons.M*Vars.v**5) + (0**3*(0*(1.0*2*Cons.E_0*Cons.M*Vars.v**4 + 1.0*0*Vars.v**5) + 4.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*2*Cons.E_0*Cons.M*Vars.v**5) - 1.0*0**5*2*Cons.E_0*Cons.M*Vars.v**5/-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff)/-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff)/-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff)/-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff)/-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff)/(Vars.Fcal_coeff*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff)
    dvdt_T5 = 1.0/dtdv
    dydt=np.zeros(10)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_2p5(Cons,Vars)[1:])
    dydt[0] = dvdt_T5[0]
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_2p5(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_2p5(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    dydt[9] = np.linalg.norm(OmegaVec_2p5(Cons,Vars)[1:])
    return dydt      

@njit(cache=True)
def Recalculate_3p0(Cons,Vars,y):
    Vars.v = np.array([y[0]])
    Vars.rfrak_chi1 = np.array([y[1],y[2]])
    Vars.rfrak_chi2 = np.array([y[3],y[4]])
    Vars.Ecc = Cons.e_i*(Cons.v_i/Vars.v)**(19/6)
    Vars.xi = y[9]+Cons.xi_i
    Vars.R = y[5]*Cons.wHat + y[6]*Cons.xHat + y[7]*Cons.yHat + y[8]*Cons.zHat
    Vars.nHat = mul(mul(Vars.R,Cons.xHat),inverse(Vars.R))
    Vars.lambdaHat = mul(mul(Vars.R,Cons.yHat),inverse(Vars.R))
    Vars.ellHat = mul(mul(Vars.R,Cons.zHat),inverse(Vars.R))
    Vars.R_S1 = exp(Vars.rfrak_chi1[0]*Cons.xHat + Vars.rfrak_chi1[1]*Cons.yHat)
    Vars.R_S2 = exp(Vars.rfrak_chi2[0]*Cons.xHat + Vars.rfrak_chi2[1]*Cons.yHat)
    Vars.chiVec1 = mul(mul(mul(Cons.S_chi1,Vars.R_S1),Cons.zHat),mul(conjugate(Vars.R_S1),conjugate(Cons.S_chi1)))
    Vars.chiVec2 = mul(mul(mul(Cons.S_chi2,Vars.R_S2),Cons.zHat),mul(conjugate(Vars.R_S2),conjugate(Cons.S_chi2)))
    Vars.chi1_n = np.array([dot(Vars.chiVec1[1:],Vars.nHat[1:])])
    Vars.chi1_lambda = np.array([dot(Vars.chiVec1[1:],Vars.lambdaHat[1:])])
    Vars.chi1_ell = np.array([dot(Vars.chiVec1[1:],Vars.ellHat[1:])])
    Vars.chi2_n = np.array([dot(Vars.chiVec2[1:],Vars.nHat[1:])])
    Vars.chi2_lambda = np.array([dot(Vars.chiVec2[1:],Vars.lambdaHat[1:])])
    Vars.chi2_ell = np.array([dot(Vars.chiVec2[1:],Vars.ellHat[1:])])
    Vars.S_ell = Cons.M1**2*Vars.chi1_ell + Cons.M2**2*Vars.chi2_ell
    Vars.S_n = Cons.M1**2*Vars.chi1_n + Cons.M2**2*Vars.chi2_n
    Vars.S_lambda = Cons.M1**2*Vars.chi1_lambda + Cons.M2**2*Vars.chi2_lambda
    Vars.Sigma_ell = Cons.M*(-Cons.M1*Vars.chi1_ell + Cons.M2*Vars.chi2_ell)
    Vars.Sigma_n = Cons.M*(-Cons.M1*Vars.chi1_n + Cons.M2*Vars.chi2_n)
    Vars.Sigma_lambda = Cons.M*(-Cons.M1*Vars.chi1_lambda + Cons.M2*Vars.chi2_lambda)
    Vars.chi_s_ell = Vars.chi1_ell/2 + Vars.chi2_ell/2
    Vars.chi_a_ell = Vars.chi1_ell/2 - Vars.chi2_ell/2
    Vars.logv = log(Vars.v)
    Vars.Fcal_coeff = 32*Cons.nu**2*Vars.v**10/5
    Vars.Fcal_SQ_4 = Cons.chi1chi1*(-89*Cons.delta/192 + 89*Cons.nu/96 - 89/192) - 103*Cons.chi1chi2*Cons.nu/48 + Cons.chi2chi2*(89*Cons.delta/192 + 89*Cons.nu/96 - 89/192) + Vars.chi_a_ell*(Vars.chi_a_ell*(287/96 - 12*Cons.nu) + 287*Vars.chi_s_ell*Cons.delta/48) + Vars.chi_s_ell**2*(Cons.nu/24 + 287/96)
    Vars.Fcal_SO_3 = (-4*Vars.S_ell - 5*Vars.Sigma_ell*Cons.delta/4)/Cons.M**2
    Vars.Fcal_SO_5 = (Vars.S_ell*(272*Cons.nu/9 - 9/2) + Vars.Sigma_ell*Cons.delta*(43*Cons.nu/4 - 13/16))/Cons.M**2
    Vars.Fcal_SO_6 = (-16*Vars.S_ell*pi - 31*Vars.Sigma_ell*Cons.delta*pi/6)/Cons.M**2
    Vars.Fecc_coeff = 157*Cons.e_i**2.0*Cons.v_i**6.33333333333333/(24*Vars.v**6.33333333333333)
    Vars.E_SQ_4 = -3*Vars.chi_a_ell**2/2 - 3*Vars.chi_s_ell**2/2 - Cons.delta*(Cons.chi2chi2/2 + 3*Vars.chi_a_ell*Vars.chi_s_ell) + Cons.nu*(Cons.chi1chi2 + 6*Vars.chi_a_ell**2) + (Cons.chi1chi1 + Cons.chi2chi2)*(Cons.delta - 2*Cons.nu + 1)/4
    Vars.E_SO_3 = (14*Vars.S_ell/3 + 2*Vars.Sigma_ell*Cons.delta)/Cons.M**2
    Vars.E_SO_5 = (Vars.S_ell*(11 - 61*Cons.nu/9) + Vars.Sigma_ell*Cons.delta*(3 - 10*Cons.nu/3))/Cons.M**2

@njit
def OmegaVec_chiVec_1_3p0(Cons,Vars):
    Omega1_coeff = Vars.v**5/Cons.M
    return Omega1_coeff*(Vars.ellHat*(-0.75*Cons.delta + 0.5*Cons.nu + Vars.v**2*(Cons.delta*(0.625*Cons.nu - 0.5625) + Cons.nu*(1.25 - 0.0416666666666667*Cons.nu) + Vars.v**2*(Cons.delta*(Cons.nu*(4.875 - 0.15625*Cons.nu) - 0.84375) + Cons.nu*(Cons.nu*(-0.0208333333333333*Cons.nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + Vars.nHat*Vars.v*(3.0*Vars.chi1_n*Cons.nu + 3.0*Cons.M2**2*Vars.chi2_n/Cons.M**2) - Cons.M2**2*Vars.chiVec2*Vars.v/Cons.M**2)

@njit
def OmegaVec_chiVec_2_3p0(Cons,Vars):
    Omega2_coeff = Vars.v**5/Cons.M
    return Omega2_coeff*(Vars.ellHat*(0.75*Cons.delta + 0.5*Cons.nu + Vars.v**2*(Cons.delta*(0.5625 - 0.625*Cons.nu) + Cons.nu*(1.25 - 0.0416666666666667*Cons.nu) + Vars.v**2*(Cons.delta*(Cons.nu*(0.15625*Cons.nu - 4.875) + 0.84375) + Cons.nu*(Cons.nu*(-0.0208333333333333*Cons.nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + Vars.nHat*Vars.v*(3.0*Vars.chi2_n*Cons.nu + 3.0*Cons.M1**2*Vars.chi1_n/Cons.M**2) - Cons.M1**2*Vars.chiVec1*Vars.v/Cons.M**2)

@njit
def OmegaVec_3p0(Cons,Vars):
    gamma_PN_0 = 1.00000000000000
    a_ell_0 = 7.0*Vars.S_n + 3.0*Vars.Sigma_n*Cons.delta
    gamma_PN_3 = (1.66666666666667*Vars.S_ell + Vars.Sigma_ell*Cons.delta)/Cons.M**2
    a_ell_4 = Vars.S_n*(5.77777777777778*Cons.nu**2 + 14.75*Cons.nu + 1.5) + Vars.Sigma_n*Cons.delta*(2.83333333333333*Cons.nu**2 + 9.125*Cons.nu + 1.5)
    a_ell_2 = Vars.S_n*(-9.66666666666667*Cons.nu - 10.0) + Vars.Sigma_n*Cons.delta*(-4.5*Cons.nu - 6.0)
    gamma_PN_5 = (Vars.S_ell*(0.888888888888889*Cons.nu + 3.33333333333333) + 2.0*Vars.Sigma_ell*Cons.delta)/Cons.M**2
    gamma_PN_4 = 1.0 - 5.41666666666667*Cons.nu
    gamma_PN_6 = 0.0123456790123457*Cons.nu**3 + 6.36111111111111*Cons.nu**2 - 2.98177812235564*Cons.nu + 1.0
    gamma_PN_2 = 1.0 - 0.333333333333333*Cons.nu
    return Vars.ellHat*Vars.v**3*np.sqrt(1.0 - Vars.Ecc**2.0)*(-Vars.Ecc*np.cos(Vars.xi) + 1.0)**(-2.0)/Cons.M + Vars.nHat*Vars.v**6*(a_ell_0 + Vars.v**2*(a_ell_2 + a_ell_4*Vars.v**2))*(gamma_PN_0 + Vars.v**2*(gamma_PN_2 + Vars.v*(gamma_PN_3 + Vars.v*(gamma_PN_4 + Vars.v*(gamma_PN_5 + gamma_PN_6*Vars.v)))))/Cons.M**3


@njit(cache=True)
def TaylorT1_3p0(Cons,Vars):
    Flux = Vars.Fcal_coeff*Vars.Fecc_coeff*(Cons.Fecc_0 + Cons.Fecc_2*Vars.v**2) + Vars.Fcal_coeff*(Cons.Fcal_0 + Vars.v**2*(Cons.Fcal_2 + Vars.v*(Cons.Fcal_3 + Vars.Fcal_SO_3 + Vars.v*(Cons.Fcal_4 + Vars.Fcal_SQ_4 + Vars.v*(Cons.Fcal_5 + Vars.Fcal_SO_5 + Vars.v*(Cons.Fcal_6 + Vars.Fcal_SO_6 + Cons.Fcal_lnv_6*Vars.logv))))))
    dEdV = 1.16666666666667*Cons.Eecc_2*Cons.M*Cons.e_i**2.0*Cons.nu*Vars.v**(-3.33333333333333)*Cons.v_i**6.33333333333333 - 0.5*Cons.M*Cons.nu*Vars.v*(2.0*Cons.E_0 + Vars.v**2*(4.0*Cons.E_2 + Vars.v*(5.0*Vars.E_SO_3 + Vars.v*(6.0*Cons.E_4 + 6.0*Vars.E_SQ_4 + Vars.v*(8.0*Cons.E_6*Vars.v + 7.0*Vars.E_SO_5)))))
    Absorption = 0
    dvdt_T1 = (-Absorption - Flux)/dEdV
    dydt=np.zeros(10)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_3p0(Cons,Vars)[1:])
    dydt[0] = dvdt_T1[0] 
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_3p0(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_3p0(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    dydt[9] = np.linalg.norm(OmegaVec_3p0(Cons,Vars)[1:])
    return dydt

@njit(cache=True)
def TaylorT4_3p0(Cons,Vars):
    dvdt_T4 = -2.0*Vars.Fcal_coeff*(1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff + 1.0*0*Vars.v + 1.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**2 + 1.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*Vars.v**3 + 1.0*-Cons.Fcal_4 - Vars.Fcal_SQ_4*Vars.v**4 + 1.0*-Cons.Fcal_5 - Vars.Fcal_SO_5*Vars.v**5 + 1.0*-Cons.Fcal_6 - Vars.Fcal_SO_6 - Cons.Fcal_lnv_6*Vars.logv*Vars.v**6 + (0*(-1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v - 1.0*0*Vars.v**2 - 1.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**3 - 1.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*Vars.v**4 - 1.0*-Cons.Fcal_4 - Vars.Fcal_SQ_4*Vars.v**5 - 1.0*-Cons.Fcal_5 - Vars.Fcal_SO_5*Vars.v**6) + 4*Cons.E_2*Cons.M*(-1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**2 - 1.0*0*Vars.v**3 - 1.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**4 - 1.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*Vars.v**5 - 1.0*-Cons.Fcal_4 - Vars.Fcal_SQ_4*Vars.v**6) + 5*Vars.E_SO_3*Cons.M*(-1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**3 - 1.0*0*Vars.v**4 - 1.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**5 - 1.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*Vars.v**6) + 6*Cons.E_4*Cons.M + 6*Vars.E_SQ_4*Cons.M*(-1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**4 - 1.0*0*Vars.v**5 - 1.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**6) + 7*Vars.E_SO_5*Cons.M*(-1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**5 - 1.0*0*Vars.v**6) - 1.0*8*Cons.E_6*Cons.M*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**6 + (0*(0*(1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**2 + 1.0*0*Vars.v**3 + 1.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**4 + 1.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*Vars.v**5 + 1.0*-Cons.Fcal_4 - Vars.Fcal_SQ_4*Vars.v**6) + 4*Cons.E_2*Cons.M*(2.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**3 + 2.0*0*Vars.v**4 + 2.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**5 + 2.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*Vars.v**6) + 5*Vars.E_SO_3*Cons.M*(2.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**4 + 2.0*0*Vars.v**5 + 2.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**6) + 6*Cons.E_4*Cons.M + 6*Vars.E_SQ_4*Cons.M*(2.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**5 + 2.0*0*Vars.v**6) + 2.0*7*Vars.E_SO_5*Cons.M*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**6) + 4*Cons.E_2*Cons.M*(4*Cons.E_2*Cons.M*(1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**4 + 1.0*0*Vars.v**5 + 1.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**6) + 5*Vars.E_SO_3*Cons.M*(2.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**5 + 2.0*0*Vars.v**6) + 2.0*6*Cons.E_4*Cons.M + 6*Vars.E_SQ_4*Cons.M*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**6) + 1.0*5*Vars.E_SO_3*Cons.M**2*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**6 + (0*(0*(0*(-1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**3 - 1.0*0*Vars.v**4 - 1.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**5 - 1.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*Vars.v**6) + 4*Cons.E_2*Cons.M*(-3.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**4 - 3.0*0*Vars.v**5 - 3.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**6) + 5*Vars.E_SO_3*Cons.M*(-3.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**5 - 3.0*0*Vars.v**6) - 3.0*6*Cons.E_4*Cons.M + 6*Vars.E_SQ_4*Cons.M*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**6) + 4*Cons.E_2*Cons.M*(4*Cons.E_2*Cons.M*(-3.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**5 - 3.0*0*Vars.v**6) - 6.0*5*Vars.E_SO_3*Cons.M*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**6)) - 1.0*4*Cons.E_2*Cons.M**3*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**6 + (0**2*(0*(0*(1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**4 + 1.0*0*Vars.v**5 + 1.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**6) + 4*Cons.E_2*Cons.M*(4.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**5 + 4.0*0*Vars.v**6) + 4.0*5*Vars.E_SO_3*Cons.M*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**6) + 6.0*4*Cons.E_2*Cons.M**2*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**6) + (0**4*(0*(-1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**5 - 1.0*0*Vars.v**6) - 5.0*4*Cons.E_2*Cons.M*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**6) + 1.0*0**6*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**6/2*Cons.E_0*Cons.M)/2*Cons.E_0*Cons.M)/2*Cons.E_0*Cons.M)/2*Cons.E_0*Cons.M)/2*Cons.E_0*Cons.M)/2*Cons.E_0*Cons.M)/(Cons.nu*Vars.v*2*Cons.E_0*Cons.M)
    dydt=np.zeros(10)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_3p0(Cons,Vars)[1:])
    dydt[0] = dvdt_T4[0]
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_3p0(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_3p0(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    dydt[9] = np.linalg.norm(OmegaVec_3p0(Cons,Vars)[1:])
    return dydt      

@njit(cache=True)
def TaylorT5_3p0(Cons,Vars):
    dtdv = -0.5*Cons.nu*Vars.v*(1.0*2*Cons.E_0*Cons.M + 1.0*0*Vars.v + 1.0*4*Cons.E_2*Cons.M*Vars.v**2 + 1.0*5*Vars.E_SO_3*Cons.M*Vars.v**3 + 1.0*6*Cons.E_4*Cons.M + 6*Vars.E_SQ_4*Cons.M*Vars.v**4 + 1.0*7*Vars.E_SO_5*Cons.M*Vars.v**5 + 1.0*8*Cons.E_6*Cons.M*Vars.v**6 + (0*(-1.0*2*Cons.E_0*Cons.M*Vars.v - 1.0*0*Vars.v**2 - 1.0*4*Cons.E_2*Cons.M*Vars.v**3 - 1.0*5*Vars.E_SO_3*Cons.M*Vars.v**4 - 1.0*6*Cons.E_4*Cons.M + 6*Vars.E_SQ_4*Cons.M*Vars.v**5 - 1.0*7*Vars.E_SO_5*Cons.M*Vars.v**6) + -Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*(-1.0*2*Cons.E_0*Cons.M*Vars.v**2 - 1.0*0*Vars.v**3 - 1.0*4*Cons.E_2*Cons.M*Vars.v**4 - 1.0*5*Vars.E_SO_3*Cons.M*Vars.v**5 - 1.0*6*Cons.E_4*Cons.M + 6*Vars.E_SQ_4*Cons.M*Vars.v**6) + -Cons.Fcal_3 - Vars.Fcal_SO_3*(-1.0*2*Cons.E_0*Cons.M*Vars.v**3 - 1.0*0*Vars.v**4 - 1.0*4*Cons.E_2*Cons.M*Vars.v**5 - 1.0*5*Vars.E_SO_3*Cons.M*Vars.v**6) + -Cons.Fcal_4 - Vars.Fcal_SQ_4*(-1.0*2*Cons.E_0*Cons.M*Vars.v**4 - 1.0*0*Vars.v**5 - 1.0*4*Cons.E_2*Cons.M*Vars.v**6) + -Cons.Fcal_5 - Vars.Fcal_SO_5*(-1.0*2*Cons.E_0*Cons.M*Vars.v**5 - 1.0*0*Vars.v**6) - 1.0*-Cons.Fcal_6 - Vars.Fcal_SO_6 - Cons.Fcal_lnv_6*Vars.logv*2*Cons.E_0*Cons.M*Vars.v**6 + (0*(0*(1.0*2*Cons.E_0*Cons.M*Vars.v**2 + 1.0*0*Vars.v**3 + 1.0*4*Cons.E_2*Cons.M*Vars.v**4 + 1.0*5*Vars.E_SO_3*Cons.M*Vars.v**5 + 1.0*6*Cons.E_4*Cons.M + 6*Vars.E_SQ_4*Cons.M*Vars.v**6) + -Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*(2.0*2*Cons.E_0*Cons.M*Vars.v**3 + 2.0*0*Vars.v**4 + 2.0*4*Cons.E_2*Cons.M*Vars.v**5 + 2.0*5*Vars.E_SO_3*Cons.M*Vars.v**6) + -Cons.Fcal_3 - Vars.Fcal_SO_3*(2.0*2*Cons.E_0*Cons.M*Vars.v**4 + 2.0*0*Vars.v**5 + 2.0*4*Cons.E_2*Cons.M*Vars.v**6) + -Cons.Fcal_4 - Vars.Fcal_SQ_4*(2.0*2*Cons.E_0*Cons.M*Vars.v**5 + 2.0*0*Vars.v**6) + 2.0*-Cons.Fcal_5 - Vars.Fcal_SO_5*2*Cons.E_0*Cons.M*Vars.v**6) + -Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*(-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*(1.0*2*Cons.E_0*Cons.M*Vars.v**4 + 1.0*0*Vars.v**5 + 1.0*4*Cons.E_2*Cons.M*Vars.v**6) + -Cons.Fcal_3 - Vars.Fcal_SO_3*(2.0*2*Cons.E_0*Cons.M*Vars.v**5 + 2.0*0*Vars.v**6) + 2.0*-Cons.Fcal_4 - Vars.Fcal_SQ_4*2*Cons.E_0*Cons.M*Vars.v**6) + 1.0*-Cons.Fcal_3 - Vars.Fcal_SO_3**2*2*Cons.E_0*Cons.M*Vars.v**6 + (0*(0*(0*(-1.0*2*Cons.E_0*Cons.M*Vars.v**3 - 1.0*0*Vars.v**4 - 1.0*4*Cons.E_2*Cons.M*Vars.v**5 - 1.0*5*Vars.E_SO_3*Cons.M*Vars.v**6) + -Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*(-3.0*2*Cons.E_0*Cons.M*Vars.v**4 - 3.0*0*Vars.v**5 - 3.0*4*Cons.E_2*Cons.M*Vars.v**6) + -Cons.Fcal_3 - Vars.Fcal_SO_3*(-3.0*2*Cons.E_0*Cons.M*Vars.v**5 - 3.0*0*Vars.v**6) - 3.0*-Cons.Fcal_4 - Vars.Fcal_SQ_4*2*Cons.E_0*Cons.M*Vars.v**6) + -Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*(-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*(-3.0*2*Cons.E_0*Cons.M*Vars.v**5 - 3.0*0*Vars.v**6) - 6.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*2*Cons.E_0*Cons.M*Vars.v**6)) - 1.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff**3*2*Cons.E_0*Cons.M*Vars.v**6 + (0**2*(0*(0*(1.0*2*Cons.E_0*Cons.M*Vars.v**4 + 1.0*0*Vars.v**5 + 1.0*4*Cons.E_2*Cons.M*Vars.v**6) + -Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*(4.0*2*Cons.E_0*Cons.M*Vars.v**5 + 4.0*0*Vars.v**6) + 4.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*2*Cons.E_0*Cons.M*Vars.v**6) + 6.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff**2*2*Cons.E_0*Cons.M*Vars.v**6) + (0**4*(0*(-1.0*2*Cons.E_0*Cons.M*Vars.v**5 - 1.0*0*Vars.v**6) - 5.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*2*Cons.E_0*Cons.M*Vars.v**6) + 1.0*0**6*2*Cons.E_0*Cons.M*Vars.v**6/-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff)/-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff)/-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff)/-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff)/-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff)/-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff)/(Vars.Fcal_coeff*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff)
    dvdt_T5 = 1.0/dtdv
    dydt=np.zeros(10)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_3p0(Cons,Vars)[1:])
    dydt[0] = dvdt_T5[0]
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_3p0(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_3p0(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    dydt[9] = np.linalg.norm(OmegaVec_3p0(Cons,Vars)[1:])
    return dydt      

@njit(cache=True)
def Recalculate_3p5(Cons,Vars,y):
    Vars.v = np.array([y[0]])
    Vars.rfrak_chi1 = np.array([y[1],y[2]])
    Vars.rfrak_chi2 = np.array([y[3],y[4]])
    Vars.Ecc = Cons.e_i*(Cons.v_i/Vars.v)**(19/6)
    Vars.xi = y[9]+Cons.xi_i
    Vars.R = y[5]*Cons.wHat + y[6]*Cons.xHat + y[7]*Cons.yHat + y[8]*Cons.zHat
    Vars.nHat = mul(mul(Vars.R,Cons.xHat),inverse(Vars.R))
    Vars.lambdaHat = mul(mul(Vars.R,Cons.yHat),inverse(Vars.R))
    Vars.ellHat = mul(mul(Vars.R,Cons.zHat),inverse(Vars.R))
    Vars.R_S1 = exp(Vars.rfrak_chi1[0]*Cons.xHat + Vars.rfrak_chi1[1]*Cons.yHat)
    Vars.R_S2 = exp(Vars.rfrak_chi2[0]*Cons.xHat + Vars.rfrak_chi2[1]*Cons.yHat)
    Vars.chiVec1 = mul(mul(mul(Cons.S_chi1,Vars.R_S1),Cons.zHat),mul(conjugate(Vars.R_S1),conjugate(Cons.S_chi1)))
    Vars.chiVec2 = mul(mul(mul(Cons.S_chi2,Vars.R_S2),Cons.zHat),mul(conjugate(Vars.R_S2),conjugate(Cons.S_chi2)))
    Vars.chi1_n = np.array([dot(Vars.chiVec1[1:],Vars.nHat[1:])])
    Vars.chi1_lambda = np.array([dot(Vars.chiVec1[1:],Vars.lambdaHat[1:])])
    Vars.chi1_ell = np.array([dot(Vars.chiVec1[1:],Vars.ellHat[1:])])
    Vars.chi2_n = np.array([dot(Vars.chiVec2[1:],Vars.nHat[1:])])
    Vars.chi2_lambda = np.array([dot(Vars.chiVec2[1:],Vars.lambdaHat[1:])])
    Vars.chi2_ell = np.array([dot(Vars.chiVec2[1:],Vars.ellHat[1:])])
    Vars.S_ell = Cons.M1**2*Vars.chi1_ell + Cons.M2**2*Vars.chi2_ell
    Vars.S_n = Cons.M1**2*Vars.chi1_n + Cons.M2**2*Vars.chi2_n
    Vars.S_lambda = Cons.M1**2*Vars.chi1_lambda + Cons.M2**2*Vars.chi2_lambda
    Vars.Sigma_ell = Cons.M*(-Cons.M1*Vars.chi1_ell + Cons.M2*Vars.chi2_ell)
    Vars.Sigma_n = Cons.M*(-Cons.M1*Vars.chi1_n + Cons.M2*Vars.chi2_n)
    Vars.Sigma_lambda = Cons.M*(-Cons.M1*Vars.chi1_lambda + Cons.M2*Vars.chi2_lambda)
    Vars.chi_s_ell = Vars.chi1_ell/2 + Vars.chi2_ell/2
    Vars.chi_a_ell = Vars.chi1_ell/2 - Vars.chi2_ell/2
    Vars.logv = log(Vars.v)
    Vars.Fcal_coeff = 32*Cons.nu**2*Vars.v**10/5
    Vars.Fcal_SQ_4 = Cons.chi1chi1*(-89*Cons.delta/192 + 89*Cons.nu/96 - 89/192) - 103*Cons.chi1chi2*Cons.nu/48 + Cons.chi2chi2*(89*Cons.delta/192 + 89*Cons.nu/96 - 89/192) + Vars.chi_a_ell*(Vars.chi_a_ell*(287/96 - 12*Cons.nu) + 287*Vars.chi_s_ell*Cons.delta/48) + Vars.chi_s_ell**2*(Cons.nu/24 + 287/96)
    Vars.Fcal_SO_3 = (-4*Vars.S_ell - 5*Vars.Sigma_ell*Cons.delta/4)/Cons.M**2
    Vars.Fcal_SO_5 = (Vars.S_ell*(272*Cons.nu/9 - 9/2) + Vars.Sigma_ell*Cons.delta*(43*Cons.nu/4 - 13/16))/Cons.M**2
    Vars.Fcal_SO_6 = (-16*Vars.S_ell*pi - 31*Vars.Sigma_ell*Cons.delta*pi/6)/Cons.M**2
    Vars.Fcal_SO_7 = (Vars.S_ell*(-2810*Cons.nu**2/27 + 6172*Cons.nu/189 + 476645/6804) + Vars.Sigma_ell*Cons.delta*(-1501*Cons.nu**2/36 + 1849*Cons.nu/126 + 9535/336))/Cons.M**2
    Vars.Fecc_coeff = 157*Cons.e_i**2.0*Cons.v_i**6.33333333333333/(24*Vars.v**6.33333333333333)
    Vars.E_SQ_4 = -3*Vars.chi_a_ell**2/2 - 3*Vars.chi_s_ell**2/2 - Cons.delta*(Cons.chi2chi2/2 + 3*Vars.chi_a_ell*Vars.chi_s_ell) + Cons.nu*(Cons.chi1chi2 + 6*Vars.chi_a_ell**2) + (Cons.chi1chi1 + Cons.chi2chi2)*(Cons.delta - 2*Cons.nu + 1)/4
    Vars.E_SO_3 = (14*Vars.S_ell/3 + 2*Vars.Sigma_ell*Cons.delta)/Cons.M**2
    Vars.E_SO_5 = (Vars.S_ell*(11 - 61*Cons.nu/9) + Vars.Sigma_ell*Cons.delta*(3 - 10*Cons.nu/3))/Cons.M**2
    Vars.E_SO_7 = (Vars.S_ell*(29*Cons.nu**2/12 - 367*Cons.nu/4 + 135/4) + Vars.Sigma_ell*Cons.delta*(5*Cons.nu**2/4 - 39*Cons.nu + 27/4))/Cons.M**2

@njit
def OmegaVec_chiVec_1_3p5(Cons,Vars):
    Omega1_coeff = Vars.v**5/Cons.M
    return Omega1_coeff*(Vars.ellHat*(-0.75*Cons.delta + 0.5*Cons.nu + Vars.v**2*(Cons.delta*(0.625*Cons.nu - 0.5625) + Cons.nu*(1.25 - 0.0416666666666667*Cons.nu) + Vars.v**2*(Cons.delta*(Cons.nu*(4.875 - 0.15625*Cons.nu) - 0.84375) + Cons.nu*(Cons.nu*(-0.0208333333333333*Cons.nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + Vars.nHat*Vars.v*(3.0*Vars.chi1_n*Cons.nu + 3.0*Cons.M2**2*Vars.chi2_n/Cons.M**2) - Cons.M2**2*Vars.chiVec2*Vars.v/Cons.M**2)

@njit
def OmegaVec_chiVec_2_3p5(Cons,Vars):
    Omega2_coeff = Vars.v**5/Cons.M
    return Omega2_coeff*(Vars.ellHat*(0.75*Cons.delta + 0.5*Cons.nu + Vars.v**2*(Cons.delta*(0.5625 - 0.625*Cons.nu) + Cons.nu*(1.25 - 0.0416666666666667*Cons.nu) + Vars.v**2*(Cons.delta*(Cons.nu*(0.15625*Cons.nu - 4.875) + 0.84375) + Cons.nu*(Cons.nu*(-0.0208333333333333*Cons.nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + Vars.nHat*Vars.v*(3.0*Vars.chi2_n*Cons.nu + 3.0*Cons.M1**2*Vars.chi1_n/Cons.M**2) - Cons.M1**2*Vars.chiVec1*Vars.v/Cons.M**2)

@njit
def OmegaVec_3p5(Cons,Vars):
    gamma_PN_0 = 1.00000000000000
    a_ell_0 = 7.0*Vars.S_n + 3.0*Vars.Sigma_n*Cons.delta
    gamma_PN_3 = (1.66666666666667*Vars.S_ell + Vars.Sigma_ell*Cons.delta)/Cons.M**2
    a_ell_4 = Vars.S_n*(5.77777777777778*Cons.nu**2 + 14.75*Cons.nu + 1.5) + Vars.Sigma_n*Cons.delta*(2.83333333333333*Cons.nu**2 + 9.125*Cons.nu + 1.5)
    a_ell_2 = Vars.S_n*(-9.66666666666667*Cons.nu - 10.0) + Vars.Sigma_n*Cons.delta*(-4.5*Cons.nu - 6.0)
    gamma_PN_5 = (Vars.S_ell*(0.888888888888889*Cons.nu + 3.33333333333333) + 2.0*Vars.Sigma_ell*Cons.delta)/Cons.M**2
    gamma_PN_7 = (Vars.S_ell*(-6.0*Cons.nu**2 - 10.5833333333333*Cons.nu + 5.0) - 2.66666666666667*Vars.Sigma_ell*Cons.delta*Cons.nu**2 + Vars.Sigma_ell*Cons.delta*(3.0 - 10.1666666666667*Cons.nu))/Cons.M**2
    gamma_PN_4 = 1.0 - 5.41666666666667*Cons.nu
    gamma_PN_6 = 0.0123456790123457*Cons.nu**3 + 6.36111111111111*Cons.nu**2 - 2.98177812235564*Cons.nu + 1.0
    gamma_PN_2 = 1.0 - 0.333333333333333*Cons.nu
    return Vars.ellHat*Vars.v**3*np.sqrt(1.0 - Vars.Ecc**2.0)*(-Vars.Ecc*np.cos(Vars.xi) + 1.0)**(-2.0)/Cons.M + Vars.nHat*Vars.v**6*(a_ell_0 + Vars.v**2*(a_ell_2 + a_ell_4*Vars.v**2))*(gamma_PN_0 + Vars.v**2*(gamma_PN_2 + Vars.v*(gamma_PN_3 + Vars.v*(gamma_PN_4 + Vars.v*(gamma_PN_5 + Vars.v*(gamma_PN_6 + gamma_PN_7*Vars.v))))))/Cons.M**3


@njit(cache=True)
def TaylorT1_3p5(Cons,Vars):
    Flux = Vars.Fcal_coeff*Vars.Fecc_coeff*(Cons.Fecc_0 + Cons.Fecc_2*Vars.v**2) + Vars.Fcal_coeff*(Cons.Fcal_0 + Vars.v**2*(Cons.Fcal_2 + Vars.v*(Cons.Fcal_3 + Vars.Fcal_SO_3 + Vars.v*(Cons.Fcal_4 + Vars.Fcal_SQ_4 + Vars.v*(Cons.Fcal_5 + Vars.Fcal_SO_5 + Vars.v*(Cons.Fcal_6 + Vars.Fcal_SO_6 + Cons.Fcal_lnv_6*Vars.logv + Vars.v*(Cons.Fcal_7 + Vars.Fcal_SO_7)))))))
    dEdV = 1.16666666666667*Cons.Eecc_2*Cons.M*Cons.e_i**2.0*Cons.nu*Vars.v**(-3.33333333333333)*Cons.v_i**6.33333333333333 - 0.5*Cons.M*Cons.nu*Vars.v*(2.0*Cons.E_0 + Vars.v**2*(4.0*Cons.E_2 + Vars.v*(5.0*Vars.E_SO_3 + Vars.v*(6.0*Cons.E_4 + 6.0*Vars.E_SQ_4 + Vars.v*(7.0*Vars.E_SO_5 + Vars.v*(8.0*Cons.E_6 + 9.0*Vars.E_SO_7*Vars.v))))))
    Absorption = 0
    dvdt_T1 = (-Absorption - Flux)/dEdV
    dydt=np.zeros(10)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_3p5(Cons,Vars)[1:])
    dydt[0] = dvdt_T1[0] 
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_3p5(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_3p5(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    dydt[9] = np.linalg.norm(OmegaVec_3p5(Cons,Vars)[1:])
    return dydt

@njit(cache=True)
def TaylorT4_3p5(Cons,Vars):
    dvdt_T4 = -2.0*Vars.Fcal_coeff*(1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff + 1.0*0*Vars.v + 1.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**2 + 1.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*Vars.v**3 + 1.0*-Cons.Fcal_4 - Vars.Fcal_SQ_4*Vars.v**4 + 1.0*-Cons.Fcal_5 - Vars.Fcal_SO_5*Vars.v**5 + 1.0*-Cons.Fcal_6 - Vars.Fcal_SO_6 - Cons.Fcal_lnv_6*Vars.logv*Vars.v**6 + 1.0*-Cons.Fcal_7 - Vars.Fcal_SO_7*Vars.v**7 + (0*(-1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v - 1.0*0*Vars.v**2 - 1.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**3 - 1.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*Vars.v**4 - 1.0*-Cons.Fcal_4 - Vars.Fcal_SQ_4*Vars.v**5 - 1.0*-Cons.Fcal_5 - Vars.Fcal_SO_5*Vars.v**6 - 1.0*-Cons.Fcal_6 - Vars.Fcal_SO_6 - Cons.Fcal_lnv_6*Vars.logv*Vars.v**7) + 4*Cons.E_2*Cons.M*(-1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**2 - 1.0*0*Vars.v**3 - 1.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**4 - 1.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*Vars.v**5 - 1.0*-Cons.Fcal_4 - Vars.Fcal_SQ_4*Vars.v**6 - 1.0*-Cons.Fcal_5 - Vars.Fcal_SO_5*Vars.v**7) + 5*Vars.E_SO_3*Cons.M*(-1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**3 - 1.0*0*Vars.v**4 - 1.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**5 - 1.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*Vars.v**6 - 1.0*-Cons.Fcal_4 - Vars.Fcal_SQ_4*Vars.v**7) + 6*Cons.E_4*Cons.M + 6*Vars.E_SQ_4*Cons.M*(-1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**4 - 1.0*0*Vars.v**5 - 1.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**6 - 1.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*Vars.v**7) + 7*Vars.E_SO_5*Cons.M*(-1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**5 - 1.0*0*Vars.v**6 - 1.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**7) + 8*Cons.E_6*Cons.M*(-1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**6 - 1.0*0*Vars.v**7) - 1.0*9*Vars.E_SO_7*Cons.M*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**7 + (0*(0*(1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**2 + 1.0*0*Vars.v**3 + 1.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**4 + 1.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*Vars.v**5 + 1.0*-Cons.Fcal_4 - Vars.Fcal_SQ_4*Vars.v**6 + 1.0*-Cons.Fcal_5 - Vars.Fcal_SO_5*Vars.v**7) + 4*Cons.E_2*Cons.M*(2.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**3 + 2.0*0*Vars.v**4 + 2.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**5 + 2.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*Vars.v**6 + 2.0*-Cons.Fcal_4 - Vars.Fcal_SQ_4*Vars.v**7) + 5*Vars.E_SO_3*Cons.M*(2.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**4 + 2.0*0*Vars.v**5 + 2.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**6 + 2.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*Vars.v**7) + 6*Cons.E_4*Cons.M + 6*Vars.E_SQ_4*Cons.M*(2.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**5 + 2.0*0*Vars.v**6 + 2.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**7) + 7*Vars.E_SO_5*Cons.M*(2.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**6 + 2.0*0*Vars.v**7) + 2.0*8*Cons.E_6*Cons.M*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**7) + 4*Cons.E_2*Cons.M*(4*Cons.E_2*Cons.M*(1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**4 + 1.0*0*Vars.v**5 + 1.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**6 + 1.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*Vars.v**7) + 5*Vars.E_SO_3*Cons.M*(2.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**5 + 2.0*0*Vars.v**6 + 2.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**7) + 6*Cons.E_4*Cons.M + 6*Vars.E_SQ_4*Cons.M*(2.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**6 + 2.0*0*Vars.v**7) + 2.0*7*Vars.E_SO_5*Cons.M*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**7) + 5*Vars.E_SO_3*Cons.M*(5*Vars.E_SO_3*Cons.M*(1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**6 + 1.0*0*Vars.v**7) + 2.0*6*Cons.E_4*Cons.M + 6*Vars.E_SQ_4*Cons.M*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**7) + (0*(0*(0*(-1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**3 - 1.0*0*Vars.v**4 - 1.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**5 - 1.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*Vars.v**6 - 1.0*-Cons.Fcal_4 - Vars.Fcal_SQ_4*Vars.v**7) + 4*Cons.E_2*Cons.M*(-3.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**4 - 3.0*0*Vars.v**5 - 3.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**6 - 3.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*Vars.v**7) + 5*Vars.E_SO_3*Cons.M*(-3.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**5 - 3.0*0*Vars.v**6 - 3.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**7) + 6*Cons.E_4*Cons.M + 6*Vars.E_SQ_4*Cons.M*(-3.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**6 - 3.0*0*Vars.v**7) - 3.0*7*Vars.E_SO_5*Cons.M*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**7) + 4*Cons.E_2*Cons.M*(4*Cons.E_2*Cons.M*(-3.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**5 - 3.0*0*Vars.v**6 - 3.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**7) + 5*Vars.E_SO_3*Cons.M*(-6.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**6 - 6.0*0*Vars.v**7) - 6.0*6*Cons.E_4*Cons.M + 6*Vars.E_SQ_4*Cons.M*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**7) - 3.0*5*Vars.E_SO_3*Cons.M**2*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**7) + 4*Cons.E_2*Cons.M**2*(4*Cons.E_2*Cons.M*(-1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**6 - 1.0*0*Vars.v**7) - 3.0*5*Vars.E_SO_3*Cons.M*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**7) + (0*(0*(0*(0*(1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**4 + 1.0*0*Vars.v**5 + 1.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**6 + 1.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*Vars.v**7) + 4*Cons.E_2*Cons.M*(4.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**5 + 4.0*0*Vars.v**6 + 4.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**7) + 5*Vars.E_SO_3*Cons.M*(4.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**6 + 4.0*0*Vars.v**7) + 4.0*6*Cons.E_4*Cons.M + 6*Vars.E_SQ_4*Cons.M*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**7) + 4*Cons.E_2*Cons.M*(4*Cons.E_2*Cons.M*(6.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**6 + 6.0*0*Vars.v**7) + 12.0*5*Vars.E_SO_3*Cons.M*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**7)) + 4.0*4*Cons.E_2*Cons.M**3*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**7) + (0**3*(0*(0*(-1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**5 - 1.0*0*Vars.v**6 - 1.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**7) + 4*Cons.E_2*Cons.M*(-5.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**6 - 5.0*0*Vars.v**7) - 5.0*5*Vars.E_SO_3*Cons.M*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**7) - 10.0*4*Cons.E_2*Cons.M**2*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**7) + (0**5*(0*(1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**6 + 1.0*0*Vars.v**7) + 6.0*4*Cons.E_2*Cons.M*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**7) - 1.0*0**7*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**7/2*Cons.E_0*Cons.M)/2*Cons.E_0*Cons.M)/2*Cons.E_0*Cons.M)/2*Cons.E_0*Cons.M)/2*Cons.E_0*Cons.M)/2*Cons.E_0*Cons.M)/2*Cons.E_0*Cons.M)/(Cons.nu*Vars.v*2*Cons.E_0*Cons.M)
    dydt=np.zeros(10)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_3p5(Cons,Vars)[1:])
    dydt[0] = dvdt_T4[0]
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_3p5(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_3p5(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    dydt[9] = np.linalg.norm(OmegaVec_3p5(Cons,Vars)[1:])
    return dydt      

@njit(cache=True)
def TaylorT5_3p5(Cons,Vars):
    dtdv = -0.5*Cons.nu*Vars.v*(1.0*2*Cons.E_0*Cons.M + 1.0*0*Vars.v + 1.0*4*Cons.E_2*Cons.M*Vars.v**2 + 1.0*5*Vars.E_SO_3*Cons.M*Vars.v**3 + 1.0*6*Cons.E_4*Cons.M + 6*Vars.E_SQ_4*Cons.M*Vars.v**4 + 1.0*7*Vars.E_SO_5*Cons.M*Vars.v**5 + 1.0*8*Cons.E_6*Cons.M*Vars.v**6 + 1.0*9*Vars.E_SO_7*Cons.M*Vars.v**7 + (0*(-1.0*2*Cons.E_0*Cons.M*Vars.v - 1.0*0*Vars.v**2 - 1.0*4*Cons.E_2*Cons.M*Vars.v**3 - 1.0*5*Vars.E_SO_3*Cons.M*Vars.v**4 - 1.0*6*Cons.E_4*Cons.M + 6*Vars.E_SQ_4*Cons.M*Vars.v**5 - 1.0*7*Vars.E_SO_5*Cons.M*Vars.v**6 - 1.0*8*Cons.E_6*Cons.M*Vars.v**7) + -Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*(-1.0*2*Cons.E_0*Cons.M*Vars.v**2 - 1.0*0*Vars.v**3 - 1.0*4*Cons.E_2*Cons.M*Vars.v**4 - 1.0*5*Vars.E_SO_3*Cons.M*Vars.v**5 - 1.0*6*Cons.E_4*Cons.M + 6*Vars.E_SQ_4*Cons.M*Vars.v**6 - 1.0*7*Vars.E_SO_5*Cons.M*Vars.v**7) + -Cons.Fcal_3 - Vars.Fcal_SO_3*(-1.0*2*Cons.E_0*Cons.M*Vars.v**3 - 1.0*0*Vars.v**4 - 1.0*4*Cons.E_2*Cons.M*Vars.v**5 - 1.0*5*Vars.E_SO_3*Cons.M*Vars.v**6 - 1.0*6*Cons.E_4*Cons.M + 6*Vars.E_SQ_4*Cons.M*Vars.v**7) + -Cons.Fcal_4 - Vars.Fcal_SQ_4*(-1.0*2*Cons.E_0*Cons.M*Vars.v**4 - 1.0*0*Vars.v**5 - 1.0*4*Cons.E_2*Cons.M*Vars.v**6 - 1.0*5*Vars.E_SO_3*Cons.M*Vars.v**7) + -Cons.Fcal_5 - Vars.Fcal_SO_5*(-1.0*2*Cons.E_0*Cons.M*Vars.v**5 - 1.0*0*Vars.v**6 - 1.0*4*Cons.E_2*Cons.M*Vars.v**7) + -Cons.Fcal_6 - Vars.Fcal_SO_6 - Cons.Fcal_lnv_6*Vars.logv*(-1.0*2*Cons.E_0*Cons.M*Vars.v**6 - 1.0*0*Vars.v**7) - 1.0*-Cons.Fcal_7 - Vars.Fcal_SO_7*2*Cons.E_0*Cons.M*Vars.v**7 + (0*(0*(1.0*2*Cons.E_0*Cons.M*Vars.v**2 + 1.0*0*Vars.v**3 + 1.0*4*Cons.E_2*Cons.M*Vars.v**4 + 1.0*5*Vars.E_SO_3*Cons.M*Vars.v**5 + 1.0*6*Cons.E_4*Cons.M + 6*Vars.E_SQ_4*Cons.M*Vars.v**6 + 1.0*7*Vars.E_SO_5*Cons.M*Vars.v**7) + -Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*(2.0*2*Cons.E_0*Cons.M*Vars.v**3 + 2.0*0*Vars.v**4 + 2.0*4*Cons.E_2*Cons.M*Vars.v**5 + 2.0*5*Vars.E_SO_3*Cons.M*Vars.v**6 + 2.0*6*Cons.E_4*Cons.M + 6*Vars.E_SQ_4*Cons.M*Vars.v**7) + -Cons.Fcal_3 - Vars.Fcal_SO_3*(2.0*2*Cons.E_0*Cons.M*Vars.v**4 + 2.0*0*Vars.v**5 + 2.0*4*Cons.E_2*Cons.M*Vars.v**6 + 2.0*5*Vars.E_SO_3*Cons.M*Vars.v**7) + -Cons.Fcal_4 - Vars.Fcal_SQ_4*(2.0*2*Cons.E_0*Cons.M*Vars.v**5 + 2.0*0*Vars.v**6 + 2.0*4*Cons.E_2*Cons.M*Vars.v**7) + -Cons.Fcal_5 - Vars.Fcal_SO_5*(2.0*2*Cons.E_0*Cons.M*Vars.v**6 + 2.0*0*Vars.v**7) + 2.0*-Cons.Fcal_6 - Vars.Fcal_SO_6 - Cons.Fcal_lnv_6*Vars.logv*2*Cons.E_0*Cons.M*Vars.v**7) + -Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*(-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*(1.0*2*Cons.E_0*Cons.M*Vars.v**4 + 1.0*0*Vars.v**5 + 1.0*4*Cons.E_2*Cons.M*Vars.v**6 + 1.0*5*Vars.E_SO_3*Cons.M*Vars.v**7) + -Cons.Fcal_3 - Vars.Fcal_SO_3*(2.0*2*Cons.E_0*Cons.M*Vars.v**5 + 2.0*0*Vars.v**6 + 2.0*4*Cons.E_2*Cons.M*Vars.v**7) + -Cons.Fcal_4 - Vars.Fcal_SQ_4*(2.0*2*Cons.E_0*Cons.M*Vars.v**6 + 2.0*0*Vars.v**7) + 2.0*-Cons.Fcal_5 - Vars.Fcal_SO_5*2*Cons.E_0*Cons.M*Vars.v**7) + -Cons.Fcal_3 - Vars.Fcal_SO_3*(-Cons.Fcal_3 - Vars.Fcal_SO_3*(1.0*2*Cons.E_0*Cons.M*Vars.v**6 + 1.0*0*Vars.v**7) + 2.0*-Cons.Fcal_4 - Vars.Fcal_SQ_4*2*Cons.E_0*Cons.M*Vars.v**7) + (0*(0*(0*(-1.0*2*Cons.E_0*Cons.M*Vars.v**3 - 1.0*0*Vars.v**4 - 1.0*4*Cons.E_2*Cons.M*Vars.v**5 - 1.0*5*Vars.E_SO_3*Cons.M*Vars.v**6 - 1.0*6*Cons.E_4*Cons.M + 6*Vars.E_SQ_4*Cons.M*Vars.v**7) + -Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*(-3.0*2*Cons.E_0*Cons.M*Vars.v**4 - 3.0*0*Vars.v**5 - 3.0*4*Cons.E_2*Cons.M*Vars.v**6 - 3.0*5*Vars.E_SO_3*Cons.M*Vars.v**7) + -Cons.Fcal_3 - Vars.Fcal_SO_3*(-3.0*2*Cons.E_0*Cons.M*Vars.v**5 - 3.0*0*Vars.v**6 - 3.0*4*Cons.E_2*Cons.M*Vars.v**7) + -Cons.Fcal_4 - Vars.Fcal_SQ_4*(-3.0*2*Cons.E_0*Cons.M*Vars.v**6 - 3.0*0*Vars.v**7) - 3.0*-Cons.Fcal_5 - Vars.Fcal_SO_5*2*Cons.E_0*Cons.M*Vars.v**7) + -Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*(-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*(-3.0*2*Cons.E_0*Cons.M*Vars.v**5 - 3.0*0*Vars.v**6 - 3.0*4*Cons.E_2*Cons.M*Vars.v**7) + -Cons.Fcal_3 - Vars.Fcal_SO_3*(-6.0*2*Cons.E_0*Cons.M*Vars.v**6 - 6.0*0*Vars.v**7) - 6.0*-Cons.Fcal_4 - Vars.Fcal_SQ_4*2*Cons.E_0*Cons.M*Vars.v**7) - 3.0*-Cons.Fcal_3 - Vars.Fcal_SO_3**2*2*Cons.E_0*Cons.M*Vars.v**7) + -Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff**2*(-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*(-1.0*2*Cons.E_0*Cons.M*Vars.v**6 - 1.0*0*Vars.v**7) - 3.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*2*Cons.E_0*Cons.M*Vars.v**7) + (0*(0*(0*(0*(1.0*2*Cons.E_0*Cons.M*Vars.v**4 + 1.0*0*Vars.v**5 + 1.0*4*Cons.E_2*Cons.M*Vars.v**6 + 1.0*5*Vars.E_SO_3*Cons.M*Vars.v**7) + -Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*(4.0*2*Cons.E_0*Cons.M*Vars.v**5 + 4.0*0*Vars.v**6 + 4.0*4*Cons.E_2*Cons.M*Vars.v**7) + -Cons.Fcal_3 - Vars.Fcal_SO_3*(4.0*2*Cons.E_0*Cons.M*Vars.v**6 + 4.0*0*Vars.v**7) + 4.0*-Cons.Fcal_4 - Vars.Fcal_SQ_4*2*Cons.E_0*Cons.M*Vars.v**7) + -Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*(-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*(6.0*2*Cons.E_0*Cons.M*Vars.v**6 + 6.0*0*Vars.v**7) + 12.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*2*Cons.E_0*Cons.M*Vars.v**7)) + 4.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff**3*2*Cons.E_0*Cons.M*Vars.v**7) + (0**3*(0*(0*(-1.0*2*Cons.E_0*Cons.M*Vars.v**5 - 1.0*0*Vars.v**6 - 1.0*4*Cons.E_2*Cons.M*Vars.v**7) + -Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*(-5.0*2*Cons.E_0*Cons.M*Vars.v**6 - 5.0*0*Vars.v**7) - 5.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*2*Cons.E_0*Cons.M*Vars.v**7) - 10.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff**2*2*Cons.E_0*Cons.M*Vars.v**7) + (0**5*(0*(1.0*2*Cons.E_0*Cons.M*Vars.v**6 + 1.0*0*Vars.v**7) + 6.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*2*Cons.E_0*Cons.M*Vars.v**7) - 1.0*0**7*2*Cons.E_0*Cons.M*Vars.v**7/-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff)/-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff)/-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff)/-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff)/-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff)/-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff)/-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff)/(Vars.Fcal_coeff*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff)
    dvdt_T5 = 1.0/dtdv
    dydt=np.zeros(10)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_3p5(Cons,Vars)[1:])
    dydt[0] = dvdt_T5[0]
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_3p5(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_3p5(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    dydt[9] = np.linalg.norm(OmegaVec_3p5(Cons,Vars)[1:])
    return dydt      

@njit(cache=True)
def Recalculate_4p0(Cons,Vars,y):
    Vars.v = np.array([y[0]])
    Vars.rfrak_chi1 = np.array([y[1],y[2]])
    Vars.rfrak_chi2 = np.array([y[3],y[4]])
    Vars.Ecc = Cons.e_i*(Cons.v_i/Vars.v)**(19/6)
    Vars.xi = y[9]+Cons.xi_i
    Vars.R = y[5]*Cons.wHat + y[6]*Cons.xHat + y[7]*Cons.yHat + y[8]*Cons.zHat
    Vars.nHat = mul(mul(Vars.R,Cons.xHat),inverse(Vars.R))
    Vars.lambdaHat = mul(mul(Vars.R,Cons.yHat),inverse(Vars.R))
    Vars.ellHat = mul(mul(Vars.R,Cons.zHat),inverse(Vars.R))
    Vars.R_S1 = exp(Vars.rfrak_chi1[0]*Cons.xHat + Vars.rfrak_chi1[1]*Cons.yHat)
    Vars.R_S2 = exp(Vars.rfrak_chi2[0]*Cons.xHat + Vars.rfrak_chi2[1]*Cons.yHat)
    Vars.chiVec1 = mul(mul(mul(Cons.S_chi1,Vars.R_S1),Cons.zHat),mul(conjugate(Vars.R_S1),conjugate(Cons.S_chi1)))
    Vars.chiVec2 = mul(mul(mul(Cons.S_chi2,Vars.R_S2),Cons.zHat),mul(conjugate(Vars.R_S2),conjugate(Cons.S_chi2)))
    Vars.chi1_n = np.array([dot(Vars.chiVec1[1:],Vars.nHat[1:])])
    Vars.chi1_lambda = np.array([dot(Vars.chiVec1[1:],Vars.lambdaHat[1:])])
    Vars.chi1_ell = np.array([dot(Vars.chiVec1[1:],Vars.ellHat[1:])])
    Vars.chi2_n = np.array([dot(Vars.chiVec2[1:],Vars.nHat[1:])])
    Vars.chi2_lambda = np.array([dot(Vars.chiVec2[1:],Vars.lambdaHat[1:])])
    Vars.chi2_ell = np.array([dot(Vars.chiVec2[1:],Vars.ellHat[1:])])
    Vars.S_ell = Cons.M1**2*Vars.chi1_ell + Cons.M2**2*Vars.chi2_ell
    Vars.S_n = Cons.M1**2*Vars.chi1_n + Cons.M2**2*Vars.chi2_n
    Vars.S_lambda = Cons.M1**2*Vars.chi1_lambda + Cons.M2**2*Vars.chi2_lambda
    Vars.Sigma_ell = Cons.M*(-Cons.M1*Vars.chi1_ell + Cons.M2*Vars.chi2_ell)
    Vars.Sigma_n = Cons.M*(-Cons.M1*Vars.chi1_n + Cons.M2*Vars.chi2_n)
    Vars.Sigma_lambda = Cons.M*(-Cons.M1*Vars.chi1_lambda + Cons.M2*Vars.chi2_lambda)
    Vars.chi_s_ell = Vars.chi1_ell/2 + Vars.chi2_ell/2
    Vars.chi_a_ell = Vars.chi1_ell/2 - Vars.chi2_ell/2
    Vars.logv = log(Vars.v)
    Vars.Fcal_coeff = 32*Cons.nu**2*Vars.v**10/5
    Vars.Fcal_SQ_4 = Cons.chi1chi1*(-89*Cons.delta/192 + 89*Cons.nu/96 - 89/192) - 103*Cons.chi1chi2*Cons.nu/48 + Cons.chi2chi2*(89*Cons.delta/192 + 89*Cons.nu/96 - 89/192) + Vars.chi_a_ell*(Vars.chi_a_ell*(287/96 - 12*Cons.nu) + 287*Vars.chi_s_ell*Cons.delta/48) + Vars.chi_s_ell**2*(Cons.nu/24 + 287/96)
    Vars.Fcal_SO_3 = (-4*Vars.S_ell - 5*Vars.Sigma_ell*Cons.delta/4)/Cons.M**2
    Vars.Fcal_SO_5 = (Vars.S_ell*(272*Cons.nu/9 - 9/2) + Vars.Sigma_ell*Cons.delta*(43*Cons.nu/4 - 13/16))/Cons.M**2
    Vars.Fcal_SO_6 = (-16*Vars.S_ell*pi - 31*Vars.Sigma_ell*Cons.delta*pi/6)/Cons.M**2
    Vars.Fcal_SO_7 = (Vars.S_ell*(-2810*Cons.nu**2/27 + 6172*Cons.nu/189 + 476645/6804) + Vars.Sigma_ell*Cons.delta*(-1501*Cons.nu**2/36 + 1849*Cons.nu/126 + 9535/336))/Cons.M**2
    Vars.Fcal_SO_8 = (Vars.S_ell*pi*(13879*Cons.nu/72 - 3485/96) + Vars.Sigma_ell*Cons.delta*pi*(130583*Cons.nu/2016 - 7163/672))/Cons.M**2
    Vars.Fecc_coeff = 157*Cons.e_i**2.0*Cons.v_i**6.33333333333333/(24*Vars.v**6.33333333333333)
    Vars.E_SQ_4 = -3*Vars.chi_a_ell**2/2 - 3*Vars.chi_s_ell**2/2 - Cons.delta*(Cons.chi2chi2/2 + 3*Vars.chi_a_ell*Vars.chi_s_ell) + Cons.nu*(Cons.chi1chi2 + 6*Vars.chi_a_ell**2) + (Cons.chi1chi1 + Cons.chi2chi2)*(Cons.delta - 2*Cons.nu + 1)/4
    Vars.E_SO_3 = (14*Vars.S_ell/3 + 2*Vars.Sigma_ell*Cons.delta)/Cons.M**2
    Vars.E_SO_5 = (Vars.S_ell*(11 - 61*Cons.nu/9) + Vars.Sigma_ell*Cons.delta*(3 - 10*Cons.nu/3))/Cons.M**2
    Vars.E_SO_7 = (Vars.S_ell*(29*Cons.nu**2/12 - 367*Cons.nu/4 + 135/4) + Vars.Sigma_ell*Cons.delta*(5*Cons.nu**2/4 - 39*Cons.nu + 27/4))/Cons.M**2

@njit
def OmegaVec_chiVec_1_4p0(Cons,Vars):
    Omega1_coeff = Vars.v**5/Cons.M
    return Omega1_coeff*(Vars.ellHat*(-0.75*Cons.delta + 0.5*Cons.nu + Vars.v**2*(Cons.delta*(0.625*Cons.nu - 0.5625) + Cons.nu*(1.25 - 0.0416666666666667*Cons.nu) + Vars.v**2*(Cons.delta*(Cons.nu*(4.875 - 0.15625*Cons.nu) - 0.84375) + Cons.nu*(Cons.nu*(-0.0208333333333333*Cons.nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + Vars.nHat*Vars.v*(3.0*Vars.chi1_n*Cons.nu + 3.0*Cons.M2**2*Vars.chi2_n/Cons.M**2) - Cons.M2**2*Vars.chiVec2*Vars.v/Cons.M**2)

@njit
def OmegaVec_chiVec_2_4p0(Cons,Vars):
    Omega2_coeff = Vars.v**5/Cons.M
    return Omega2_coeff*(Vars.ellHat*(0.75*Cons.delta + 0.5*Cons.nu + Vars.v**2*(Cons.delta*(0.5625 - 0.625*Cons.nu) + Cons.nu*(1.25 - 0.0416666666666667*Cons.nu) + Vars.v**2*(Cons.delta*(Cons.nu*(0.15625*Cons.nu - 4.875) + 0.84375) + Cons.nu*(Cons.nu*(-0.0208333333333333*Cons.nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75) + Vars.nHat*Vars.v*(3.0*Vars.chi2_n*Cons.nu + 3.0*Cons.M1**2*Vars.chi1_n/Cons.M**2) - Cons.M1**2*Vars.chiVec1*Vars.v/Cons.M**2)

@njit
def OmegaVec_4p0(Cons,Vars):
    gamma_PN_0 = 1.00000000000000
    a_ell_0 = 7.0*Vars.S_n + 3.0*Vars.Sigma_n*Cons.delta
    gamma_PN_3 = (1.66666666666667*Vars.S_ell + Vars.Sigma_ell*Cons.delta)/Cons.M**2
    a_ell_4 = Vars.S_n*(5.77777777777778*Cons.nu**2 + 14.75*Cons.nu + 1.5) + Vars.Sigma_n*Cons.delta*(2.83333333333333*Cons.nu**2 + 9.125*Cons.nu + 1.5)
    a_ell_2 = Vars.S_n*(-9.66666666666667*Cons.nu - 10.0) + Vars.Sigma_n*Cons.delta*(-4.5*Cons.nu - 6.0)
    gamma_PN_5 = (Vars.S_ell*(0.888888888888889*Cons.nu + 3.33333333333333) + 2.0*Vars.Sigma_ell*Cons.delta)/Cons.M**2
    gamma_PN_7 = (Vars.S_ell*(-6.0*Cons.nu**2 - 10.5833333333333*Cons.nu + 5.0) - 2.66666666666667*Vars.Sigma_ell*Cons.delta*Cons.nu**2 + Vars.Sigma_ell*Cons.delta*(3.0 - 10.1666666666667*Cons.nu))/Cons.M**2
    gamma_PN_4 = 1.0 - 5.41666666666667*Cons.nu
    gamma_PN_6 = 0.0123456790123457*Cons.nu**3 + 6.36111111111111*Cons.nu**2 - 2.98177812235564*Cons.nu + 1.0
    gamma_PN_2 = 1.0 - 0.333333333333333*Cons.nu
    return Vars.ellHat*Vars.v**3*np.sqrt(1.0 - Vars.Ecc**2.0)*(-Vars.Ecc*np.cos(Vars.xi) + 1.0)**(-2.0)/Cons.M + Vars.nHat*Vars.v**6*(a_ell_0 + Vars.v**2*(a_ell_2 + a_ell_4*Vars.v**2))*(gamma_PN_0 + Vars.v**2*(gamma_PN_2 + Vars.v*(gamma_PN_3 + Vars.v*(gamma_PN_4 + Vars.v*(gamma_PN_5 + Vars.v*(gamma_PN_6 + gamma_PN_7*Vars.v))))))/Cons.M**3


@njit(cache=True)
def TaylorT1_4p0(Cons,Vars):
    Flux = Vars.Fcal_coeff*Vars.Fecc_coeff*(Cons.Fecc_0 + Cons.Fecc_2*Vars.v**2) + Vars.Fcal_coeff*(Cons.Fcal_0 + Vars.v**2*(Cons.Fcal_2 + Vars.v*(Cons.Fcal_3 + Vars.Fcal_SO_3 + Vars.v*(Cons.Fcal_4 + Vars.Fcal_SQ_4 + Vars.v*(Cons.Fcal_5 + Vars.Fcal_SO_5 + Vars.v*(Cons.Fcal_6 + Vars.Fcal_SO_6 + Cons.Fcal_lnv_6*Vars.logv + Vars.v*(Cons.Fcal_7 + Vars.Fcal_SO_7 + Vars.v*(Cons.Fcal_8 + Vars.Fcal_SO_8 + Cons.Fcal_lnv_8*Vars.logv))))))))
    dEdV = 1.16666666666667*Cons.Eecc_2*Cons.M*Cons.e_i**2.0*Cons.nu*Vars.v**(-3.33333333333333)*Cons.v_i**6.33333333333333 - 0.5*Cons.M*Cons.nu*Vars.v*(2.0*Cons.E_0 + Vars.v**2*(4.0*Cons.E_2 + Vars.v*(5.0*Vars.E_SO_3 + Vars.v*(6.0*Cons.E_4 + 6.0*Vars.E_SQ_4 + Vars.v*(7.0*Vars.E_SO_5 + Vars.v*(8.0*Cons.E_6 + Vars.v*(9.0*Vars.E_SO_7 + Vars.v*(10.0*Cons.E_8 + Cons.E_lnv_8*(10.0*Vars.logv + 1.0)))))))))
    Absorption = 0
    dvdt_T1 = (-Absorption - Flux)/dEdV
    dydt=np.zeros(10)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_4p0(Cons,Vars)[1:])
    dydt[0] = dvdt_T1[0] 
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_4p0(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_4p0(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    dydt[9] = np.linalg.norm(OmegaVec_4p0(Cons,Vars)[1:])
    return dydt

@njit(cache=True)
def TaylorT4_4p0(Cons,Vars):
    dvdt_T4 = -2.0*Vars.Fcal_coeff*(1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff + 1.0*0*Vars.v + 1.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**2 + 1.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*Vars.v**3 + 1.0*-Cons.Fcal_4 - Vars.Fcal_SQ_4*Vars.v**4 + 1.0*-Cons.Fcal_5 - Vars.Fcal_SO_5*Vars.v**5 + 1.0*-Cons.Fcal_6 - Vars.Fcal_SO_6 - Cons.Fcal_lnv_6*Vars.logv*Vars.v**6 + 1.0*-Cons.Fcal_7 - Vars.Fcal_SO_7*Vars.v**7 + 1.0*-Cons.Fcal_8 - Vars.Fcal_SO_8 - Cons.Fcal_lnv_8*Vars.logv*Vars.v**8 + (0*(-1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v - 1.0*0*Vars.v**2 - 1.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**3 - 1.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*Vars.v**4 - 1.0*-Cons.Fcal_4 - Vars.Fcal_SQ_4*Vars.v**5 - 1.0*-Cons.Fcal_5 - Vars.Fcal_SO_5*Vars.v**6 - 1.0*-Cons.Fcal_6 - Vars.Fcal_SO_6 - Cons.Fcal_lnv_6*Vars.logv*Vars.v**7 - 1.0*-Cons.Fcal_7 - Vars.Fcal_SO_7*Vars.v**8) + 4*Cons.E_2*Cons.M*(-1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**2 - 1.0*0*Vars.v**3 - 1.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**4 - 1.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*Vars.v**5 - 1.0*-Cons.Fcal_4 - Vars.Fcal_SQ_4*Vars.v**6 - 1.0*-Cons.Fcal_5 - Vars.Fcal_SO_5*Vars.v**7 - 1.0*-Cons.Fcal_6 - Vars.Fcal_SO_6 - Cons.Fcal_lnv_6*Vars.logv*Vars.v**8) + 5*Vars.E_SO_3*Cons.M*(-1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**3 - 1.0*0*Vars.v**4 - 1.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**5 - 1.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*Vars.v**6 - 1.0*-Cons.Fcal_4 - Vars.Fcal_SQ_4*Vars.v**7 - 1.0*-Cons.Fcal_5 - Vars.Fcal_SO_5*Vars.v**8) + 6*Cons.E_4*Cons.M + 6*Vars.E_SQ_4*Cons.M*(-1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**4 - 1.0*0*Vars.v**5 - 1.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**6 - 1.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*Vars.v**7 - 1.0*-Cons.Fcal_4 - Vars.Fcal_SQ_4*Vars.v**8) + 7*Vars.E_SO_5*Cons.M*(-1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**5 - 1.0*0*Vars.v**6 - 1.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**7 - 1.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*Vars.v**8) + 8*Cons.E_6*Cons.M*(-1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**6 - 1.0*0*Vars.v**7 - 1.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**8) + 9*Vars.E_SO_7*Cons.M*(-1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**7 - 1.0*0*Vars.v**8) - 1.0*10*Cons.E_8*Cons.M + Cons.E_lnv_8*Cons.M*(10*Vars.logv + 1)*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**8 + (0*(0*(1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**2 + 1.0*0*Vars.v**3 + 1.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**4 + 1.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*Vars.v**5 + 1.0*-Cons.Fcal_4 - Vars.Fcal_SQ_4*Vars.v**6 + 1.0*-Cons.Fcal_5 - Vars.Fcal_SO_5*Vars.v**7 + 1.0*-Cons.Fcal_6 - Vars.Fcal_SO_6 - Cons.Fcal_lnv_6*Vars.logv*Vars.v**8) + 4*Cons.E_2*Cons.M*(2.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**3 + 2.0*0*Vars.v**4 + 2.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**5 + 2.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*Vars.v**6 + 2.0*-Cons.Fcal_4 - Vars.Fcal_SQ_4*Vars.v**7 + 2.0*-Cons.Fcal_5 - Vars.Fcal_SO_5*Vars.v**8) + 5*Vars.E_SO_3*Cons.M*(2.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**4 + 2.0*0*Vars.v**5 + 2.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**6 + 2.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*Vars.v**7 + 2.0*-Cons.Fcal_4 - Vars.Fcal_SQ_4*Vars.v**8) + 6*Cons.E_4*Cons.M + 6*Vars.E_SQ_4*Cons.M*(2.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**5 + 2.0*0*Vars.v**6 + 2.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**7 + 2.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*Vars.v**8) + 7*Vars.E_SO_5*Cons.M*(2.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**6 + 2.0*0*Vars.v**7 + 2.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**8) + 8*Cons.E_6*Cons.M*(2.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**7 + 2.0*0*Vars.v**8) + 2.0*9*Vars.E_SO_7*Cons.M*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**8) + 4*Cons.E_2*Cons.M*(4*Cons.E_2*Cons.M*(1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**4 + 1.0*0*Vars.v**5 + 1.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**6 + 1.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*Vars.v**7 + 1.0*-Cons.Fcal_4 - Vars.Fcal_SQ_4*Vars.v**8) + 5*Vars.E_SO_3*Cons.M*(2.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**5 + 2.0*0*Vars.v**6 + 2.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**7 + 2.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*Vars.v**8) + 6*Cons.E_4*Cons.M + 6*Vars.E_SQ_4*Cons.M*(2.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**6 + 2.0*0*Vars.v**7 + 2.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**8) + 7*Vars.E_SO_5*Cons.M*(2.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**7 + 2.0*0*Vars.v**8) + 2.0*8*Cons.E_6*Cons.M*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**8) + 5*Vars.E_SO_3*Cons.M*(5*Vars.E_SO_3*Cons.M*(1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**6 + 1.0*0*Vars.v**7 + 1.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**8) + 6*Cons.E_4*Cons.M + 6*Vars.E_SQ_4*Cons.M*(2.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**7 + 2.0*0*Vars.v**8) + 2.0*7*Vars.E_SO_5*Cons.M*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**8) + 1.0*6*Cons.E_4*Cons.M + 6*Vars.E_SQ_4*Cons.M**2*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**8 + (0*(0*(0*(-1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**3 - 1.0*0*Vars.v**4 - 1.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**5 - 1.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*Vars.v**6 - 1.0*-Cons.Fcal_4 - Vars.Fcal_SQ_4*Vars.v**7 - 1.0*-Cons.Fcal_5 - Vars.Fcal_SO_5*Vars.v**8) + 4*Cons.E_2*Cons.M*(-3.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**4 - 3.0*0*Vars.v**5 - 3.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**6 - 3.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*Vars.v**7 - 3.0*-Cons.Fcal_4 - Vars.Fcal_SQ_4*Vars.v**8) + 5*Vars.E_SO_3*Cons.M*(-3.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**5 - 3.0*0*Vars.v**6 - 3.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**7 - 3.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*Vars.v**8) + 6*Cons.E_4*Cons.M + 6*Vars.E_SQ_4*Cons.M*(-3.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**6 - 3.0*0*Vars.v**7 - 3.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**8) + 7*Vars.E_SO_5*Cons.M*(-3.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**7 - 3.0*0*Vars.v**8) - 3.0*8*Cons.E_6*Cons.M*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**8) + 4*Cons.E_2*Cons.M*(4*Cons.E_2*Cons.M*(-3.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**5 - 3.0*0*Vars.v**6 - 3.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**7 - 3.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*Vars.v**8) + 5*Vars.E_SO_3*Cons.M*(-6.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**6 - 6.0*0*Vars.v**7 - 6.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**8) + 6*Cons.E_4*Cons.M + 6*Vars.E_SQ_4*Cons.M*(-6.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**7 - 6.0*0*Vars.v**8) - 6.0*7*Vars.E_SO_5*Cons.M*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**8) + 5*Vars.E_SO_3*Cons.M*(5*Vars.E_SO_3*Cons.M*(-3.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**7 - 3.0*0*Vars.v**8) - 6.0*6*Cons.E_4*Cons.M + 6*Vars.E_SQ_4*Cons.M*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**8)) + 4*Cons.E_2*Cons.M*(4*Cons.E_2*Cons.M*(4*Cons.E_2*Cons.M*(-1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**6 - 1.0*0*Vars.v**7 - 1.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**8) + 5*Vars.E_SO_3*Cons.M*(-3.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**7 - 3.0*0*Vars.v**8) - 3.0*6*Cons.E_4*Cons.M + 6*Vars.E_SQ_4*Cons.M*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**8) - 3.0*5*Vars.E_SO_3*Cons.M**2*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**8) + (0*(0*(0*(0*(1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**4 + 1.0*0*Vars.v**5 + 1.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**6 + 1.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*Vars.v**7 + 1.0*-Cons.Fcal_4 - Vars.Fcal_SQ_4*Vars.v**8) + 4*Cons.E_2*Cons.M*(4.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**5 + 4.0*0*Vars.v**6 + 4.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**7 + 4.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*Vars.v**8) + 5*Vars.E_SO_3*Cons.M*(4.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**6 + 4.0*0*Vars.v**7 + 4.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**8) + 6*Cons.E_4*Cons.M + 6*Vars.E_SQ_4*Cons.M*(4.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**7 + 4.0*0*Vars.v**8) + 4.0*7*Vars.E_SO_5*Cons.M*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**8) + 4*Cons.E_2*Cons.M*(4*Cons.E_2*Cons.M*(6.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**6 + 6.0*0*Vars.v**7 + 6.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**8) + 5*Vars.E_SO_3*Cons.M*(12.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**7 + 12.0*0*Vars.v**8) + 12.0*6*Cons.E_4*Cons.M + 6*Vars.E_SQ_4*Cons.M*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**8) + 6.0*5*Vars.E_SO_3*Cons.M**2*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**8) + 4*Cons.E_2*Cons.M**2*(4*Cons.E_2*Cons.M*(4.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**7 + 4.0*0*Vars.v**8) + 12.0*5*Vars.E_SO_3*Cons.M*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**8)) + 1.0*4*Cons.E_2*Cons.M**4*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**8 + (0**2*(0*(0*(0*(-1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**5 - 1.0*0*Vars.v**6 - 1.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**7 - 1.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*Vars.v**8) + 4*Cons.E_2*Cons.M*(-5.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**6 - 5.0*0*Vars.v**7 - 5.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**8) + 5*Vars.E_SO_3*Cons.M*(-5.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**7 - 5.0*0*Vars.v**8) - 5.0*6*Cons.E_4*Cons.M + 6*Vars.E_SQ_4*Cons.M*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**8) + 4*Cons.E_2*Cons.M*(4*Cons.E_2*Cons.M*(-10.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**7 - 10.0*0*Vars.v**8) - 20.0*5*Vars.E_SO_3*Cons.M*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**8)) - 10.0*4*Cons.E_2*Cons.M**3*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**8) + (0**4*(0*(0*(1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**6 + 1.0*0*Vars.v**7 + 1.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*Vars.v**8) + 4*Cons.E_2*Cons.M*(6.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**7 + 6.0*0*Vars.v**8) + 6.0*5*Vars.E_SO_3*Cons.M*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**8) + 15.0*4*Cons.E_2*Cons.M**2*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**8) + (0**6*(0*(-1.0*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**7 - 1.0*0*Vars.v**8) - 7.0*4*Cons.E_2*Cons.M*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**8) + 1.0*0**8*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff*Vars.v**8/2*Cons.E_0*Cons.M)/2*Cons.E_0*Cons.M)/2*Cons.E_0*Cons.M)/2*Cons.E_0*Cons.M)/2*Cons.E_0*Cons.M)/2*Cons.E_0*Cons.M)/2*Cons.E_0*Cons.M)/2*Cons.E_0*Cons.M)/(Cons.nu*Vars.v*2*Cons.E_0*Cons.M)
    dydt=np.zeros(10)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_4p0(Cons,Vars)[1:])
    dydt[0] = dvdt_T4[0]
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_4p0(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_4p0(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    dydt[9] = np.linalg.norm(OmegaVec_4p0(Cons,Vars)[1:])
    return dydt      

@njit(cache=True)
def TaylorT5_4p0(Cons,Vars):
    dtdv = -0.5*Cons.nu*Vars.v*(1.0*2*Cons.E_0*Cons.M + 1.0*0*Vars.v + 1.0*4*Cons.E_2*Cons.M*Vars.v**2 + 1.0*5*Vars.E_SO_3*Cons.M*Vars.v**3 + 1.0*6*Cons.E_4*Cons.M + 6*Vars.E_SQ_4*Cons.M*Vars.v**4 + 1.0*7*Vars.E_SO_5*Cons.M*Vars.v**5 + 1.0*8*Cons.E_6*Cons.M*Vars.v**6 + 1.0*9*Vars.E_SO_7*Cons.M*Vars.v**7 + 1.0*10*Cons.E_8*Cons.M + Cons.E_lnv_8*Cons.M*(10*Vars.logv + 1)*Vars.v**8 + (0*(-1.0*2*Cons.E_0*Cons.M*Vars.v - 1.0*0*Vars.v**2 - 1.0*4*Cons.E_2*Cons.M*Vars.v**3 - 1.0*5*Vars.E_SO_3*Cons.M*Vars.v**4 - 1.0*6*Cons.E_4*Cons.M + 6*Vars.E_SQ_4*Cons.M*Vars.v**5 - 1.0*7*Vars.E_SO_5*Cons.M*Vars.v**6 - 1.0*8*Cons.E_6*Cons.M*Vars.v**7 - 1.0*9*Vars.E_SO_7*Cons.M*Vars.v**8) + -Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*(-1.0*2*Cons.E_0*Cons.M*Vars.v**2 - 1.0*0*Vars.v**3 - 1.0*4*Cons.E_2*Cons.M*Vars.v**4 - 1.0*5*Vars.E_SO_3*Cons.M*Vars.v**5 - 1.0*6*Cons.E_4*Cons.M + 6*Vars.E_SQ_4*Cons.M*Vars.v**6 - 1.0*7*Vars.E_SO_5*Cons.M*Vars.v**7 - 1.0*8*Cons.E_6*Cons.M*Vars.v**8) + -Cons.Fcal_3 - Vars.Fcal_SO_3*(-1.0*2*Cons.E_0*Cons.M*Vars.v**3 - 1.0*0*Vars.v**4 - 1.0*4*Cons.E_2*Cons.M*Vars.v**5 - 1.0*5*Vars.E_SO_3*Cons.M*Vars.v**6 - 1.0*6*Cons.E_4*Cons.M + 6*Vars.E_SQ_4*Cons.M*Vars.v**7 - 1.0*7*Vars.E_SO_5*Cons.M*Vars.v**8) + -Cons.Fcal_4 - Vars.Fcal_SQ_4*(-1.0*2*Cons.E_0*Cons.M*Vars.v**4 - 1.0*0*Vars.v**5 - 1.0*4*Cons.E_2*Cons.M*Vars.v**6 - 1.0*5*Vars.E_SO_3*Cons.M*Vars.v**7 - 1.0*6*Cons.E_4*Cons.M + 6*Vars.E_SQ_4*Cons.M*Vars.v**8) + -Cons.Fcal_5 - Vars.Fcal_SO_5*(-1.0*2*Cons.E_0*Cons.M*Vars.v**5 - 1.0*0*Vars.v**6 - 1.0*4*Cons.E_2*Cons.M*Vars.v**7 - 1.0*5*Vars.E_SO_3*Cons.M*Vars.v**8) + -Cons.Fcal_6 - Vars.Fcal_SO_6 - Cons.Fcal_lnv_6*Vars.logv*(-1.0*2*Cons.E_0*Cons.M*Vars.v**6 - 1.0*0*Vars.v**7 - 1.0*4*Cons.E_2*Cons.M*Vars.v**8) + -Cons.Fcal_7 - Vars.Fcal_SO_7*(-1.0*2*Cons.E_0*Cons.M*Vars.v**7 - 1.0*0*Vars.v**8) - 1.0*-Cons.Fcal_8 - Vars.Fcal_SO_8 - Cons.Fcal_lnv_8*Vars.logv*2*Cons.E_0*Cons.M*Vars.v**8 + (0*(0*(1.0*2*Cons.E_0*Cons.M*Vars.v**2 + 1.0*0*Vars.v**3 + 1.0*4*Cons.E_2*Cons.M*Vars.v**4 + 1.0*5*Vars.E_SO_3*Cons.M*Vars.v**5 + 1.0*6*Cons.E_4*Cons.M + 6*Vars.E_SQ_4*Cons.M*Vars.v**6 + 1.0*7*Vars.E_SO_5*Cons.M*Vars.v**7 + 1.0*8*Cons.E_6*Cons.M*Vars.v**8) + -Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*(2.0*2*Cons.E_0*Cons.M*Vars.v**3 + 2.0*0*Vars.v**4 + 2.0*4*Cons.E_2*Cons.M*Vars.v**5 + 2.0*5*Vars.E_SO_3*Cons.M*Vars.v**6 + 2.0*6*Cons.E_4*Cons.M + 6*Vars.E_SQ_4*Cons.M*Vars.v**7 + 2.0*7*Vars.E_SO_5*Cons.M*Vars.v**8) + -Cons.Fcal_3 - Vars.Fcal_SO_3*(2.0*2*Cons.E_0*Cons.M*Vars.v**4 + 2.0*0*Vars.v**5 + 2.0*4*Cons.E_2*Cons.M*Vars.v**6 + 2.0*5*Vars.E_SO_3*Cons.M*Vars.v**7 + 2.0*6*Cons.E_4*Cons.M + 6*Vars.E_SQ_4*Cons.M*Vars.v**8) + -Cons.Fcal_4 - Vars.Fcal_SQ_4*(2.0*2*Cons.E_0*Cons.M*Vars.v**5 + 2.0*0*Vars.v**6 + 2.0*4*Cons.E_2*Cons.M*Vars.v**7 + 2.0*5*Vars.E_SO_3*Cons.M*Vars.v**8) + -Cons.Fcal_5 - Vars.Fcal_SO_5*(2.0*2*Cons.E_0*Cons.M*Vars.v**6 + 2.0*0*Vars.v**7 + 2.0*4*Cons.E_2*Cons.M*Vars.v**8) + -Cons.Fcal_6 - Vars.Fcal_SO_6 - Cons.Fcal_lnv_6*Vars.logv*(2.0*2*Cons.E_0*Cons.M*Vars.v**7 + 2.0*0*Vars.v**8) + 2.0*-Cons.Fcal_7 - Vars.Fcal_SO_7*2*Cons.E_0*Cons.M*Vars.v**8) + -Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*(-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*(1.0*2*Cons.E_0*Cons.M*Vars.v**4 + 1.0*0*Vars.v**5 + 1.0*4*Cons.E_2*Cons.M*Vars.v**6 + 1.0*5*Vars.E_SO_3*Cons.M*Vars.v**7 + 1.0*6*Cons.E_4*Cons.M + 6*Vars.E_SQ_4*Cons.M*Vars.v**8) + -Cons.Fcal_3 - Vars.Fcal_SO_3*(2.0*2*Cons.E_0*Cons.M*Vars.v**5 + 2.0*0*Vars.v**6 + 2.0*4*Cons.E_2*Cons.M*Vars.v**7 + 2.0*5*Vars.E_SO_3*Cons.M*Vars.v**8) + -Cons.Fcal_4 - Vars.Fcal_SQ_4*(2.0*2*Cons.E_0*Cons.M*Vars.v**6 + 2.0*0*Vars.v**7 + 2.0*4*Cons.E_2*Cons.M*Vars.v**8) + -Cons.Fcal_5 - Vars.Fcal_SO_5*(2.0*2*Cons.E_0*Cons.M*Vars.v**7 + 2.0*0*Vars.v**8) + 2.0*-Cons.Fcal_6 - Vars.Fcal_SO_6 - Cons.Fcal_lnv_6*Vars.logv*2*Cons.E_0*Cons.M*Vars.v**8) + -Cons.Fcal_3 - Vars.Fcal_SO_3*(-Cons.Fcal_3 - Vars.Fcal_SO_3*(1.0*2*Cons.E_0*Cons.M*Vars.v**6 + 1.0*0*Vars.v**7 + 1.0*4*Cons.E_2*Cons.M*Vars.v**8) + -Cons.Fcal_4 - Vars.Fcal_SQ_4*(2.0*2*Cons.E_0*Cons.M*Vars.v**7 + 2.0*0*Vars.v**8) + 2.0*-Cons.Fcal_5 - Vars.Fcal_SO_5*2*Cons.E_0*Cons.M*Vars.v**8) + 1.0*-Cons.Fcal_4 - Vars.Fcal_SQ_4**2*2*Cons.E_0*Cons.M*Vars.v**8 + (0*(0*(0*(-1.0*2*Cons.E_0*Cons.M*Vars.v**3 - 1.0*0*Vars.v**4 - 1.0*4*Cons.E_2*Cons.M*Vars.v**5 - 1.0*5*Vars.E_SO_3*Cons.M*Vars.v**6 - 1.0*6*Cons.E_4*Cons.M + 6*Vars.E_SQ_4*Cons.M*Vars.v**7 - 1.0*7*Vars.E_SO_5*Cons.M*Vars.v**8) + -Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*(-3.0*2*Cons.E_0*Cons.M*Vars.v**4 - 3.0*0*Vars.v**5 - 3.0*4*Cons.E_2*Cons.M*Vars.v**6 - 3.0*5*Vars.E_SO_3*Cons.M*Vars.v**7 - 3.0*6*Cons.E_4*Cons.M + 6*Vars.E_SQ_4*Cons.M*Vars.v**8) + -Cons.Fcal_3 - Vars.Fcal_SO_3*(-3.0*2*Cons.E_0*Cons.M*Vars.v**5 - 3.0*0*Vars.v**6 - 3.0*4*Cons.E_2*Cons.M*Vars.v**7 - 3.0*5*Vars.E_SO_3*Cons.M*Vars.v**8) + -Cons.Fcal_4 - Vars.Fcal_SQ_4*(-3.0*2*Cons.E_0*Cons.M*Vars.v**6 - 3.0*0*Vars.v**7 - 3.0*4*Cons.E_2*Cons.M*Vars.v**8) + -Cons.Fcal_5 - Vars.Fcal_SO_5*(-3.0*2*Cons.E_0*Cons.M*Vars.v**7 - 3.0*0*Vars.v**8) - 3.0*-Cons.Fcal_6 - Vars.Fcal_SO_6 - Cons.Fcal_lnv_6*Vars.logv*2*Cons.E_0*Cons.M*Vars.v**8) + -Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*(-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*(-3.0*2*Cons.E_0*Cons.M*Vars.v**5 - 3.0*0*Vars.v**6 - 3.0*4*Cons.E_2*Cons.M*Vars.v**7 - 3.0*5*Vars.E_SO_3*Cons.M*Vars.v**8) + -Cons.Fcal_3 - Vars.Fcal_SO_3*(-6.0*2*Cons.E_0*Cons.M*Vars.v**6 - 6.0*0*Vars.v**7 - 6.0*4*Cons.E_2*Cons.M*Vars.v**8) + -Cons.Fcal_4 - Vars.Fcal_SQ_4*(-6.0*2*Cons.E_0*Cons.M*Vars.v**7 - 6.0*0*Vars.v**8) - 6.0*-Cons.Fcal_5 - Vars.Fcal_SO_5*2*Cons.E_0*Cons.M*Vars.v**8) + -Cons.Fcal_3 - Vars.Fcal_SO_3*(-Cons.Fcal_3 - Vars.Fcal_SO_3*(-3.0*2*Cons.E_0*Cons.M*Vars.v**7 - 3.0*0*Vars.v**8) - 6.0*-Cons.Fcal_4 - Vars.Fcal_SQ_4*2*Cons.E_0*Cons.M*Vars.v**8)) + -Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*(-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*(-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*(-1.0*2*Cons.E_0*Cons.M*Vars.v**6 - 1.0*0*Vars.v**7 - 1.0*4*Cons.E_2*Cons.M*Vars.v**8) + -Cons.Fcal_3 - Vars.Fcal_SO_3*(-3.0*2*Cons.E_0*Cons.M*Vars.v**7 - 3.0*0*Vars.v**8) - 3.0*-Cons.Fcal_4 - Vars.Fcal_SQ_4*2*Cons.E_0*Cons.M*Vars.v**8) - 3.0*-Cons.Fcal_3 - Vars.Fcal_SO_3**2*2*Cons.E_0*Cons.M*Vars.v**8) + (0*(0*(0*(0*(1.0*2*Cons.E_0*Cons.M*Vars.v**4 + 1.0*0*Vars.v**5 + 1.0*4*Cons.E_2*Cons.M*Vars.v**6 + 1.0*5*Vars.E_SO_3*Cons.M*Vars.v**7 + 1.0*6*Cons.E_4*Cons.M + 6*Vars.E_SQ_4*Cons.M*Vars.v**8) + -Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*(4.0*2*Cons.E_0*Cons.M*Vars.v**5 + 4.0*0*Vars.v**6 + 4.0*4*Cons.E_2*Cons.M*Vars.v**7 + 4.0*5*Vars.E_SO_3*Cons.M*Vars.v**8) + -Cons.Fcal_3 - Vars.Fcal_SO_3*(4.0*2*Cons.E_0*Cons.M*Vars.v**6 + 4.0*0*Vars.v**7 + 4.0*4*Cons.E_2*Cons.M*Vars.v**8) + -Cons.Fcal_4 - Vars.Fcal_SQ_4*(4.0*2*Cons.E_0*Cons.M*Vars.v**7 + 4.0*0*Vars.v**8) + 4.0*-Cons.Fcal_5 - Vars.Fcal_SO_5*2*Cons.E_0*Cons.M*Vars.v**8) + -Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*(-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*(6.0*2*Cons.E_0*Cons.M*Vars.v**6 + 6.0*0*Vars.v**7 + 6.0*4*Cons.E_2*Cons.M*Vars.v**8) + -Cons.Fcal_3 - Vars.Fcal_SO_3*(12.0*2*Cons.E_0*Cons.M*Vars.v**7 + 12.0*0*Vars.v**8) + 12.0*-Cons.Fcal_4 - Vars.Fcal_SQ_4*2*Cons.E_0*Cons.M*Vars.v**8) + 6.0*-Cons.Fcal_3 - Vars.Fcal_SO_3**2*2*Cons.E_0*Cons.M*Vars.v**8) + -Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff**2*(-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*(4.0*2*Cons.E_0*Cons.M*Vars.v**7 + 4.0*0*Vars.v**8) + 12.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*2*Cons.E_0*Cons.M*Vars.v**8)) + 1.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff**4*2*Cons.E_0*Cons.M*Vars.v**8 + (0**2*(0*(0*(0*(-1.0*2*Cons.E_0*Cons.M*Vars.v**5 - 1.0*0*Vars.v**6 - 1.0*4*Cons.E_2*Cons.M*Vars.v**7 - 1.0*5*Vars.E_SO_3*Cons.M*Vars.v**8) + -Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*(-5.0*2*Cons.E_0*Cons.M*Vars.v**6 - 5.0*0*Vars.v**7 - 5.0*4*Cons.E_2*Cons.M*Vars.v**8) + -Cons.Fcal_3 - Vars.Fcal_SO_3*(-5.0*2*Cons.E_0*Cons.M*Vars.v**7 - 5.0*0*Vars.v**8) - 5.0*-Cons.Fcal_4 - Vars.Fcal_SQ_4*2*Cons.E_0*Cons.M*Vars.v**8) + -Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*(-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*(-10.0*2*Cons.E_0*Cons.M*Vars.v**7 - 10.0*0*Vars.v**8) - 20.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*2*Cons.E_0*Cons.M*Vars.v**8)) - 10.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff**3*2*Cons.E_0*Cons.M*Vars.v**8) + (0**4*(0*(0*(1.0*2*Cons.E_0*Cons.M*Vars.v**6 + 1.0*0*Vars.v**7 + 1.0*4*Cons.E_2*Cons.M*Vars.v**8) + -Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*(6.0*2*Cons.E_0*Cons.M*Vars.v**7 + 6.0*0*Vars.v**8) + 6.0*-Cons.Fcal_3 - Vars.Fcal_SO_3*2*Cons.E_0*Cons.M*Vars.v**8) + 15.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff**2*2*Cons.E_0*Cons.M*Vars.v**8) + (0**6*(0*(-1.0*2*Cons.E_0*Cons.M*Vars.v**7 - 1.0*0*Vars.v**8) - 7.0*-Cons.Fcal_2 - Cons.Fecc_2*Vars.Fecc_coeff*2*Cons.E_0*Cons.M*Vars.v**8) + 1.0*0**8*2*Cons.E_0*Cons.M*Vars.v**8/-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff)/-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff)/-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff)/-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff)/-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff)/-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff)/-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff)/-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff)/(Vars.Fcal_coeff*-Cons.Fcal_0 - Cons.Fecc_0*Vars.Fecc_coeff)
    dvdt_T5 = 1.0/dtdv
    dydt=np.zeros(10)
    [dydt[5],dydt[6],dydt[7],dydt[8]] = FrameFromAngularVelocityIntegrand(Vars.R, OmegaVec_4p0(Cons,Vars)[1:])
    dydt[0] = dvdt_T5[0]
    if(Cons.EvolveSpin1):
        dydt[1], dydt[2]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi1[0], Vars.rfrak_chi1[1],(mul(mul(inverse(Cons.S_chi1),OmegaVec_chiVec_1_4p0(Cons,Vars)),Cons.S_chi1))[1:])
    else:
        dydt[1], dydt[2] = 0.0, 0.0
    if(Cons.EvolveSpin2):
        dydt[3], dydt[4]=FrameFromAngularVelocity_2D_Integrand(Vars.rfrak_chi2[0], Vars.rfrak_chi2[1],(mul(mul(inverse(Cons.S_chi2),OmegaVec_chiVec_2_4p0(Cons,Vars)),Cons.S_chi2))[1:])
    else:
        dydt[3], dydt[4] = 0.0, 0.0
    dydt[9] = np.linalg.norm(OmegaVec_4p0(Cons,Vars)[1:])
    return dydt      

class PNEv:
    def Integrand(t,y):
        PNEv.Recalculate.get(2*PNEv.PNEvolutionOrder)(PNEv.Cons,PNEv.Vars,y)
        dydt=PNEv.Taylor.get(PNEv.TaylorTn+20*PNEv.PNEvolutionOrder)(PNEv.Cons,PNEv.Vars)
        if PNEv.Vars.v>=1.0 and PNEv.NotForward:
            print("Beyond domain of PN validity, this is a good way to terminate.")
            PNEv.terminal1=False
        #if dydt[0]<1.0e-12 and PNEv.NotForward:
        #    print("v is decreasing, which is not an uncommon way to stop.")
        #    PNEv.terminal2=False
        return dydt
        
    def Evolution(wHat_i, xHat_i, yHat_i, zHat_i, M1_i, M2_i, e_i, xi_i, v_i, S_chi1_i, S_chi2_i, R_i,
        t_PNStart=False, t_PNEnd=False, PNEvolutionOrder=3.5, TaylorTn=1, StepsPerOrbit=32, ForwardInTime=True, tol=1e-8, MinStep=1e-7): 
        # Initialization of constants
        PNEv.terminal1=True
        PNEv.terminal2=True
        PNEv.NotForward=True
        PNEv.PNEvolutionOrder=PNEvolutionOrder
        PNEv.TaylorTn=TaylorTn
        PNEv.Recalculate={            0:Recalculate_0,
            1:Recalculate_0p50,
            2:Recalculate_1p0,
            3:Recalculate_1p5,
            4:Recalculate_2p0,
            5:Recalculate_2p5,
            6:Recalculate_3p0,
            7:Recalculate_3p5,
            8:Recalculate_4p0}
        PNEv.Taylor={
            1:TaylorT1_0,
            11:TaylorT1_0p50,
            21:TaylorT1_1p0,
            31:TaylorT1_1p5,
            41:TaylorT1_2p0,
            51:TaylorT1_2p5,
            61:TaylorT1_3p0,
            71:TaylorT1_3p5,
            81:TaylorT1_4p0,
            4:TaylorT4_0,
            14:TaylorT4_0p50,
            24:TaylorT4_1p0,
            34:TaylorT4_1p5,
            44:TaylorT4_2p0,
            54:TaylorT4_2p5,
            64:TaylorT4_3p0,
            74:TaylorT4_3p5,
            84:TaylorT4_4p0,
            5:TaylorT5_0,
            15:TaylorT5_0p50,
            25:TaylorT5_1p0,
            35:TaylorT5_1p5,
            45:TaylorT5_2p0,
            55:TaylorT5_2p5,
            65:TaylorT5_3p0,
            75:TaylorT5_3p5,
            85:TaylorT5_4p0}
        z=np.array([0.0])
        PNEv.Cons=Cons(z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,True,True)
        PNEv.Vars=Vars(z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z)
        Initialization(PNEv.Cons, wHat_i, xHat_i, yHat_i, zHat_i, M1_i, M2_i, e_i, xi_i, v_i, v_i, S_chi1_i, S_chi2_i)
    
        def terminate(t,y):
            return 1.0*PNEv.terminal1*PNEv.terminal2
        terminate.terminal=True
        TMerger=5.0/(256.0*PNEv.Cons.nu*v_i**8)
        TEnd=TMerger
        if t_PNEnd:
            TEnd=t_PNEnd
        time=[0.0]
        while time[-1]<TEnd and 2*PNEv.Cons.M*(256*PNEv.Cons.nu*(TMerger-time[-1])/5)**(3/8)/StepsPerOrbit>MinStep:
            time.append(time[-1]+(2*PNEv.Cons.M*(256*PNEv.Cons.nu*(TMerger-time[-1])/5)**(3/8)/StepsPerOrbit)[0])
        time=np.delete(time, -1)
       
        # Integrate
        try:
            yy=solve_ivp(PNEv.Integrand, [time[0],time[-1]], [v_i,0.0,
                0.0,0.0,0.0,R_i[0],R_i[1],R_i[2],R_i[3],0.0], method='DOP853',
                t_eval=time, dense_output=True, events=terminate, rtol=tol, atol=tol)
        except:
            yy=solve_ivp(PNEv.Integrand, [time[0],time[-1]], [v_i,0.0,
                0.0,0.0,0.0,R_i[0],R_i[1],R_i[2],R_i[3],0.0], method='DOP853',
                dense_output=True, events=terminate, rtol=tol, atol=tol)
            time=time[time<yy.t[-1]]
            yy=solve_ivp(PNEv.Integrand, [time[0],time[-1]], [v_i,0.0,
                0.0,0.0,0.0,R_i[0],R_i[1],R_i[2],R_i[3],0.0], method='DOP853',
                t_eval=time, dense_output=True, events=terminate, rtol=tol, atol=tol)
        if ForwardInTime:
            PNEv.NotForward=False
            time=[0.0]
            TStart=-3*TMerger
            if t_PNStart:
                TStart=t_PNStart
            while time[-1]>TStart:
                time.append(time[-1]-(2*PNEv.Cons.M*(256*PNEv.Cons.nu*(TMerger-time[-1])/5)**(3/8)/StepsPerOrbit)[0])
            yyForward=solve_ivp(PNEv.Integrand, [time[0],time[-1]], [v_i,0.0,
                0.0,0.0,0.0,R_i[0],R_i[1],R_i[2],R_i[3],0.0], method='DOP853',
                t_eval=time, dense_output=True, rtol=tol, atol=tol)
            yy.t=np.append(yyForward.t[1:][::-1],yy.t)
            data=np.empty((10,len(yy.t)))
            for i in range(10):
                data[i]=np.append(yyForward.y[i][1:][::-1],yy.y[i])
            yy.y=data
             
        return yy
