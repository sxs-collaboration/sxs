# File produced automatically by PNCodeGen.ipynb
import numpy as np
from numpy import dot, log, sqrt, pi
from numpy import euler_gamma as EulerGamma
from numba import jit, njit, float64, complex128
from numba.experimental import jitclass
import quaternionic

qmul = njit(quaternionic.algebra.multiply)
qexp=njit(quaternionic.algebra.exp)
qconj=njit(quaternionic.algebra.conj)
I=1j

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
def normalized(A):
    return A/np.linalg.norm(A)

ConsSpec=[('wHat', float64[:]),('xHat', float64[:]),('yHat', float64[:]),('zHat', float64[:]),('M1', float64[:]),('M2', float64[:]),('S_chi1', float64[:]),('S_chi2', float64[:]),('M', float64[:]),('delta', float64[:]),('nu', float64[:]),('hHat_0_0_0', complex128[:]),('hHat_0_0_2', complex128[:]),('hHat_0_0_4', complex128[:]),('hHat_0_0_6', complex128[:]),('hHat_1_1_6', complex128[:]),('hHat_2_0_0', complex128[:]),('hHat_2_0_2', complex128[:]),('hHat_2_0_4', complex128[:]),('hHat_2_0_5', complex128[:]),('hHat_2_0_6', complex128[:]),('hHat_3_3_6', complex128[:]),('hHat_3_1_6', complex128[:]),('hHat_4_4_5', complex128[:]),('hHat_4_0_0', complex128[:]),('hHat_4_0_2', complex128[:]),('hHat_4_0_4', complex128[:]),('hHat_4_0_5', complex128[:]),('hHat_4_0_6', complex128[:]),('hHat_5_5_6', complex128[:]),('hHat_5_3_6', complex128[:]),('hHat_5_1_6', complex128[:]),('hHat_6_0_0', complex128[:]),('hHat_6_0_2', complex128[:]),('hHat_6_0_4', complex128[:]),('hHat_6_0_5', complex128[:]),('hHat_6_0_6', complex128[:]),('hHat_8_0_0', complex128[:]),('hHat_8_0_2', complex128[:]),('hHat_8_0_4', complex128[:]),('hHat_8_0_6', complex128[:])]
@jitclass(ConsSpec)
class Constants:
    def __init__(self,wHat,xHat,yHat,zHat,M1,M2,S_chi1,S_chi2,M,delta,nu,hHat_0_0_0,hHat_0_0_2,hHat_0_0_4,hHat_0_0_6,hHat_1_1_6,hHat_2_0_0,hHat_2_0_2,hHat_2_0_4,hHat_2_0_5,hHat_2_0_6,hHat_3_3_6,hHat_3_1_6,hHat_4_4_5,hHat_4_0_0,hHat_4_0_2,hHat_4_0_4,hHat_4_0_5,hHat_4_0_6,hHat_5_5_6,hHat_5_3_6,hHat_5_1_6,hHat_6_0_0,hHat_6_0_2,hHat_6_0_4,hHat_6_0_5,hHat_6_0_6,hHat_8_0_0,hHat_8_0_2,hHat_8_0_4,hHat_8_0_6):
        self.wHat=wHat
        self.xHat=xHat
        self.yHat=yHat
        self.zHat=zHat
        self.M1=M1
        self.M2=M2
        self.S_chi1=S_chi1
        self.S_chi2=S_chi2
        self.M=M
        self.delta=delta
        self.nu=nu
        self.hHat_0_0_0=hHat_0_0_0
        self.hHat_0_0_2=hHat_0_0_2
        self.hHat_0_0_4=hHat_0_0_4
        self.hHat_0_0_6=hHat_0_0_6
        self.hHat_1_1_6=hHat_1_1_6
        self.hHat_2_0_0=hHat_2_0_0
        self.hHat_2_0_2=hHat_2_0_2
        self.hHat_2_0_4=hHat_2_0_4
        self.hHat_2_0_5=hHat_2_0_5
        self.hHat_2_0_6=hHat_2_0_6
        self.hHat_3_3_6=hHat_3_3_6
        self.hHat_3_1_6=hHat_3_1_6
        self.hHat_4_4_5=hHat_4_4_5
        self.hHat_4_0_0=hHat_4_0_0
        self.hHat_4_0_2=hHat_4_0_2
        self.hHat_4_0_4=hHat_4_0_4
        self.hHat_4_0_5=hHat_4_0_5
        self.hHat_4_0_6=hHat_4_0_6
        self.hHat_5_5_6=hHat_5_5_6
        self.hHat_5_3_6=hHat_5_3_6
        self.hHat_5_1_6=hHat_5_1_6
        self.hHat_6_0_0=hHat_6_0_0
        self.hHat_6_0_2=hHat_6_0_2
        self.hHat_6_0_4=hHat_6_0_4
        self.hHat_6_0_5=hHat_6_0_5
        self.hHat_6_0_6=hHat_6_0_6
        self.hHat_8_0_0=hHat_8_0_0
        self.hHat_8_0_2=hHat_8_0_2
        self.hHat_8_0_4=hHat_8_0_4
        self.hHat_8_0_6=hHat_8_0_6

VarsSpec=[('v', float64[:]),('rfrak_chi1', float64[:]),('rfrak_chi2', float64[:]),('R', float64[:]),('R_S1', float64[:]),('R_S2', float64[:]),('chiVec1', float64[:]),('chiVec2', float64[:]),('chi1_n', float64[:]),('chi1_lambda', float64[:]),('chi1_ell', float64[:]),('chi2_n', float64[:]),('chi2_lambda', float64[:]),('chi2_ell', float64[:]),('S_ell', float64[:]),('S_n', float64[:]),('S_lambda', float64[:]),('Sigma_ell', float64[:]),('Sigma_n', float64[:]),('Sigma_lambda', float64[:]),('x', float64[:]),('rhOverM_coeff', float64[:]),('hHat_spin_Symm_0_0_3', complex128[:]),('hHat_spin_Symm_0_0_4', complex128[:]),('hHat_spin_Symm_2_0_3', complex128[:]),('hHat_spin_Symm_2_0_4', complex128[:]),('hHat_spin_Symm_2_0_5', complex128[:]),('hHat_spin_Symm_2_0_6', complex128[:]),('hHat_spin_Symm_2_0_7', complex128[:]),('hHat_spin_Symm_3_3_7', complex128[:]),('hHat_spin_Symm_3_1_7', complex128[:]),('hHat_spin_Symm_4_0_3', complex128[:]),('hHat_spin_Symm_4_0_4', complex128[:]),('hHat_spin_Symm_4_0_5', complex128[:]),('hHat_spin_Symm_4_0_6', complex128[:]),('hHat_spin_Symm_4_0_7', complex128[:]),('hHat_spin_Symm_6_0_5', complex128[:]),('hHat_spin_Symm_6_0_6', complex128[:]),('hHat_spin_Symm_6_0_7', complex128[:]),('hHat_spin_Symm_8_0_7', complex128[:]),('hHat_spin_Asymm_2_1_3', complex128[:]),('hHat_spin_Asymm_2_1_4', complex128[:]),('hHat_spin_Asymm_4_1_3', complex128[:]),('hHat_spin_Asymm_4_1_4', complex128[:])]
@jitclass(VarsSpec)
class Variables:
    def __init__(self,v,rfrak_chi1,rfrak_chi2,R,R_S1,R_S2,chiVec1,chiVec2,chi1_n,chi1_lambda,chi1_ell,chi2_n,chi2_lambda,chi2_ell,S_ell,S_n,S_lambda,Sigma_ell,Sigma_n,Sigma_lambda,x,rhOverM_coeff,hHat_spin_Symm_0_0_3,hHat_spin_Symm_0_0_4,hHat_spin_Symm_2_0_3,hHat_spin_Symm_2_0_4,hHat_spin_Symm_2_0_5,hHat_spin_Symm_2_0_6,hHat_spin_Symm_2_0_7,hHat_spin_Symm_3_3_7,hHat_spin_Symm_3_1_7,hHat_spin_Symm_4_0_3,hHat_spin_Symm_4_0_4,hHat_spin_Symm_4_0_5,hHat_spin_Symm_4_0_6,hHat_spin_Symm_4_0_7,hHat_spin_Symm_6_0_5,hHat_spin_Symm_6_0_6,hHat_spin_Symm_6_0_7,hHat_spin_Symm_8_0_7,hHat_spin_Asymm_2_1_3,hHat_spin_Asymm_2_1_4,hHat_spin_Asymm_4_1_3,hHat_spin_Asymm_4_1_4):
        self.v=v
        self.rfrak_chi1=rfrak_chi1
        self.rfrak_chi2=rfrak_chi2
        self.R=R
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
        self.x=x
        self.rhOverM_coeff=rhOverM_coeff
        self.hHat_spin_Symm_0_0_3=hHat_spin_Symm_0_0_3
        self.hHat_spin_Symm_0_0_4=hHat_spin_Symm_0_0_4
        self.hHat_spin_Symm_2_0_3=hHat_spin_Symm_2_0_3
        self.hHat_spin_Symm_2_0_4=hHat_spin_Symm_2_0_4
        self.hHat_spin_Symm_2_0_5=hHat_spin_Symm_2_0_5
        self.hHat_spin_Symm_2_0_6=hHat_spin_Symm_2_0_6
        self.hHat_spin_Symm_2_0_7=hHat_spin_Symm_2_0_7
        self.hHat_spin_Symm_3_3_7=hHat_spin_Symm_3_3_7
        self.hHat_spin_Symm_3_1_7=hHat_spin_Symm_3_1_7
        self.hHat_spin_Symm_4_0_3=hHat_spin_Symm_4_0_3
        self.hHat_spin_Symm_4_0_4=hHat_spin_Symm_4_0_4
        self.hHat_spin_Symm_4_0_5=hHat_spin_Symm_4_0_5
        self.hHat_spin_Symm_4_0_6=hHat_spin_Symm_4_0_6
        self.hHat_spin_Symm_4_0_7=hHat_spin_Symm_4_0_7
        self.hHat_spin_Symm_6_0_5=hHat_spin_Symm_6_0_5
        self.hHat_spin_Symm_6_0_6=hHat_spin_Symm_6_0_6
        self.hHat_spin_Symm_6_0_7=hHat_spin_Symm_6_0_7
        self.hHat_spin_Symm_8_0_7=hHat_spin_Symm_8_0_7
        self.hHat_spin_Asymm_2_1_3=hHat_spin_Asymm_2_1_3
        self.hHat_spin_Asymm_2_1_4=hHat_spin_Asymm_2_1_4
        self.hHat_spin_Asymm_4_1_3=hHat_spin_Asymm_4_1_3
        self.hHat_spin_Asymm_4_1_4=hHat_spin_Asymm_4_1_4

@njit(cache=True)
def Initialization(Cons, wHat_i, xHat_i, yHat_i, zHat_i, M1_i, M2_i, v_i, S_chi1_i, S_chi2_i): 
    Cons.wHat=wHat_i
    Cons.xHat=xHat_i
    Cons.yHat=yHat_i
    Cons.zHat=zHat_i
    Cons.M1=np.array([M1_i])
    Cons.M2=np.array([M2_i])
    Cons.S_chi1=S_chi1_i
    Cons.S_chi2=S_chi2_i
    Cons.M=Cons.M1 + Cons.M2
    Cons.delta=(Cons.M1 - Cons.M2)/Cons.M
    Cons.nu=Cons.M1*Cons.M2/Cons.M**2
    Cons.hHat_0_0_0=np.array([1 + 0.0*I])
    Cons.hHat_0_0_2=-Cons.nu/12 - 0.75 + 0.0*I
    Cons.hHat_0_0_4=-Cons.nu**2/24 + 19*Cons.nu/8 - 3.375 + 0.0*I
    Cons.hHat_0_0_6=-35*Cons.nu**3/5184 - 155*Cons.nu**2/96 + Cons.nu*(59.8003472222222 - 205*pi**2/96) - 10.546875 + 0.0*I
    Cons.hHat_1_1_6=10.8244118157276*Cons.delta*Cons.nu + 0.0*I
    Cons.hHat_2_0_0=np.array([2*sqrt(5)/7 + 0.0*I])
    Cons.hHat_2_0_2=335*Cons.nu/168 - 4075*sqrt(5)/14112 + 0.0*I
    Cons.hHat_2_0_4=205*sqrt(5)*Cons.nu**2/1232 - 123815*sqrt(5)*Cons.nu/155232 - 151877213*sqrt(5)/234710784 + 0.0*I
    Cons.hHat_2_0_5=0.860544217687075*sqrt(5)*Cons.nu*pi - 253*sqrt(5)*pi/1176 + 0.0*I
    Cons.hHat_2_0_6=1321981*sqrt(5)*Cons.nu**3/20756736 + 69527951*sqrt(5)*Cons.nu**2/581188608 + Cons.nu*(-205*sqrt(5)*pi**2/336 + 700464542023*sqrt(5)/48819843072) - 4397711103307*sqrt(5)/1864030371840 + 0.0*I
    Cons.hHat_3_3_6=-44*sqrt(35)*Cons.delta*Cons.nu/945 + 0.0*I
    Cons.hHat_3_1_6=484*sqrt(21)*Cons.delta*Cons.nu/315 + 0.0*I
    Cons.hHat_4_4_5=-0.318727629155838*I*Cons.nu
    Cons.hHat_4_0_0=np.array([0.023809523809523808 + 0.0*I])
    Cons.hHat_4_0_2=27227*Cons.nu/44352 - 0.14502567125335 + 0.0*I
    Cons.hHat_4_0_4=844951*Cons.nu**2/1153152 - 34829479*Cons.nu/18162144 + 0.330678707853645 + 0.0*I
    Cons.hHat_4_0_5=13565*Cons.nu*pi/12936 - 13565*pi/51744 + 0.0*I
    Cons.hHat_4_0_6=221405645*Cons.nu**3/498161664 - 4174614175*Cons.nu**2/1549836288 + Cons.nu*(-0.878864423125315 - 205*pi**2/4032) + 0.464550058346122 + 0.0*I
    Cons.hHat_5_5_6=-36*sqrt(77)*Cons.delta*Cons.nu/385 + 0.0*I
    Cons.hHat_5_3_6=4*sqrt(385)*Cons.delta*Cons.nu/10395 + 0.0*I
    Cons.hHat_5_1_6=0.27261959898367*Cons.delta*Cons.nu + 0.0*I
    Cons.hHat_6_0_0=np.array([0 + 0.0*I])
    Cons.hHat_6_0_2=215*sqrt(13)*Cons.nu/27456 - 4195*sqrt(13)/2306304 + 0.0*I
    Cons.hHat_6_0_4=17185*sqrt(13)*Cons.nu**2/164736 - 253535*sqrt(13)*Cons.nu/3459456 + 4151051*sqrt(13)/317011968 + 0.0*I
    Cons.hHat_6_0_5=0.0108225108225108*sqrt(13)*Cons.nu*pi - 5*sqrt(13)*pi/1848 + 0.0*I
    Cons.hHat_6_0_6=17709755*sqrt(13)*Cons.nu**3/134424576 - 941405305*sqrt(13)*Cons.nu**2/1613094912 + 27653500031*sqrt(13)*Cons.nu/105388867584 - 3012132889099*sqrt(13)/79673983893504 + 0.0*I
    Cons.hHat_8_0_0=np.array([0 + 0.0*I])
    Cons.hHat_8_0_2=np.array([0 + 0.0*I])
    Cons.hHat_8_0_4=3395*sqrt(17)*Cons.nu**2/700128 - 25115*sqrt(17)*Cons.nu/8401536 + 75601*sqrt(17)/151227648 + 0.0*I
    Cons.hHat_8_0_6=9677185*sqrt(17)*Cons.nu**3/212838912 - 80280125*sqrt(17)*Cons.nu**2/1277033472 + 18177898147*sqrt(17)*Cons.nu/643624869888 - 265361599*sqrt(17)/67749986304 + 0.0*I

@njit(cache=True)
def Recalculate_0(Cons,Vars,y):
    Vars.v = np.array([y[0]])
    Vars.rfrak_chi1 = np.array([y[1],y[2]])
    Vars.rfrak_chi2 = np.array([y[3],y[4]])
    Vars.R_S1 = exp(Vars.rfrak_chi1[0]*Cons.xHat + Vars.rfrak_chi1[1]*Cons.yHat)
    Vars.R_S2 = exp(Vars.rfrak_chi2[0]*Cons.xHat + Vars.rfrak_chi2[1]*Cons.yHat)
    Vars.rhOverM_coeff = Cons.nu*Vars.v**2*sqrt(pi)

@njit(cache=True)
def Modes_0(modes,Cons,Vars):
    # (ell, m) = (0, +/- 0)
    Symm = Cons.hHat_0_0_0*Vars.rhOverM_coeff
    Asymm = 0
    modes[0] = (Symm + Asymm)[0]
    # (ell, m) = (2, +/- 0)
    Symm = Cons.hHat_2_0_0*Vars.rhOverM_coeff
    Asymm = 0
    modes[6] = (Symm + Asymm)[0]
    # (ell, m) = (4, +/- 0)
    Symm = Cons.hHat_4_0_0*Vars.rhOverM_coeff
    Asymm = 0
    modes[20] = (Symm + Asymm)[0]
    # (ell, m) = (6, +/- 0)
    Symm = Cons.hHat_6_0_0*Vars.rhOverM_coeff
    Asymm = 0
    modes[42] = (Symm + Asymm)[0]
    # (ell, m) = (8, +/- 0)
    Symm = Cons.hHat_8_0_0*Vars.rhOverM_coeff
    Asymm = 0
    modes[72] = (Symm + Asymm)[0]

@njit(cache=True)
def Recalculate_0p50(Cons,Vars,y):
    Vars.v = np.array([y[0]])
    Vars.rfrak_chi1 = np.array([y[1],y[2]])
    Vars.rfrak_chi2 = np.array([y[3],y[4]])
    Vars.R_S1 = exp(Vars.rfrak_chi1[0]*Cons.xHat + Vars.rfrak_chi1[1]*Cons.yHat)
    Vars.R_S2 = exp(Vars.rfrak_chi2[0]*Cons.xHat + Vars.rfrak_chi2[1]*Cons.yHat)
    Vars.rhOverM_coeff = Cons.nu*Vars.v**2*sqrt(pi)

@njit(cache=True)
def Modes_0p50(modes,Cons,Vars):
    # (ell, m) = (0, +/- 0)
    Symm = Cons.hHat_0_0_0*Vars.rhOverM_coeff
    Asymm = 0
    modes[0] = (Symm + Asymm)[0]
    # (ell, m) = (2, +/- 0)
    Symm = Cons.hHat_2_0_0*Vars.rhOverM_coeff
    Asymm = 0
    modes[6] = (Symm + Asymm)[0]
    # (ell, m) = (4, +/- 0)
    Symm = Cons.hHat_4_0_0*Vars.rhOverM_coeff
    Asymm = 0
    modes[20] = (Symm + Asymm)[0]
    # (ell, m) = (6, +/- 0)
    Symm = Cons.hHat_6_0_0*Vars.rhOverM_coeff
    Asymm = 0
    modes[42] = (Symm + Asymm)[0]
    # (ell, m) = (8, +/- 0)
    Symm = Cons.hHat_8_0_0*Vars.rhOverM_coeff
    Asymm = 0
    modes[72] = (Symm + Asymm)[0]

@njit(cache=True)
def Recalculate_1p0(Cons,Vars,y):
    Vars.v = np.array([y[0]])
    Vars.rfrak_chi1 = np.array([y[1],y[2]])
    Vars.rfrak_chi2 = np.array([y[3],y[4]])
    Vars.R_S1 = exp(Vars.rfrak_chi1[0]*Cons.xHat + Vars.rfrak_chi1[1]*Cons.yHat)
    Vars.R_S2 = exp(Vars.rfrak_chi2[0]*Cons.xHat + Vars.rfrak_chi2[1]*Cons.yHat)
    Vars.rhOverM_coeff = Cons.nu*Vars.v**2*sqrt(pi)

@njit(cache=True)
def Modes_1p0(modes,Cons,Vars):
    # (ell, m) = (0, +/- 0)
    Symm = Vars.rhOverM_coeff*(Cons.hHat_0_0_0 + Cons.hHat_0_0_2*Vars.v**2)
    Asymm = 0
    modes[0] = (Symm + Asymm)[0]
    # (ell, m) = (2, +/- 0)
    Symm = Vars.rhOverM_coeff*(Cons.hHat_2_0_0 + Cons.hHat_2_0_2*Vars.v**2)
    Asymm = 0
    modes[6] = (Symm + Asymm)[0]
    # (ell, m) = (4, +/- 0)
    Symm = Vars.rhOverM_coeff*(Cons.hHat_4_0_0 + Cons.hHat_4_0_2*Vars.v**2)
    Asymm = 0
    modes[20] = (Symm + Asymm)[0]
    # (ell, m) = (6, +/- 0)
    Symm = Vars.rhOverM_coeff*(Cons.hHat_6_0_0 + Cons.hHat_6_0_2*Vars.v**2)
    Asymm = 0
    modes[42] = (Symm + Asymm)[0]
    # (ell, m) = (8, +/- 0)
    Symm = Vars.rhOverM_coeff*(Cons.hHat_8_0_0 + Cons.hHat_8_0_2*Vars.v**2)
    Asymm = 0
    modes[72] = (Symm + Asymm)[0]

@njit(cache=True)
def Recalculate_1p5(Cons,Vars,y):
    Vars.v = np.array([y[0]])
    Vars.rfrak_chi1 = np.array([y[1],y[2]])
    Vars.rfrak_chi2 = np.array([y[3],y[4]])
    Vars.R = y[5]*Cons.wHat + y[6]*Cons.xHat + y[7]*Cons.yHat + y[8]*Cons.zHat
    Vars.R_S1 = exp(Vars.rfrak_chi1[0]*Cons.xHat + Vars.rfrak_chi1[1]*Cons.yHat)
    Vars.R_S2 = exp(Vars.rfrak_chi2[0]*Cons.xHat + Vars.rfrak_chi2[1]*Cons.yHat)
    Vars.chiVec1 = mul(mul(conjugate(normalized(y[5]*Cons.wHat + y[6]*Cons.xHat + y[7]*Cons.yHat + y[8]*Cons.zHat)),mul(mul(mul(Cons.S_chi1,Vars.R_S1),Cons.zHat),mul(conjugate(Vars.R_S1),conjugate(Cons.S_chi1)))),normalized(y[5]*Cons.wHat + y[6]*Cons.xHat + y[7]*Cons.yHat + y[8]*Cons.zHat))
    Vars.chiVec2 = mul(mul(conjugate(normalized(y[5]*Cons.wHat + y[6]*Cons.xHat + y[7]*Cons.yHat + y[8]*Cons.zHat)),mul(mul(mul(Cons.S_chi2,Vars.R_S2),Cons.zHat),mul(conjugate(Vars.R_S2),conjugate(Cons.S_chi2)))),normalized(y[5]*Cons.wHat + y[6]*Cons.xHat + y[7]*Cons.yHat + y[8]*Cons.zHat))
    Vars.chi1_n = np.array([Vars.chiVec1[1]])
    Vars.chi1_lambda = np.array([Vars.chiVec1[2]])
    Vars.chi1_ell = np.array([Vars.chiVec1[3]])
    Vars.chi2_n = np.array([Vars.chiVec2[1]])
    Vars.chi2_lambda = np.array([Vars.chiVec2[2]])
    Vars.chi2_ell = np.array([Vars.chiVec2[3]])
    Vars.S_ell = Cons.M1**2*Vars.chi1_ell + Cons.M2**2*Vars.chi2_ell
    Vars.S_n = Cons.M1**2*Vars.chi1_n + Cons.M2**2*Vars.chi2_n
    Vars.S_lambda = Cons.M1**2*Vars.chi1_lambda + Cons.M2**2*Vars.chi2_lambda
    Vars.Sigma_ell = Cons.M*(-Cons.M1*Vars.chi1_ell + Cons.M2*Vars.chi2_ell)
    Vars.Sigma_n = Cons.M*(-Cons.M1*Vars.chi1_n + Cons.M2*Vars.chi2_n)
    Vars.Sigma_lambda = Cons.M*(-Cons.M1*Vars.chi1_lambda + Cons.M2*Vars.chi2_lambda)
    Vars.rhOverM_coeff = Cons.nu*Vars.v**2*sqrt(pi)
    Vars.hHat_spin_Symm_0_0_3 = 14*Vars.S_ell/(3*Cons.M**2) + 2*Vars.Sigma_ell*Cons.delta/Cons.M**2 + 0.0*I
    Vars.hHat_spin_Symm_2_0_3 = 32*sqrt(5)*Vars.S_ell/(21*Cons.M**2) + 419*sqrt(5)*Vars.Sigma_ell*Cons.delta/(560*Cons.M**2) + 0.0*I
    Vars.hHat_spin_Symm_4_0_3 = 5*Vars.S_ell/(21*Cons.M**2) + 19*Vars.Sigma_ell*Cons.delta/(112*Cons.M**2) + 0.0*I
    Vars.hHat_spin_Asymm_2_1_3 = -0.611297497215587*Cons.delta*(-I*Vars.Sigma_lambda + Vars.Sigma_n)/Cons.M**2 - sqrt(30)*(-61*I*Vars.S_lambda + 61*Vars.S_n)/(420*Cons.M**2)
    Vars.hHat_spin_Asymm_4_1_3 = -17*sqrt(5)*Cons.delta*(-I*Vars.Sigma_lambda + Vars.Sigma_n)/(420*Cons.M**2) - sqrt(5)*(-13*I*Vars.S_lambda + 13*Vars.S_n)/(280*Cons.M**2)

@njit(cache=True)
def Modes_1p5(modes,Cons,Vars):
    # (ell, m) = (0, +/- 0)
    Symm = Vars.rhOverM_coeff*(Cons.hHat_0_0_0 + Vars.v**2*(Cons.hHat_0_0_2 + Vars.hHat_spin_Symm_0_0_3*Vars.v))
    Asymm = 0
    modes[0] = (Symm + Asymm)[0]
    # (ell, m) = (2, +/- 0)
    Symm = Vars.rhOverM_coeff*(Cons.hHat_2_0_0 + Vars.v**2*(Cons.hHat_2_0_2 + Vars.hHat_spin_Symm_2_0_3*Vars.v))
    Asymm = 0
    modes[6] = (Symm + Asymm)[0]
    # (ell, m) = (2, +/- 1)
    Symm = 0
    Asymm = Vars.hHat_spin_Asymm_2_1_3*Vars.rhOverM_coeff*Vars.v**3
    modes[7] = (Symm + Asymm)[0]
    modes[5] = (np.conjugate(Symm - Asymm))[0]
    # (ell, m) = (4, +/- 0)
    Symm = Vars.rhOverM_coeff*(Cons.hHat_4_0_0 + Vars.v**2*(Cons.hHat_4_0_2 + Vars.hHat_spin_Symm_4_0_3*Vars.v))
    Asymm = 0
    modes[20] = (Symm + Asymm)[0]
    # (ell, m) = (4, +/- 1)
    Symm = 0
    Asymm = Vars.hHat_spin_Asymm_4_1_3*Vars.rhOverM_coeff*Vars.v**3
    modes[21] = (Symm + Asymm)[0]
    modes[19] = (np.conjugate(Symm - Asymm))[0]
    # (ell, m) = (6, +/- 0)
    Symm = Vars.rhOverM_coeff*(Cons.hHat_6_0_0 + Cons.hHat_6_0_2*Vars.v**2)
    Asymm = 0
    modes[42] = (Symm + Asymm)[0]
    # (ell, m) = (8, +/- 0)
    Symm = Vars.rhOverM_coeff*(Cons.hHat_8_0_0 + Cons.hHat_8_0_2*Vars.v**2)
    Asymm = 0
    modes[72] = (Symm + Asymm)[0]

@njit(cache=True)
def Recalculate_2p0(Cons,Vars,y):
    Vars.v = np.array([y[0]])
    Vars.rfrak_chi1 = np.array([y[1],y[2]])
    Vars.rfrak_chi2 = np.array([y[3],y[4]])
    Vars.R = y[5]*Cons.wHat + y[6]*Cons.xHat + y[7]*Cons.yHat + y[8]*Cons.zHat
    Vars.R_S1 = exp(Vars.rfrak_chi1[0]*Cons.xHat + Vars.rfrak_chi1[1]*Cons.yHat)
    Vars.R_S2 = exp(Vars.rfrak_chi2[0]*Cons.xHat + Vars.rfrak_chi2[1]*Cons.yHat)
    Vars.chiVec1 = mul(mul(conjugate(normalized(y[5]*Cons.wHat + y[6]*Cons.xHat + y[7]*Cons.yHat + y[8]*Cons.zHat)),mul(mul(mul(Cons.S_chi1,Vars.R_S1),Cons.zHat),mul(conjugate(Vars.R_S1),conjugate(Cons.S_chi1)))),normalized(y[5]*Cons.wHat + y[6]*Cons.xHat + y[7]*Cons.yHat + y[8]*Cons.zHat))
    Vars.chiVec2 = mul(mul(conjugate(normalized(y[5]*Cons.wHat + y[6]*Cons.xHat + y[7]*Cons.yHat + y[8]*Cons.zHat)),mul(mul(mul(Cons.S_chi2,Vars.R_S2),Cons.zHat),mul(conjugate(Vars.R_S2),conjugate(Cons.S_chi2)))),normalized(y[5]*Cons.wHat + y[6]*Cons.xHat + y[7]*Cons.yHat + y[8]*Cons.zHat))
    Vars.chi1_n = np.array([Vars.chiVec1[1]])
    Vars.chi1_lambda = np.array([Vars.chiVec1[2]])
    Vars.chi1_ell = np.array([Vars.chiVec1[3]])
    Vars.chi2_n = np.array([Vars.chiVec2[1]])
    Vars.chi2_lambda = np.array([Vars.chiVec2[2]])
    Vars.chi2_ell = np.array([Vars.chiVec2[3]])
    Vars.S_ell = Cons.M1**2*Vars.chi1_ell + Cons.M2**2*Vars.chi2_ell
    Vars.S_n = Cons.M1**2*Vars.chi1_n + Cons.M2**2*Vars.chi2_n
    Vars.S_lambda = Cons.M1**2*Vars.chi1_lambda + Cons.M2**2*Vars.chi2_lambda
    Vars.Sigma_ell = Cons.M*(-Cons.M1*Vars.chi1_ell + Cons.M2*Vars.chi2_ell)
    Vars.Sigma_n = Cons.M*(-Cons.M1*Vars.chi1_n + Cons.M2*Vars.chi2_n)
    Vars.Sigma_lambda = Cons.M*(-Cons.M1*Vars.chi1_lambda + Cons.M2*Vars.chi2_lambda)
    Vars.rhOverM_coeff = Cons.nu*Vars.v**2*sqrt(pi)
    Vars.hHat_spin_Symm_0_0_3 = 14*Vars.S_ell/(3*Cons.M**2) + 2*Vars.Sigma_ell*Cons.delta/Cons.M**2 + 0.0*I
    Vars.hHat_spin_Symm_0_0_4 = -4*Cons.delta*(3*Vars.S_ell*Vars.Sigma_ell + Vars.S_lambda*Vars.Sigma_n + Vars.S_n*Vars.Sigma_n)/(3*Cons.M**4) + Cons.nu*(12*Vars.Sigma_ell**2 + 4*Vars.Sigma_lambda**2 + 4*Vars.Sigma_n**2)/(3*Cons.M**4) + (-48*Vars.S_ell**2 - 16*Vars.S_lambda**2 - 16*Vars.S_n**2 - 12*Vars.Sigma_ell**2 - 3*Vars.Sigma_lambda**2 - 3*Vars.Sigma_n**2)/(12*Cons.M**4) + 0.0*I
    Vars.hHat_spin_Symm_2_0_3 = 32*sqrt(5)*Vars.S_ell/(21*Cons.M**2) + 419*sqrt(5)*Vars.Sigma_ell*Cons.delta/(560*Cons.M**2) + 0.0*I
    Vars.hHat_spin_Symm_2_0_4 = -8*sqrt(5)*Cons.delta*(3*Vars.S_ell*Vars.Sigma_ell + Vars.S_lambda*Vars.Sigma_lambda + Vars.S_n*Vars.Sigma_n)/(21*Cons.M**4) + 8*sqrt(5)*Cons.nu*(3*Vars.Sigma_ell**2 + Vars.Sigma_lambda**2 + Vars.Sigma_n**2)/(21*Cons.M**4) - sqrt(5)*(384*Vars.S_ell**2 + 128*Vars.S_lambda**2 + 128*Vars.S_n**2 + 99*Vars.Sigma_ell**2 + 24*Vars.Sigma_lambda**2 + 24*Vars.Sigma_n**2)/(336*Cons.M**4) + 0.0*I
    Vars.hHat_spin_Symm_4_0_3 = 5*Vars.S_ell/(21*Cons.M**2) + 19*Vars.Sigma_ell*Cons.delta/(112*Cons.M**2) + 0.0*I
    Vars.hHat_spin_Symm_4_0_4 = -2*Cons.delta*(3*Vars.S_ell*Vars.Sigma_ell + Vars.S_lambda*Vars.Sigma_lambda + Vars.S_n*Vars.Sigma_n)/(63*Cons.M**4) + Cons.nu*(6*Vars.Sigma_ell**2 + 2*Vars.Sigma_lambda**2 + 2*Vars.Sigma_n**2)/(63*Cons.M**4) + (-192*Vars.S_ell**2 - 64*Vars.S_lambda**2 - 64*Vars.S_n**2 - 53*Vars.Sigma_ell**2 - 12*Vars.Sigma_lambda**2 - 12*Vars.Sigma_n**2)/(2016*Cons.M**4) + 0.0*I
    Vars.hHat_spin_Asymm_2_1_3 = -0.611297497215587*Cons.delta*(-I*Vars.Sigma_lambda + Vars.Sigma_n)/Cons.M**2 - sqrt(30)*(-61*I*Vars.S_lambda + 61*Vars.S_n)/(420*Cons.M**2)
    Vars.hHat_spin_Asymm_2_1_4 = -0.391230398217976*Vars.Sigma_ell*Cons.nu*(-I*Vars.Sigma_lambda + Vars.Sigma_n)/Cons.M**4 + 0.195615199108988*Cons.delta*(Vars.S_ell*(-I*Vars.Sigma_lambda + Vars.Sigma_n) + Vars.Sigma_ell*(-I*Vars.S_lambda + Vars.S_n))/Cons.M**4 + (2.73861278752583*Vars.S_ell*(-I*Vars.S_lambda + Vars.S_n) + 0.912870929175277*Vars.Sigma_ell*(-I*Vars.Sigma_lambda + Vars.Sigma_n))/(7*Cons.M**4)
    Vars.hHat_spin_Asymm_4_1_3 = -17*sqrt(5)*Cons.delta*(-I*Vars.Sigma_lambda + Vars.Sigma_n)/(420*Cons.M**2) - sqrt(5)*(-13*I*Vars.S_lambda + 13*Vars.S_n)/(280*Cons.M**2)
    Vars.hHat_spin_Asymm_4_1_4 = -sqrt(5)*Vars.Sigma_ell*Cons.nu*(-I*Vars.Sigma_lambda + Vars.Sigma_n)/(84*Cons.M**4) + sqrt(5)*Cons.delta*(Vars.S_ell*(-I*Vars.Sigma_lambda + Vars.Sigma_n) + Vars.Sigma_ell*(-I*Vars.S_lambda + Vars.S_n))/(168*Cons.M**4) + sqrt(5)*(3*Vars.S_ell*(-I*Vars.S_lambda + Vars.S_n) + Vars.Sigma_ell*(-I*Vars.Sigma_lambda + Vars.Sigma_n))/(252*Cons.M**4)

@njit(cache=True)
def Modes_2p0(modes,Cons,Vars):
    # (ell, m) = (0, +/- 0)
    Symm = Vars.rhOverM_coeff*(Cons.hHat_0_0_0 + Vars.v**2*(Cons.hHat_0_0_2 + Vars.v*(Vars.hHat_spin_Symm_0_0_3 + Vars.v*(Cons.hHat_0_0_4 + Vars.hHat_spin_Symm_0_0_4))))
    Asymm = 0
    modes[0] = (Symm + Asymm)[0]
    # (ell, m) = (2, +/- 0)
    Symm = Vars.rhOverM_coeff*(Cons.hHat_2_0_0 + Vars.v**2*(Cons.hHat_2_0_2 + Vars.v*(Vars.hHat_spin_Symm_2_0_3 + Vars.v*(Cons.hHat_2_0_4 + Vars.hHat_spin_Symm_2_0_4))))
    Asymm = 0
    modes[6] = (Symm + Asymm)[0]
    # (ell, m) = (2, +/- 1)
    Symm = 0
    Asymm = Vars.rhOverM_coeff*Vars.v**3*(Vars.hHat_spin_Asymm_2_1_3 + Vars.hHat_spin_Asymm_2_1_4*Vars.v)
    modes[7] = (Symm + Asymm)[0]
    modes[5] = (np.conjugate(Symm - Asymm))[0]
    # (ell, m) = (4, +/- 0)
    Symm = Vars.rhOverM_coeff*(Cons.hHat_4_0_0 + Vars.v**2*(Cons.hHat_4_0_2 + Vars.v*(Vars.hHat_spin_Symm_4_0_3 + Vars.v*(Cons.hHat_4_0_4 + Vars.hHat_spin_Symm_4_0_4))))
    Asymm = 0
    modes[20] = (Symm + Asymm)[0]
    # (ell, m) = (4, +/- 1)
    Symm = 0
    Asymm = Vars.rhOverM_coeff*Vars.v**3*(Vars.hHat_spin_Asymm_4_1_3 + Vars.hHat_spin_Asymm_4_1_4*Vars.v)
    modes[21] = (Symm + Asymm)[0]
    modes[19] = (np.conjugate(Symm - Asymm))[0]
    # (ell, m) = (6, +/- 0)
    Symm = Vars.rhOverM_coeff*(Cons.hHat_6_0_0 + Vars.v**2*(Cons.hHat_6_0_2 + Cons.hHat_6_0_4*Vars.v**2))
    Asymm = 0
    modes[42] = (Symm + Asymm)[0]
    # (ell, m) = (8, +/- 0)
    Symm = Vars.rhOverM_coeff*(Cons.hHat_8_0_0 + Vars.v**2*(Cons.hHat_8_0_2 + Cons.hHat_8_0_4*Vars.v**2))
    Asymm = 0
    modes[72] = (Symm + Asymm)[0]

@njit(cache=True)
def Recalculate_2p5(Cons,Vars,y):
    Vars.v = np.array([y[0]])
    Vars.rfrak_chi1 = np.array([y[1],y[2]])
    Vars.rfrak_chi2 = np.array([y[3],y[4]])
    Vars.R = y[5]*Cons.wHat + y[6]*Cons.xHat + y[7]*Cons.yHat + y[8]*Cons.zHat
    Vars.R_S1 = exp(Vars.rfrak_chi1[0]*Cons.xHat + Vars.rfrak_chi1[1]*Cons.yHat)
    Vars.R_S2 = exp(Vars.rfrak_chi2[0]*Cons.xHat + Vars.rfrak_chi2[1]*Cons.yHat)
    Vars.chiVec1 = mul(mul(conjugate(normalized(y[5]*Cons.wHat + y[6]*Cons.xHat + y[7]*Cons.yHat + y[8]*Cons.zHat)),mul(mul(mul(Cons.S_chi1,Vars.R_S1),Cons.zHat),mul(conjugate(Vars.R_S1),conjugate(Cons.S_chi1)))),normalized(y[5]*Cons.wHat + y[6]*Cons.xHat + y[7]*Cons.yHat + y[8]*Cons.zHat))
    Vars.chiVec2 = mul(mul(conjugate(normalized(y[5]*Cons.wHat + y[6]*Cons.xHat + y[7]*Cons.yHat + y[8]*Cons.zHat)),mul(mul(mul(Cons.S_chi2,Vars.R_S2),Cons.zHat),mul(conjugate(Vars.R_S2),conjugate(Cons.S_chi2)))),normalized(y[5]*Cons.wHat + y[6]*Cons.xHat + y[7]*Cons.yHat + y[8]*Cons.zHat))
    Vars.chi1_n = np.array([Vars.chiVec1[1]])
    Vars.chi1_lambda = np.array([Vars.chiVec1[2]])
    Vars.chi1_ell = np.array([Vars.chiVec1[3]])
    Vars.chi2_n = np.array([Vars.chiVec2[1]])
    Vars.chi2_lambda = np.array([Vars.chiVec2[2]])
    Vars.chi2_ell = np.array([Vars.chiVec2[3]])
    Vars.S_ell = Cons.M1**2*Vars.chi1_ell + Cons.M2**2*Vars.chi2_ell
    Vars.S_n = Cons.M1**2*Vars.chi1_n + Cons.M2**2*Vars.chi2_n
    Vars.S_lambda = Cons.M1**2*Vars.chi1_lambda + Cons.M2**2*Vars.chi2_lambda
    Vars.Sigma_ell = Cons.M*(-Cons.M1*Vars.chi1_ell + Cons.M2*Vars.chi2_ell)
    Vars.Sigma_n = Cons.M*(-Cons.M1*Vars.chi1_n + Cons.M2*Vars.chi2_n)
    Vars.Sigma_lambda = Cons.M*(-Cons.M1*Vars.chi1_lambda + Cons.M2*Vars.chi2_lambda)
    Vars.rhOverM_coeff = Cons.nu*Vars.v**2*sqrt(pi)
    Vars.hHat_spin_Symm_0_0_3 = 14*Vars.S_ell/(3*Cons.M**2) + 2*Vars.Sigma_ell*Cons.delta/Cons.M**2 + 0.0*I
    Vars.hHat_spin_Symm_0_0_4 = -4*Cons.delta*(3*Vars.S_ell*Vars.Sigma_ell + Vars.S_lambda*Vars.Sigma_n + Vars.S_n*Vars.Sigma_n)/(3*Cons.M**4) + Cons.nu*(12*Vars.Sigma_ell**2 + 4*Vars.Sigma_lambda**2 + 4*Vars.Sigma_n**2)/(3*Cons.M**4) + (-48*Vars.S_ell**2 - 16*Vars.S_lambda**2 - 16*Vars.S_n**2 - 12*Vars.Sigma_ell**2 - 3*Vars.Sigma_lambda**2 - 3*Vars.Sigma_n**2)/(12*Cons.M**4) + 0.0*I
    Vars.hHat_spin_Symm_2_0_3 = 32*sqrt(5)*Vars.S_ell/(21*Cons.M**2) + 419*sqrt(5)*Vars.Sigma_ell*Cons.delta/(560*Cons.M**2) + 0.0*I
    Vars.hHat_spin_Symm_2_0_4 = -8*sqrt(5)*Cons.delta*(3*Vars.S_ell*Vars.Sigma_ell + Vars.S_lambda*Vars.Sigma_lambda + Vars.S_n*Vars.Sigma_n)/(21*Cons.M**4) + 8*sqrt(5)*Cons.nu*(3*Vars.Sigma_ell**2 + Vars.Sigma_lambda**2 + Vars.Sigma_n**2)/(21*Cons.M**4) - sqrt(5)*(384*Vars.S_ell**2 + 128*Vars.S_lambda**2 + 128*Vars.S_n**2 + 99*Vars.Sigma_ell**2 + 24*Vars.Sigma_lambda**2 + 24*Vars.Sigma_n**2)/(336*Cons.M**4) + 0.0*I
    Vars.hHat_spin_Symm_2_0_5 = Cons.nu*(319*sqrt(5)*Vars.S_ell/(882*Cons.M**2) + 1775*sqrt(5)*Vars.Sigma_ell*Cons.delta/(9408*Cons.M**2)) + 179869*sqrt(5)*Vars.S_ell/(74088*Cons.M**2) + 289127*sqrt(5)*Vars.Sigma_ell*Cons.delta/(790272*Cons.M**2) + 0.0*I
    Vars.hHat_spin_Symm_4_0_3 = 5*Vars.S_ell/(21*Cons.M**2) + 19*Vars.Sigma_ell*Cons.delta/(112*Cons.M**2) + 0.0*I
    Vars.hHat_spin_Symm_4_0_4 = -2*Cons.delta*(3*Vars.S_ell*Vars.Sigma_ell + Vars.S_lambda*Vars.Sigma_lambda + Vars.S_n*Vars.Sigma_n)/(63*Cons.M**4) + Cons.nu*(6*Vars.Sigma_ell**2 + 2*Vars.Sigma_lambda**2 + 2*Vars.Sigma_n**2)/(63*Cons.M**4) + (-192*Vars.S_ell**2 - 64*Vars.S_lambda**2 - 64*Vars.S_n**2 - 53*Vars.Sigma_ell**2 - 12*Vars.Sigma_lambda**2 - 12*Vars.Sigma_n**2)/(2016*Cons.M**4) + 0.0*I
    Vars.hHat_spin_Symm_4_0_5 = Cons.nu*(989419*Vars.S_ell/(232848*Cons.M**2) + 361055*Vars.Sigma_ell*Cons.delta/(155232*Cons.M**2)) - 626077*Vars.S_ell/(724416*Cons.M**2) - 179791*Vars.Sigma_ell*Cons.delta/(271656*Cons.M**2) + 0.0*I
    Vars.hHat_spin_Symm_6_0_5 = Cons.nu*(19135*sqrt(13)*Vars.S_ell/(144144*Cons.M**2) + 16655*sqrt(13)*Vars.Sigma_ell*Cons.delta/(192192*Cons.M**2)) - 373355*sqrt(13)*Vars.S_ell/(12108096*Cons.M**2) - 407135*sqrt(13)*Vars.Sigma_ell*Cons.delta/(16144128*Cons.M**2) + 0.0*I
    Vars.hHat_spin_Asymm_2_1_3 = -0.611297497215587*Cons.delta*(-I*Vars.Sigma_lambda + Vars.Sigma_n)/Cons.M**2 - sqrt(30)*(-61*I*Vars.S_lambda + 61*Vars.S_n)/(420*Cons.M**2)
    Vars.hHat_spin_Asymm_2_1_4 = -0.391230398217976*Vars.Sigma_ell*Cons.nu*(-I*Vars.Sigma_lambda + Vars.Sigma_n)/Cons.M**4 + 0.195615199108988*Cons.delta*(Vars.S_ell*(-I*Vars.Sigma_lambda + Vars.Sigma_n) + Vars.Sigma_ell*(-I*Vars.S_lambda + Vars.S_n))/Cons.M**4 + (2.73861278752583*Vars.S_ell*(-I*Vars.S_lambda + Vars.S_n) + 0.912870929175277*Vars.Sigma_ell*(-I*Vars.Sigma_lambda + Vars.Sigma_n))/(7*Cons.M**4)
    Vars.hHat_spin_Asymm_4_1_3 = -17*sqrt(5)*Cons.delta*(-I*Vars.Sigma_lambda + Vars.Sigma_n)/(420*Cons.M**2) - sqrt(5)*(-13*I*Vars.S_lambda + 13*Vars.S_n)/(280*Cons.M**2)
    Vars.hHat_spin_Asymm_4_1_4 = -sqrt(5)*Vars.Sigma_ell*Cons.nu*(-I*Vars.Sigma_lambda + Vars.Sigma_n)/(84*Cons.M**4) + sqrt(5)*Cons.delta*(Vars.S_ell*(-I*Vars.Sigma_lambda + Vars.Sigma_n) + Vars.Sigma_ell*(-I*Vars.S_lambda + Vars.S_n))/(168*Cons.M**4) + sqrt(5)*(3*Vars.S_ell*(-I*Vars.S_lambda + Vars.S_n) + Vars.Sigma_ell*(-I*Vars.Sigma_lambda + Vars.Sigma_n))/(252*Cons.M**4)

@njit(cache=True)
def Modes_2p5(modes,Cons,Vars):
    # (ell, m) = (0, +/- 0)
    Symm = Vars.rhOverM_coeff*(Cons.hHat_0_0_0 + Vars.v**2*(Cons.hHat_0_0_2 + Vars.v*(Vars.hHat_spin_Symm_0_0_3 + Vars.v*(Cons.hHat_0_0_4 + Vars.hHat_spin_Symm_0_0_4))))
    Asymm = 0
    modes[0] = (Symm + Asymm)[0]
    # (ell, m) = (2, +/- 0)
    Symm = Vars.rhOverM_coeff*(Cons.hHat_2_0_0 + Vars.v**2*(Cons.hHat_2_0_2 + Vars.v*(Vars.hHat_spin_Symm_2_0_3 + Vars.v*(Cons.hHat_2_0_4 + Vars.hHat_spin_Symm_2_0_4 + Vars.v*(Cons.hHat_2_0_5 + Vars.hHat_spin_Symm_2_0_5)))))
    Asymm = 0
    modes[6] = (Symm + Asymm)[0]
    # (ell, m) = (2, +/- 1)
    Symm = 0
    Asymm = Vars.rhOverM_coeff*Vars.v**3*(Vars.hHat_spin_Asymm_2_1_3 + Vars.hHat_spin_Asymm_2_1_4*Vars.v)
    modes[7] = (Symm + Asymm)[0]
    modes[5] = (np.conjugate(Symm - Asymm))[0]
    # (ell, m) = (4, +/- 0)
    Symm = Vars.rhOverM_coeff*(Cons.hHat_4_0_0 + Vars.v**2*(Cons.hHat_4_0_2 + Vars.v*(Vars.hHat_spin_Symm_4_0_3 + Vars.v*(Cons.hHat_4_0_4 + Vars.hHat_spin_Symm_4_0_4 + Vars.v*(Cons.hHat_4_0_5 + Vars.hHat_spin_Symm_4_0_5)))))
    Asymm = 0
    modes[20] = (Symm + Asymm)[0]
    # (ell, m) = (4, +/- 1)
    Symm = 0
    Asymm = Vars.rhOverM_coeff*Vars.v**3*(Vars.hHat_spin_Asymm_4_1_3 + Vars.hHat_spin_Asymm_4_1_4*Vars.v)
    modes[21] = (Symm + Asymm)[0]
    modes[19] = (np.conjugate(Symm - Asymm))[0]
    # (ell, m) = (4, +/- 4)
    Symm = Cons.hHat_4_4_5*Vars.rhOverM_coeff*Vars.v**5
    Asymm = 0
    modes[24] = (Symm + Asymm)[0]
    modes[16] = (np.conjugate(Symm - Asymm))[0]
    # (ell, m) = (6, +/- 0)
    Symm = Vars.rhOverM_coeff*(Cons.hHat_6_0_0 + Vars.v**2*(Cons.hHat_6_0_2 + Vars.v**2*(Cons.hHat_6_0_4 + Vars.v*(Cons.hHat_6_0_5 + Vars.hHat_spin_Symm_6_0_5))))
    Asymm = 0
    modes[42] = (Symm + Asymm)[0]
    # (ell, m) = (8, +/- 0)
    Symm = Vars.rhOverM_coeff*(Cons.hHat_8_0_0 + Vars.v**2*(Cons.hHat_8_0_2 + Cons.hHat_8_0_4*Vars.v**2))
    Asymm = 0
    modes[72] = (Symm + Asymm)[0]

@njit(cache=True)
def Recalculate_3p0(Cons,Vars,y):
    Vars.v = np.array([y[0]])
    Vars.rfrak_chi1 = np.array([y[1],y[2]])
    Vars.rfrak_chi2 = np.array([y[3],y[4]])
    Vars.R = y[5]*Cons.wHat + y[6]*Cons.xHat + y[7]*Cons.yHat + y[8]*Cons.zHat
    Vars.R_S1 = exp(Vars.rfrak_chi1[0]*Cons.xHat + Vars.rfrak_chi1[1]*Cons.yHat)
    Vars.R_S2 = exp(Vars.rfrak_chi2[0]*Cons.xHat + Vars.rfrak_chi2[1]*Cons.yHat)
    Vars.chiVec1 = mul(mul(conjugate(normalized(y[5]*Cons.wHat + y[6]*Cons.xHat + y[7]*Cons.yHat + y[8]*Cons.zHat)),mul(mul(mul(Cons.S_chi1,Vars.R_S1),Cons.zHat),mul(conjugate(Vars.R_S1),conjugate(Cons.S_chi1)))),normalized(y[5]*Cons.wHat + y[6]*Cons.xHat + y[7]*Cons.yHat + y[8]*Cons.zHat))
    Vars.chiVec2 = mul(mul(conjugate(normalized(y[5]*Cons.wHat + y[6]*Cons.xHat + y[7]*Cons.yHat + y[8]*Cons.zHat)),mul(mul(mul(Cons.S_chi2,Vars.R_S2),Cons.zHat),mul(conjugate(Vars.R_S2),conjugate(Cons.S_chi2)))),normalized(y[5]*Cons.wHat + y[6]*Cons.xHat + y[7]*Cons.yHat + y[8]*Cons.zHat))
    Vars.chi1_n = np.array([Vars.chiVec1[1]])
    Vars.chi1_lambda = np.array([Vars.chiVec1[2]])
    Vars.chi1_ell = np.array([Vars.chiVec1[3]])
    Vars.chi2_n = np.array([Vars.chiVec2[1]])
    Vars.chi2_lambda = np.array([Vars.chiVec2[2]])
    Vars.chi2_ell = np.array([Vars.chiVec2[3]])
    Vars.S_ell = Cons.M1**2*Vars.chi1_ell + Cons.M2**2*Vars.chi2_ell
    Vars.S_n = Cons.M1**2*Vars.chi1_n + Cons.M2**2*Vars.chi2_n
    Vars.S_lambda = Cons.M1**2*Vars.chi1_lambda + Cons.M2**2*Vars.chi2_lambda
    Vars.Sigma_ell = Cons.M*(-Cons.M1*Vars.chi1_ell + Cons.M2*Vars.chi2_ell)
    Vars.Sigma_n = Cons.M*(-Cons.M1*Vars.chi1_n + Cons.M2*Vars.chi2_n)
    Vars.Sigma_lambda = Cons.M*(-Cons.M1*Vars.chi1_lambda + Cons.M2*Vars.chi2_lambda)
    Vars.rhOverM_coeff = Cons.nu*Vars.v**2*sqrt(pi)
    Vars.hHat_spin_Symm_0_0_3 = 14*Vars.S_ell/(3*Cons.M**2) + 2*Vars.Sigma_ell*Cons.delta/Cons.M**2 + 0.0*I
    Vars.hHat_spin_Symm_0_0_4 = -4*Cons.delta*(3*Vars.S_ell*Vars.Sigma_ell + Vars.S_lambda*Vars.Sigma_n + Vars.S_n*Vars.Sigma_n)/(3*Cons.M**4) + Cons.nu*(12*Vars.Sigma_ell**2 + 4*Vars.Sigma_lambda**2 + 4*Vars.Sigma_n**2)/(3*Cons.M**4) + (-48*Vars.S_ell**2 - 16*Vars.S_lambda**2 - 16*Vars.S_n**2 - 12*Vars.Sigma_ell**2 - 3*Vars.Sigma_lambda**2 - 3*Vars.Sigma_n**2)/(12*Cons.M**4) + 0.0*I
    Vars.hHat_spin_Symm_2_0_3 = 32*sqrt(5)*Vars.S_ell/(21*Cons.M**2) + 419*sqrt(5)*Vars.Sigma_ell*Cons.delta/(560*Cons.M**2) + 0.0*I
    Vars.hHat_spin_Symm_2_0_4 = -8*sqrt(5)*Cons.delta*(3*Vars.S_ell*Vars.Sigma_ell + Vars.S_lambda*Vars.Sigma_lambda + Vars.S_n*Vars.Sigma_n)/(21*Cons.M**4) + 8*sqrt(5)*Cons.nu*(3*Vars.Sigma_ell**2 + Vars.Sigma_lambda**2 + Vars.Sigma_n**2)/(21*Cons.M**4) - sqrt(5)*(384*Vars.S_ell**2 + 128*Vars.S_lambda**2 + 128*Vars.S_n**2 + 99*Vars.Sigma_ell**2 + 24*Vars.Sigma_lambda**2 + 24*Vars.Sigma_n**2)/(336*Cons.M**4) + 0.0*I
    Vars.hHat_spin_Symm_2_0_5 = Cons.nu*(319*sqrt(5)*Vars.S_ell/(882*Cons.M**2) + 1775*sqrt(5)*Vars.Sigma_ell*Cons.delta/(9408*Cons.M**2)) + 179869*sqrt(5)*Vars.S_ell/(74088*Cons.M**2) + 289127*sqrt(5)*Vars.Sigma_ell*Cons.delta/(790272*Cons.M**2) + 0.0*I
    Vars.hHat_spin_Symm_2_0_6 = Cons.delta*(23*sqrt(5)*Vars.Sigma_ell*pi/(1344*Cons.M**2) + 173191*sqrt(5)*Vars.S_ell*Vars.Sigma_ell/(56448*Cons.M**4)) + Cons.nu*(-31*sqrt(5)*Vars.S_ell*Vars.Sigma_ell*Cons.delta/(42*Cons.M**4) - sqrt(5)*(20832*Vars.S_ell**2 + 142529*Vars.Sigma_ell**2)/(28224*Cons.M**4)) + 31*sqrt(5)*Vars.Sigma_ell**2*Cons.nu**2/(42*Cons.M**4) + sqrt(5)*(63808*Vars.S_ell**2 + 99197*Vars.Sigma_ell**2)/(75264*Cons.M**4) + 0.0*I
    Vars.hHat_spin_Symm_4_0_3 = 5*Vars.S_ell/(21*Cons.M**2) + 19*Vars.Sigma_ell*Cons.delta/(112*Cons.M**2) + 0.0*I
    Vars.hHat_spin_Symm_4_0_4 = -2*Cons.delta*(3*Vars.S_ell*Vars.Sigma_ell + Vars.S_lambda*Vars.Sigma_lambda + Vars.S_n*Vars.Sigma_n)/(63*Cons.M**4) + Cons.nu*(6*Vars.Sigma_ell**2 + 2*Vars.Sigma_lambda**2 + 2*Vars.Sigma_n**2)/(63*Cons.M**4) + (-192*Vars.S_ell**2 - 64*Vars.S_lambda**2 - 64*Vars.S_n**2 - 53*Vars.Sigma_ell**2 - 12*Vars.Sigma_lambda**2 - 12*Vars.Sigma_n**2)/(2016*Cons.M**4) + 0.0*I
    Vars.hHat_spin_Symm_4_0_5 = Cons.nu*(989419*Vars.S_ell/(232848*Cons.M**2) + 361055*Vars.Sigma_ell*Cons.delta/(155232*Cons.M**2)) - 626077*Vars.S_ell/(724416*Cons.M**2) - 179791*Vars.Sigma_ell*Cons.delta/(271656*Cons.M**2) + 0.0*I
    Vars.hHat_spin_Symm_4_0_6 = Cons.delta*(25*Vars.Sigma_ell*pi/(4032*Cons.M**2) + 1303525*Vars.S_ell*Vars.Sigma_ell/(620928*Cons.M**4)) + Cons.nu*(-26435*Vars.S_ell*Vars.Sigma_ell*Cons.delta/(11088*Cons.M**4) - (23685760*Vars.S_ell**2 + 32586265*Vars.Sigma_ell**2)/(9934848*Cons.M**4)) + 26435*Vars.Sigma_ell**2*Cons.nu**2/(11088*Cons.M**4) + (1380815*Vars.S_ell**2/931392 + 352565*Vars.Sigma_ell**2/516096)/Cons.M**4 + 0.0*I
    Vars.hHat_spin_Symm_6_0_5 = Cons.nu*(19135*sqrt(13)*Vars.S_ell/(144144*Cons.M**2) + 16655*sqrt(13)*Vars.Sigma_ell*Cons.delta/(192192*Cons.M**2)) - 373355*sqrt(13)*Vars.S_ell/(12108096*Cons.M**2) - 407135*sqrt(13)*Vars.Sigma_ell*Cons.delta/(16144128*Cons.M**2) + 0.0*I
    Vars.hHat_spin_Symm_6_0_6 = Cons.nu*(-215*sqrt(13)*Vars.S_ell*Vars.Sigma_ell*Cons.delta/(6864*Cons.M**4) - sqrt(13)*(577920*Vars.S_ell**2 + 495875*Vars.Sigma_ell**2)/(18450432*Cons.M**4)) + 8005*sqrt(13)*Vars.S_ell*Vars.Sigma_ell*Cons.delta/(576576*Cons.M**4) + 215*sqrt(13)*Vars.Sigma_ell**2*Cons.nu**2/(6864*Cons.M**4) + sqrt(13)*(782720*Vars.S_ell**2 + 381315*Vars.Sigma_ell**2)/(73801728*Cons.M**4) + 0.0*I
    Vars.hHat_spin_Asymm_2_1_3 = -0.611297497215587*Cons.delta*(-I*Vars.Sigma_lambda + Vars.Sigma_n)/Cons.M**2 - sqrt(30)*(-61*I*Vars.S_lambda + 61*Vars.S_n)/(420*Cons.M**2)
    Vars.hHat_spin_Asymm_2_1_4 = -0.391230398217976*Vars.Sigma_ell*Cons.nu*(-I*Vars.Sigma_lambda + Vars.Sigma_n)/Cons.M**4 + 0.195615199108988*Cons.delta*(Vars.S_ell*(-I*Vars.Sigma_lambda + Vars.Sigma_n) + Vars.Sigma_ell*(-I*Vars.S_lambda + Vars.S_n))/Cons.M**4 + (2.73861278752583*Vars.S_ell*(-I*Vars.S_lambda + Vars.S_n) + 0.912870929175277*Vars.Sigma_ell*(-I*Vars.Sigma_lambda + Vars.Sigma_n))/(7*Cons.M**4)
    Vars.hHat_spin_Asymm_4_1_3 = -17*sqrt(5)*Cons.delta*(-I*Vars.Sigma_lambda + Vars.Sigma_n)/(420*Cons.M**2) - sqrt(5)*(-13*I*Vars.S_lambda + 13*Vars.S_n)/(280*Cons.M**2)
    Vars.hHat_spin_Asymm_4_1_4 = -sqrt(5)*Vars.Sigma_ell*Cons.nu*(-I*Vars.Sigma_lambda + Vars.Sigma_n)/(84*Cons.M**4) + sqrt(5)*Cons.delta*(Vars.S_ell*(-I*Vars.Sigma_lambda + Vars.Sigma_n) + Vars.Sigma_ell*(-I*Vars.S_lambda + Vars.S_n))/(168*Cons.M**4) + sqrt(5)*(3*Vars.S_ell*(-I*Vars.S_lambda + Vars.S_n) + Vars.Sigma_ell*(-I*Vars.Sigma_lambda + Vars.Sigma_n))/(252*Cons.M**4)

@njit(cache=True)
def Modes_3p0(modes,Cons,Vars):
    # (ell, m) = (0, +/- 0)
    Symm = Vars.rhOverM_coeff*(Cons.hHat_0_0_0 + Vars.v**2*(Cons.hHat_0_0_2 + Vars.v*(Vars.hHat_spin_Symm_0_0_3 + Vars.v*(Cons.hHat_0_0_4 + Cons.hHat_0_0_6*Vars.v**2 + Vars.hHat_spin_Symm_0_0_4))))
    Asymm = 0
    modes[0] = (Symm + Asymm)[0]
    # (ell, m) = (1, +/- 1)
    Symm = Cons.hHat_1_1_6*Vars.rhOverM_coeff*Vars.v**6
    Asymm = 0
    modes[3] = (Symm + Asymm)[0]
    modes[1] = (-np.conjugate(Symm - Asymm))[0]
    # (ell, m) = (2, +/- 0)
    Symm = Vars.rhOverM_coeff*(Cons.hHat_2_0_0 + Vars.v**2*(Cons.hHat_2_0_2 + Vars.v*(Vars.hHat_spin_Symm_2_0_3 + Vars.v*(Cons.hHat_2_0_4 + Vars.hHat_spin_Symm_2_0_4 + Vars.v*(Cons.hHat_2_0_5 + Vars.hHat_spin_Symm_2_0_5 + Vars.v*(Cons.hHat_2_0_6 + Vars.hHat_spin_Symm_2_0_6))))))
    Asymm = 0
    modes[6] = (Symm + Asymm)[0]
    # (ell, m) = (2, +/- 1)
    Symm = 0
    Asymm = Vars.rhOverM_coeff*Vars.v**3*(Vars.hHat_spin_Asymm_2_1_3 + Vars.hHat_spin_Asymm_2_1_4*Vars.v)
    modes[7] = (Symm + Asymm)[0]
    modes[5] = (np.conjugate(Symm - Asymm))[0]
    # (ell, m) = (3, +/- 1)
    Symm = Cons.hHat_3_1_6*Vars.rhOverM_coeff*Vars.v**6
    Asymm = 0
    modes[13] = (Symm + Asymm)[0]
    modes[11] = (-np.conjugate(Symm - Asymm))[0]
    # (ell, m) = (3, +/- 3)
    Symm = Cons.hHat_3_3_6*Vars.rhOverM_coeff*Vars.v**6
    Asymm = 0
    modes[15] = (Symm + Asymm)[0]
    modes[9] = (-np.conjugate(Symm - Asymm))[0]
    # (ell, m) = (4, +/- 0)
    Symm = Vars.rhOverM_coeff*(Cons.hHat_4_0_0 + Vars.v**2*(Cons.hHat_4_0_2 + Vars.v*(Vars.hHat_spin_Symm_4_0_3 + Vars.v*(Cons.hHat_4_0_4 + Vars.hHat_spin_Symm_4_0_4 + Vars.v*(Cons.hHat_4_0_5 + Vars.hHat_spin_Symm_4_0_5 + Vars.v*(Cons.hHat_4_0_6 + Vars.hHat_spin_Symm_4_0_6))))))
    Asymm = 0
    modes[20] = (Symm + Asymm)[0]
    # (ell, m) = (4, +/- 1)
    Symm = 0
    Asymm = Vars.rhOverM_coeff*Vars.v**3*(Vars.hHat_spin_Asymm_4_1_3 + Vars.hHat_spin_Asymm_4_1_4*Vars.v)
    modes[21] = (Symm + Asymm)[0]
    modes[19] = (np.conjugate(Symm - Asymm))[0]
    # (ell, m) = (4, +/- 4)
    Symm = Cons.hHat_4_4_5*Vars.rhOverM_coeff*Vars.v**5
    Asymm = 0
    modes[24] = (Symm + Asymm)[0]
    modes[16] = (np.conjugate(Symm - Asymm))[0]
    # (ell, m) = (5, +/- 1)
    Symm = Cons.hHat_5_1_6*Vars.rhOverM_coeff*Vars.v**6
    Asymm = 0
    modes[31] = (Symm + Asymm)[0]
    modes[29] = (-np.conjugate(Symm - Asymm))[0]
    # (ell, m) = (5, +/- 3)
    Symm = Cons.hHat_5_3_6*Vars.rhOverM_coeff*Vars.v**6
    Asymm = 0
    modes[33] = (Symm + Asymm)[0]
    modes[27] = (-np.conjugate(Symm - Asymm))[0]
    # (ell, m) = (5, +/- 5)
    Symm = Cons.hHat_5_5_6*Vars.rhOverM_coeff*Vars.v**6
    Asymm = 0
    modes[35] = (Symm + Asymm)[0]
    modes[25] = (-np.conjugate(Symm - Asymm))[0]
    # (ell, m) = (6, +/- 0)
    Symm = Vars.rhOverM_coeff*(Cons.hHat_6_0_0 + Vars.v**2*(Cons.hHat_6_0_2 + Vars.v**2*(Cons.hHat_6_0_4 + Vars.v*(Cons.hHat_6_0_5 + Vars.hHat_spin_Symm_6_0_5 + Vars.v*(Cons.hHat_6_0_6 + Vars.hHat_spin_Symm_6_0_6)))))
    Asymm = 0
    modes[42] = (Symm + Asymm)[0]
    # (ell, m) = (8, +/- 0)
    Symm = Vars.rhOverM_coeff*(Cons.hHat_8_0_0 + Vars.v**2*(Cons.hHat_8_0_2 + Vars.v**2*(Cons.hHat_8_0_4 + Cons.hHat_8_0_6*Vars.v**2)))
    Asymm = 0
    modes[72] = (Symm + Asymm)[0]

@njit(cache=True)
def Recalculate_3p5(Cons,Vars,y):
    Vars.v = np.array([y[0]])
    Vars.rfrak_chi1 = np.array([y[1],y[2]])
    Vars.rfrak_chi2 = np.array([y[3],y[4]])
    Vars.R = y[5]*Cons.wHat + y[6]*Cons.xHat + y[7]*Cons.yHat + y[8]*Cons.zHat
    Vars.R_S1 = exp(Vars.rfrak_chi1[0]*Cons.xHat + Vars.rfrak_chi1[1]*Cons.yHat)
    Vars.R_S2 = exp(Vars.rfrak_chi2[0]*Cons.xHat + Vars.rfrak_chi2[1]*Cons.yHat)
    Vars.chiVec1 = mul(mul(conjugate(normalized(y[5]*Cons.wHat + y[6]*Cons.xHat + y[7]*Cons.yHat + y[8]*Cons.zHat)),mul(mul(mul(Cons.S_chi1,Vars.R_S1),Cons.zHat),mul(conjugate(Vars.R_S1),conjugate(Cons.S_chi1)))),normalized(y[5]*Cons.wHat + y[6]*Cons.xHat + y[7]*Cons.yHat + y[8]*Cons.zHat))
    Vars.chiVec2 = mul(mul(conjugate(normalized(y[5]*Cons.wHat + y[6]*Cons.xHat + y[7]*Cons.yHat + y[8]*Cons.zHat)),mul(mul(mul(Cons.S_chi2,Vars.R_S2),Cons.zHat),mul(conjugate(Vars.R_S2),conjugate(Cons.S_chi2)))),normalized(y[5]*Cons.wHat + y[6]*Cons.xHat + y[7]*Cons.yHat + y[8]*Cons.zHat))
    Vars.chi1_n = np.array([Vars.chiVec1[1]])
    Vars.chi1_lambda = np.array([Vars.chiVec1[2]])
    Vars.chi1_ell = np.array([Vars.chiVec1[3]])
    Vars.chi2_n = np.array([Vars.chiVec2[1]])
    Vars.chi2_lambda = np.array([Vars.chiVec2[2]])
    Vars.chi2_ell = np.array([Vars.chiVec2[3]])
    Vars.S_ell = Cons.M1**2*Vars.chi1_ell + Cons.M2**2*Vars.chi2_ell
    Vars.S_n = Cons.M1**2*Vars.chi1_n + Cons.M2**2*Vars.chi2_n
    Vars.S_lambda = Cons.M1**2*Vars.chi1_lambda + Cons.M2**2*Vars.chi2_lambda
    Vars.Sigma_ell = Cons.M*(-Cons.M1*Vars.chi1_ell + Cons.M2*Vars.chi2_ell)
    Vars.Sigma_n = Cons.M*(-Cons.M1*Vars.chi1_n + Cons.M2*Vars.chi2_n)
    Vars.Sigma_lambda = Cons.M*(-Cons.M1*Vars.chi1_lambda + Cons.M2*Vars.chi2_lambda)
    Vars.x = Vars.v**2
    Vars.rhOverM_coeff = Cons.nu*Vars.v**2*sqrt(pi)
    Vars.hHat_spin_Symm_0_0_3 = 14*Vars.S_ell/(3*Cons.M**2) + 2*Vars.Sigma_ell*Cons.delta/Cons.M**2 + 0.0*I
    Vars.hHat_spin_Symm_0_0_4 = -4*Cons.delta*(3*Vars.S_ell*Vars.Sigma_ell + Vars.S_lambda*Vars.Sigma_n + Vars.S_n*Vars.Sigma_n)/(3*Cons.M**4) + Cons.nu*(12*Vars.Sigma_ell**2 + 4*Vars.Sigma_lambda**2 + 4*Vars.Sigma_n**2)/(3*Cons.M**4) + (-48*Vars.S_ell**2 - 16*Vars.S_lambda**2 - 16*Vars.S_n**2 - 12*Vars.Sigma_ell**2 - 3*Vars.Sigma_lambda**2 - 3*Vars.Sigma_n**2)/(12*Cons.M**4) + 0.0*I
    Vars.hHat_spin_Symm_2_0_3 = 32*sqrt(5)*Vars.S_ell/(21*Cons.M**2) + 419*sqrt(5)*Vars.Sigma_ell*Cons.delta/(560*Cons.M**2) + 0.0*I
    Vars.hHat_spin_Symm_2_0_4 = -8*sqrt(5)*Cons.delta*(3*Vars.S_ell*Vars.Sigma_ell + Vars.S_lambda*Vars.Sigma_lambda + Vars.S_n*Vars.Sigma_n)/(21*Cons.M**4) + 8*sqrt(5)*Cons.nu*(3*Vars.Sigma_ell**2 + Vars.Sigma_lambda**2 + Vars.Sigma_n**2)/(21*Cons.M**4) - sqrt(5)*(384*Vars.S_ell**2 + 128*Vars.S_lambda**2 + 128*Vars.S_n**2 + 99*Vars.Sigma_ell**2 + 24*Vars.Sigma_lambda**2 + 24*Vars.Sigma_n**2)/(336*Cons.M**4) + 0.0*I
    Vars.hHat_spin_Symm_2_0_5 = Cons.nu*(319*sqrt(5)*Vars.S_ell/(882*Cons.M**2) + 1775*sqrt(5)*Vars.Sigma_ell*Cons.delta/(9408*Cons.M**2)) + 179869*sqrt(5)*Vars.S_ell/(74088*Cons.M**2) + 289127*sqrt(5)*Vars.Sigma_ell*Cons.delta/(790272*Cons.M**2) + 0.0*I
    Vars.hHat_spin_Symm_2_0_6 = Cons.delta*(23*sqrt(5)*Vars.Sigma_ell*pi/(1344*Cons.M**2) + 173191*sqrt(5)*Vars.S_ell*Vars.Sigma_ell/(56448*Cons.M**4)) + Cons.nu*(-31*sqrt(5)*Vars.S_ell*Vars.Sigma_ell*Cons.delta/(42*Cons.M**4) - sqrt(5)*(20832*Vars.S_ell**2 + 142529*Vars.Sigma_ell**2)/(28224*Cons.M**4)) + 31*sqrt(5)*Vars.Sigma_ell**2*Cons.nu**2/(42*Cons.M**4) + sqrt(5)*(63808*Vars.S_ell**2 + 99197*Vars.Sigma_ell**2)/(75264*Cons.M**4) + 0.0*I
    Vars.hHat_spin_Symm_2_0_7 = Cons.delta*(14713254113*sqrt(5)*Vars.Sigma_ell/(11266117632*Cons.M**2) - 256*sqrt(5)*Vars.S_ell*Vars.Sigma_ell*pi/(21*Cons.M**4) + sqrt(5)*(1605952*Vars.S_ell**2*Vars.Sigma_ell + 106463*Vars.Sigma_ell**3)/(48384*Cons.M**6)) + Cons.nu**2*(-34739*sqrt(5)*Vars.S_ell/(49896*Cons.M**2) - 66403*sqrt(5)*Vars.Sigma_ell*Cons.delta/(177408*Cons.M**2)) + Cons.nu*(Cons.delta*(-66021037*sqrt(5)*Vars.Sigma_ell/(5588352*Cons.M**2) - 733*sqrt(5)*Vars.Sigma_ell**3/(84*Cons.M**6)) - 53090267*sqrt(5)*Vars.S_ell/(2095632*Cons.M**2) + 256*sqrt(5)*Vars.Sigma_ell**2*pi/(21*Cons.M**4) - 1603*sqrt(5)*Vars.S_ell*Vars.Sigma_ell**2/(27*Cons.M**6)) + 9724899179*sqrt(5)*Vars.S_ell/(1056198528*Cons.M**2) - sqrt(5)*pi*(3072*Vars.S_ell**2 + 779*Vars.Sigma_ell**2)/(252*Cons.M**4) + sqrt(5)*(36992*Vars.S_ell**3 + 22539*Vars.S_ell*Vars.Sigma_ell**2)/(1512*Cons.M**6) + 0.0*I
    Vars.hHat_spin_Symm_3_3_7 = -8*sqrt(35)*Vars.Sigma_ell*Cons.nu*Vars.x**3.5/(105*Cons.M**2) + 0.0*I
    Vars.hHat_spin_Symm_3_1_7 = -1.04744587313276*Vars.Sigma_ell*Cons.nu*Vars.x**3.5/Cons.M**2 + 0.0*I
    Vars.hHat_spin_Symm_4_0_3 = 5*Vars.S_ell/(21*Cons.M**2) + 19*Vars.Sigma_ell*Cons.delta/(112*Cons.M**2) + 0.0*I
    Vars.hHat_spin_Symm_4_0_4 = -2*Cons.delta*(3*Vars.S_ell*Vars.Sigma_ell + Vars.S_lambda*Vars.Sigma_lambda + Vars.S_n*Vars.Sigma_n)/(63*Cons.M**4) + Cons.nu*(6*Vars.Sigma_ell**2 + 2*Vars.Sigma_lambda**2 + 2*Vars.Sigma_n**2)/(63*Cons.M**4) + (-192*Vars.S_ell**2 - 64*Vars.S_lambda**2 - 64*Vars.S_n**2 - 53*Vars.Sigma_ell**2 - 12*Vars.Sigma_lambda**2 - 12*Vars.Sigma_n**2)/(2016*Cons.M**4) + 0.0*I
    Vars.hHat_spin_Symm_4_0_5 = Cons.nu*(989419*Vars.S_ell/(232848*Cons.M**2) + 361055*Vars.Sigma_ell*Cons.delta/(155232*Cons.M**2)) - 626077*Vars.S_ell/(724416*Cons.M**2) - 179791*Vars.Sigma_ell*Cons.delta/(271656*Cons.M**2) + 0.0*I
    Vars.hHat_spin_Symm_4_0_6 = Cons.delta*(25*Vars.Sigma_ell*pi/(4032*Cons.M**2) + 1303525*Vars.S_ell*Vars.Sigma_ell/(620928*Cons.M**4)) + Cons.nu*(-26435*Vars.S_ell*Vars.Sigma_ell*Cons.delta/(11088*Cons.M**4) - (23685760*Vars.S_ell**2 + 32586265*Vars.Sigma_ell**2)/(9934848*Cons.M**4)) + 26435*Vars.Sigma_ell**2*Cons.nu**2/(11088*Cons.M**4) + (1380815*Vars.S_ell**2/931392 + 352565*Vars.Sigma_ell**2/516096)/Cons.M**4 + 0.0*I
    Vars.hHat_spin_Symm_4_0_7 = Cons.delta*(25969272883*Vars.Sigma_ell/(39943507968*Cons.M**2) - 64*Vars.S_ell*Vars.Sigma_ell*pi/(63*Cons.M**4) + (2971*Vars.S_ell**2*Vars.Sigma_ell/2268 - 55*Vars.Sigma_ell**3/6912)/Cons.M**6) + Cons.nu**2*(4773259*Vars.S_ell/(5189184*Cons.M**2) + 5099915*Vars.Sigma_ell*Cons.delta/(6918912*Cons.M**2)) + Cons.nu*(Cons.delta*(-160373197*Vars.Sigma_ell/(32288256*Cons.M**2) - Vars.Sigma_ell**3/(84*Cons.M**6)) - 712823899*Vars.S_ell/(163459296*Cons.M**2) + 64*Vars.Sigma_ell**2*pi/(63*Cons.M**4) - 109*Vars.S_ell*Vars.Sigma_ell**2/(81*Cons.M**6)) + 348470641109*Vars.S_ell/(329533940736*Cons.M**2) + pi*(-768*Vars.S_ell**2 - 193*Vars.Sigma_ell**2)/(756*Cons.M**4) + (11776*Vars.S_ell**3 + 2887*Vars.S_ell*Vars.Sigma_ell**2)/(9072*Cons.M**6) + 0.0*I
    Vars.hHat_spin_Symm_6_0_5 = Cons.nu*(19135*sqrt(13)*Vars.S_ell/(144144*Cons.M**2) + 16655*sqrt(13)*Vars.Sigma_ell*Cons.delta/(192192*Cons.M**2)) - 373355*sqrt(13)*Vars.S_ell/(12108096*Cons.M**2) - 407135*sqrt(13)*Vars.Sigma_ell*Cons.delta/(16144128*Cons.M**2) + 0.0*I
    Vars.hHat_spin_Symm_6_0_6 = Cons.nu*(-215*sqrt(13)*Vars.S_ell*Vars.Sigma_ell*Cons.delta/(6864*Cons.M**4) - sqrt(13)*(577920*Vars.S_ell**2 + 495875*Vars.Sigma_ell**2)/(18450432*Cons.M**4)) + 8005*sqrt(13)*Vars.S_ell*Vars.Sigma_ell*Cons.delta/(576576*Cons.M**4) + 215*sqrt(13)*Vars.Sigma_ell**2*Cons.nu**2/(6864*Cons.M**4) + sqrt(13)*(782720*Vars.S_ell**2 + 381315*Vars.Sigma_ell**2)/(73801728*Cons.M**4) + 0.0*I
    Vars.hHat_spin_Symm_6_0_7 = Cons.nu**2*(198335*sqrt(13)*Vars.S_ell/(247104*Cons.M**2) + 10165*sqrt(13)*Vars.Sigma_ell*Cons.delta/(20592*Cons.M**2)) + Cons.nu*(-10654115*sqrt(13)*Vars.S_ell/(15567552*Cons.M**2) - 7026595*sqrt(13)*Vars.Sigma_ell*Cons.delta/(13837824*Cons.M**2)) + 1986903733*sqrt(13)*Vars.S_ell/(15692092416*Cons.M**2) + 106983697*sqrt(13)*Vars.Sigma_ell*Cons.delta/(951035904*Cons.M**2) + 0.0*I
    Vars.hHat_spin_Symm_8_0_7 = Cons.nu**2*(349685*sqrt(17)*Vars.S_ell/(3150576*Cons.M**2) + 99085*sqrt(17)*Vars.Sigma_ell*Cons.delta/(1400256*Cons.M**2)) + Cons.nu*(-2586845*sqrt(17)*Vars.S_ell/(37806912*Cons.M**2) - 2640475*sqrt(17)*Vars.Sigma_ell*Cons.delta/(50409216*Cons.M**2)) + 7786903*sqrt(17)*Vars.S_ell/(680524416*Cons.M**2) + 9147661*sqrt(17)*Vars.Sigma_ell*Cons.delta/(907365888*Cons.M**2) + 0.0*I
    Vars.hHat_spin_Asymm_2_1_3 = -0.611297497215587*Cons.delta*(-I*Vars.Sigma_lambda + Vars.Sigma_n)/Cons.M**2 - sqrt(30)*(-61*I*Vars.S_lambda + 61*Vars.S_n)/(420*Cons.M**2)
    Vars.hHat_spin_Asymm_2_1_4 = -0.391230398217976*Vars.Sigma_ell*Cons.nu*(-I*Vars.Sigma_lambda + Vars.Sigma_n)/Cons.M**4 + 0.195615199108988*Cons.delta*(Vars.S_ell*(-I*Vars.Sigma_lambda + Vars.Sigma_n) + Vars.Sigma_ell*(-I*Vars.S_lambda + Vars.S_n))/Cons.M**4 + (2.73861278752583*Vars.S_ell*(-I*Vars.S_lambda + Vars.S_n) + 0.912870929175277*Vars.Sigma_ell*(-I*Vars.Sigma_lambda + Vars.Sigma_n))/(7*Cons.M**4)
    Vars.hHat_spin_Asymm_4_1_3 = -17*sqrt(5)*Cons.delta*(-I*Vars.Sigma_lambda + Vars.Sigma_n)/(420*Cons.M**2) - sqrt(5)*(-13*I*Vars.S_lambda + 13*Vars.S_n)/(280*Cons.M**2)
    Vars.hHat_spin_Asymm_4_1_4 = -sqrt(5)*Vars.Sigma_ell*Cons.nu*(-I*Vars.Sigma_lambda + Vars.Sigma_n)/(84*Cons.M**4) + sqrt(5)*Cons.delta*(Vars.S_ell*(-I*Vars.Sigma_lambda + Vars.Sigma_n) + Vars.Sigma_ell*(-I*Vars.S_lambda + Vars.S_n))/(168*Cons.M**4) + sqrt(5)*(3*Vars.S_ell*(-I*Vars.S_lambda + Vars.S_n) + Vars.Sigma_ell*(-I*Vars.Sigma_lambda + Vars.Sigma_n))/(252*Cons.M**4)

@njit(cache=True)
def Modes_3p5(modes,Cons,Vars):
    # (ell, m) = (0, +/- 0)
    Symm = Vars.rhOverM_coeff*(Cons.hHat_0_0_0 + Vars.v**2*(Cons.hHat_0_0_2 + Vars.v*(Vars.hHat_spin_Symm_0_0_3 + Vars.v*(Cons.hHat_0_0_4 + Cons.hHat_0_0_6*Vars.v**2 + Vars.hHat_spin_Symm_0_0_4))))
    Asymm = 0
    modes[0] = (Symm + Asymm)[0]
    # (ell, m) = (1, +/- 1)
    Symm = Cons.hHat_1_1_6*Vars.rhOverM_coeff*Vars.v**6
    Asymm = 0
    modes[3] = (Symm + Asymm)[0]
    modes[1] = (-np.conjugate(Symm - Asymm))[0]
    # (ell, m) = (2, +/- 0)
    Symm = Vars.rhOverM_coeff*(Cons.hHat_2_0_0 + Vars.v**2*(Cons.hHat_2_0_2 + Vars.v*(Vars.hHat_spin_Symm_2_0_3 + Vars.v*(Cons.hHat_2_0_4 + Vars.hHat_spin_Symm_2_0_4 + Vars.v*(Cons.hHat_2_0_5 + Vars.hHat_spin_Symm_2_0_5 + Vars.v*(Cons.hHat_2_0_6 + Vars.hHat_spin_Symm_2_0_6 + Vars.hHat_spin_Symm_2_0_7*Vars.v))))))
    Asymm = 0
    modes[6] = (Symm + Asymm)[0]
    # (ell, m) = (2, +/- 1)
    Symm = 0
    Asymm = Vars.rhOverM_coeff*Vars.v**3*(Vars.hHat_spin_Asymm_2_1_3 + Vars.hHat_spin_Asymm_2_1_4*Vars.v)
    modes[7] = (Symm + Asymm)[0]
    modes[5] = (np.conjugate(Symm - Asymm))[0]
    # (ell, m) = (3, +/- 1)
    Symm = Vars.rhOverM_coeff*Vars.v**6*(Cons.hHat_3_1_6 + Vars.hHat_spin_Symm_3_1_7*Vars.v)
    Asymm = 0
    modes[13] = (Symm + Asymm)[0]
    modes[11] = (-np.conjugate(Symm - Asymm))[0]
    # (ell, m) = (3, +/- 3)
    Symm = Vars.rhOverM_coeff*Vars.v**6*(Cons.hHat_3_3_6 + Vars.hHat_spin_Symm_3_3_7*Vars.v)
    Asymm = 0
    modes[15] = (Symm + Asymm)[0]
    modes[9] = (-np.conjugate(Symm - Asymm))[0]
    # (ell, m) = (4, +/- 0)
    Symm = Vars.rhOverM_coeff*(Cons.hHat_4_0_0 + Vars.v**2*(Cons.hHat_4_0_2 + Vars.v*(Vars.hHat_spin_Symm_4_0_3 + Vars.v*(Cons.hHat_4_0_4 + Vars.hHat_spin_Symm_4_0_4 + Vars.v*(Cons.hHat_4_0_5 + Vars.hHat_spin_Symm_4_0_5 + Vars.v*(Cons.hHat_4_0_6 + Vars.hHat_spin_Symm_4_0_6 + Vars.hHat_spin_Symm_4_0_7*Vars.v))))))
    Asymm = 0
    modes[20] = (Symm + Asymm)[0]
    # (ell, m) = (4, +/- 1)
    Symm = 0
    Asymm = Vars.rhOverM_coeff*Vars.v**3*(Vars.hHat_spin_Asymm_4_1_3 + Vars.hHat_spin_Asymm_4_1_4*Vars.v)
    modes[21] = (Symm + Asymm)[0]
    modes[19] = (np.conjugate(Symm - Asymm))[0]
    # (ell, m) = (4, +/- 4)
    Symm = Cons.hHat_4_4_5*Vars.rhOverM_coeff*Vars.v**5
    Asymm = 0
    modes[24] = (Symm + Asymm)[0]
    modes[16] = (np.conjugate(Symm - Asymm))[0]
    # (ell, m) = (5, +/- 1)
    Symm = Cons.hHat_5_1_6*Vars.rhOverM_coeff*Vars.v**6
    Asymm = 0
    modes[31] = (Symm + Asymm)[0]
    modes[29] = (-np.conjugate(Symm - Asymm))[0]
    # (ell, m) = (5, +/- 3)
    Symm = Cons.hHat_5_3_6*Vars.rhOverM_coeff*Vars.v**6
    Asymm = 0
    modes[33] = (Symm + Asymm)[0]
    modes[27] = (-np.conjugate(Symm - Asymm))[0]
    # (ell, m) = (5, +/- 5)
    Symm = Cons.hHat_5_5_6*Vars.rhOverM_coeff*Vars.v**6
    Asymm = 0
    modes[35] = (Symm + Asymm)[0]
    modes[25] = (-np.conjugate(Symm - Asymm))[0]
    # (ell, m) = (6, +/- 0)
    Symm = Vars.rhOverM_coeff*(Cons.hHat_6_0_0 + Vars.v**2*(Cons.hHat_6_0_2 + Vars.v**2*(Cons.hHat_6_0_4 + Vars.v*(Cons.hHat_6_0_5 + Vars.hHat_spin_Symm_6_0_5 + Vars.v*(Cons.hHat_6_0_6 + Vars.hHat_spin_Symm_6_0_6 + Vars.hHat_spin_Symm_6_0_7*Vars.v)))))
    Asymm = 0
    modes[42] = (Symm + Asymm)[0]
    # (ell, m) = (8, +/- 0)
    Symm = Vars.rhOverM_coeff*(Cons.hHat_8_0_0 + Vars.v**2*(Cons.hHat_8_0_2 + Vars.v**2*(Cons.hHat_8_0_4 + Vars.v**2*(Cons.hHat_8_0_6 + Vars.hHat_spin_Symm_8_0_7*Vars.v))))
    Asymm = 0
    modes[72] = (Symm + Asymm)[0]

@njit(cache=True)
def Recalculate_4p0(Cons,Vars,y):
    Vars.v = np.array([y[0]])
    Vars.rfrak_chi1 = np.array([y[1],y[2]])
    Vars.rfrak_chi2 = np.array([y[3],y[4]])
    Vars.R = y[5]*Cons.wHat + y[6]*Cons.xHat + y[7]*Cons.yHat + y[8]*Cons.zHat
    Vars.R_S1 = exp(Vars.rfrak_chi1[0]*Cons.xHat + Vars.rfrak_chi1[1]*Cons.yHat)
    Vars.R_S2 = exp(Vars.rfrak_chi2[0]*Cons.xHat + Vars.rfrak_chi2[1]*Cons.yHat)
    Vars.chiVec1 = mul(mul(conjugate(normalized(y[5]*Cons.wHat + y[6]*Cons.xHat + y[7]*Cons.yHat + y[8]*Cons.zHat)),mul(mul(mul(Cons.S_chi1,Vars.R_S1),Cons.zHat),mul(conjugate(Vars.R_S1),conjugate(Cons.S_chi1)))),normalized(y[5]*Cons.wHat + y[6]*Cons.xHat + y[7]*Cons.yHat + y[8]*Cons.zHat))
    Vars.chiVec2 = mul(mul(conjugate(normalized(y[5]*Cons.wHat + y[6]*Cons.xHat + y[7]*Cons.yHat + y[8]*Cons.zHat)),mul(mul(mul(Cons.S_chi2,Vars.R_S2),Cons.zHat),mul(conjugate(Vars.R_S2),conjugate(Cons.S_chi2)))),normalized(y[5]*Cons.wHat + y[6]*Cons.xHat + y[7]*Cons.yHat + y[8]*Cons.zHat))
    Vars.chi1_n = np.array([Vars.chiVec1[1]])
    Vars.chi1_lambda = np.array([Vars.chiVec1[2]])
    Vars.chi1_ell = np.array([Vars.chiVec1[3]])
    Vars.chi2_n = np.array([Vars.chiVec2[1]])
    Vars.chi2_lambda = np.array([Vars.chiVec2[2]])
    Vars.chi2_ell = np.array([Vars.chiVec2[3]])
    Vars.S_ell = Cons.M1**2*Vars.chi1_ell + Cons.M2**2*Vars.chi2_ell
    Vars.S_n = Cons.M1**2*Vars.chi1_n + Cons.M2**2*Vars.chi2_n
    Vars.S_lambda = Cons.M1**2*Vars.chi1_lambda + Cons.M2**2*Vars.chi2_lambda
    Vars.Sigma_ell = Cons.M*(-Cons.M1*Vars.chi1_ell + Cons.M2*Vars.chi2_ell)
    Vars.Sigma_n = Cons.M*(-Cons.M1*Vars.chi1_n + Cons.M2*Vars.chi2_n)
    Vars.Sigma_lambda = Cons.M*(-Cons.M1*Vars.chi1_lambda + Cons.M2*Vars.chi2_lambda)
    Vars.x = Vars.v**2
    Vars.rhOverM_coeff = Cons.nu*Vars.v**2*sqrt(pi)
    Vars.hHat_spin_Symm_0_0_3 = 14*Vars.S_ell/(3*Cons.M**2) + 2*Vars.Sigma_ell*Cons.delta/Cons.M**2 + 0.0*I
    Vars.hHat_spin_Symm_0_0_4 = -4*Cons.delta*(3*Vars.S_ell*Vars.Sigma_ell + Vars.S_lambda*Vars.Sigma_n + Vars.S_n*Vars.Sigma_n)/(3*Cons.M**4) + Cons.nu*(12*Vars.Sigma_ell**2 + 4*Vars.Sigma_lambda**2 + 4*Vars.Sigma_n**2)/(3*Cons.M**4) + (-48*Vars.S_ell**2 - 16*Vars.S_lambda**2 - 16*Vars.S_n**2 - 12*Vars.Sigma_ell**2 - 3*Vars.Sigma_lambda**2 - 3*Vars.Sigma_n**2)/(12*Cons.M**4) + 0.0*I
    Vars.hHat_spin_Symm_2_0_3 = 32*sqrt(5)*Vars.S_ell/(21*Cons.M**2) + 419*sqrt(5)*Vars.Sigma_ell*Cons.delta/(560*Cons.M**2) + 0.0*I
    Vars.hHat_spin_Symm_2_0_4 = -8*sqrt(5)*Cons.delta*(3*Vars.S_ell*Vars.Sigma_ell + Vars.S_lambda*Vars.Sigma_lambda + Vars.S_n*Vars.Sigma_n)/(21*Cons.M**4) + 8*sqrt(5)*Cons.nu*(3*Vars.Sigma_ell**2 + Vars.Sigma_lambda**2 + Vars.Sigma_n**2)/(21*Cons.M**4) - sqrt(5)*(384*Vars.S_ell**2 + 128*Vars.S_lambda**2 + 128*Vars.S_n**2 + 99*Vars.Sigma_ell**2 + 24*Vars.Sigma_lambda**2 + 24*Vars.Sigma_n**2)/(336*Cons.M**4) + 0.0*I
    Vars.hHat_spin_Symm_2_0_5 = Cons.nu*(319*sqrt(5)*Vars.S_ell/(882*Cons.M**2) + 1775*sqrt(5)*Vars.Sigma_ell*Cons.delta/(9408*Cons.M**2)) + 179869*sqrt(5)*Vars.S_ell/(74088*Cons.M**2) + 289127*sqrt(5)*Vars.Sigma_ell*Cons.delta/(790272*Cons.M**2) + 0.0*I
    Vars.hHat_spin_Symm_2_0_6 = Cons.delta*(23*sqrt(5)*Vars.Sigma_ell*pi/(1344*Cons.M**2) + 173191*sqrt(5)*Vars.S_ell*Vars.Sigma_ell/(56448*Cons.M**4)) + Cons.nu*(-31*sqrt(5)*Vars.S_ell*Vars.Sigma_ell*Cons.delta/(42*Cons.M**4) - sqrt(5)*(20832*Vars.S_ell**2 + 142529*Vars.Sigma_ell**2)/(28224*Cons.M**4)) + 31*sqrt(5)*Vars.Sigma_ell**2*Cons.nu**2/(42*Cons.M**4) + sqrt(5)*(63808*Vars.S_ell**2 + 99197*Vars.Sigma_ell**2)/(75264*Cons.M**4) + 0.0*I
    Vars.hHat_spin_Symm_2_0_7 = Cons.delta*(14713254113*sqrt(5)*Vars.Sigma_ell/(11266117632*Cons.M**2) - 256*sqrt(5)*Vars.S_ell*Vars.Sigma_ell*pi/(21*Cons.M**4) + sqrt(5)*(1605952*Vars.S_ell**2*Vars.Sigma_ell + 106463*Vars.Sigma_ell**3)/(48384*Cons.M**6)) + Cons.nu**2*(-34739*sqrt(5)*Vars.S_ell/(49896*Cons.M**2) - 66403*sqrt(5)*Vars.Sigma_ell*Cons.delta/(177408*Cons.M**2)) + Cons.nu*(Cons.delta*(-66021037*sqrt(5)*Vars.Sigma_ell/(5588352*Cons.M**2) - 733*sqrt(5)*Vars.Sigma_ell**3/(84*Cons.M**6)) - 53090267*sqrt(5)*Vars.S_ell/(2095632*Cons.M**2) + 256*sqrt(5)*Vars.Sigma_ell**2*pi/(21*Cons.M**4) - 1603*sqrt(5)*Vars.S_ell*Vars.Sigma_ell**2/(27*Cons.M**6)) + 9724899179*sqrt(5)*Vars.S_ell/(1056198528*Cons.M**2) - sqrt(5)*pi*(3072*Vars.S_ell**2 + 779*Vars.Sigma_ell**2)/(252*Cons.M**4) + sqrt(5)*(36992*Vars.S_ell**3 + 22539*Vars.S_ell*Vars.Sigma_ell**2)/(1512*Cons.M**6) + 0.0*I
    Vars.hHat_spin_Symm_3_3_7 = -8*sqrt(35)*Vars.Sigma_ell*Cons.nu*Vars.x**3.5/(105*Cons.M**2) + 0.0*I
    Vars.hHat_spin_Symm_3_1_7 = -1.04744587313276*Vars.Sigma_ell*Cons.nu*Vars.x**3.5/Cons.M**2 + 0.0*I
    Vars.hHat_spin_Symm_4_0_3 = 5*Vars.S_ell/(21*Cons.M**2) + 19*Vars.Sigma_ell*Cons.delta/(112*Cons.M**2) + 0.0*I
    Vars.hHat_spin_Symm_4_0_4 = -2*Cons.delta*(3*Vars.S_ell*Vars.Sigma_ell + Vars.S_lambda*Vars.Sigma_lambda + Vars.S_n*Vars.Sigma_n)/(63*Cons.M**4) + Cons.nu*(6*Vars.Sigma_ell**2 + 2*Vars.Sigma_lambda**2 + 2*Vars.Sigma_n**2)/(63*Cons.M**4) + (-192*Vars.S_ell**2 - 64*Vars.S_lambda**2 - 64*Vars.S_n**2 - 53*Vars.Sigma_ell**2 - 12*Vars.Sigma_lambda**2 - 12*Vars.Sigma_n**2)/(2016*Cons.M**4) + 0.0*I
    Vars.hHat_spin_Symm_4_0_5 = Cons.nu*(989419*Vars.S_ell/(232848*Cons.M**2) + 361055*Vars.Sigma_ell*Cons.delta/(155232*Cons.M**2)) - 626077*Vars.S_ell/(724416*Cons.M**2) - 179791*Vars.Sigma_ell*Cons.delta/(271656*Cons.M**2) + 0.0*I
    Vars.hHat_spin_Symm_4_0_6 = Cons.delta*(25*Vars.Sigma_ell*pi/(4032*Cons.M**2) + 1303525*Vars.S_ell*Vars.Sigma_ell/(620928*Cons.M**4)) + Cons.nu*(-26435*Vars.S_ell*Vars.Sigma_ell*Cons.delta/(11088*Cons.M**4) - (23685760*Vars.S_ell**2 + 32586265*Vars.Sigma_ell**2)/(9934848*Cons.M**4)) + 26435*Vars.Sigma_ell**2*Cons.nu**2/(11088*Cons.M**4) + (1380815*Vars.S_ell**2/931392 + 352565*Vars.Sigma_ell**2/516096)/Cons.M**4 + 0.0*I
    Vars.hHat_spin_Symm_4_0_7 = Cons.delta*(25969272883*Vars.Sigma_ell/(39943507968*Cons.M**2) - 64*Vars.S_ell*Vars.Sigma_ell*pi/(63*Cons.M**4) + (2971*Vars.S_ell**2*Vars.Sigma_ell/2268 - 55*Vars.Sigma_ell**3/6912)/Cons.M**6) + Cons.nu**2*(4773259*Vars.S_ell/(5189184*Cons.M**2) + 5099915*Vars.Sigma_ell*Cons.delta/(6918912*Cons.M**2)) + Cons.nu*(Cons.delta*(-160373197*Vars.Sigma_ell/(32288256*Cons.M**2) - Vars.Sigma_ell**3/(84*Cons.M**6)) - 712823899*Vars.S_ell/(163459296*Cons.M**2) + 64*Vars.Sigma_ell**2*pi/(63*Cons.M**4) - 109*Vars.S_ell*Vars.Sigma_ell**2/(81*Cons.M**6)) + 348470641109*Vars.S_ell/(329533940736*Cons.M**2) + pi*(-768*Vars.S_ell**2 - 193*Vars.Sigma_ell**2)/(756*Cons.M**4) + (11776*Vars.S_ell**3 + 2887*Vars.S_ell*Vars.Sigma_ell**2)/(9072*Cons.M**6) + 0.0*I
    Vars.hHat_spin_Symm_6_0_5 = Cons.nu*(19135*sqrt(13)*Vars.S_ell/(144144*Cons.M**2) + 16655*sqrt(13)*Vars.Sigma_ell*Cons.delta/(192192*Cons.M**2)) - 373355*sqrt(13)*Vars.S_ell/(12108096*Cons.M**2) - 407135*sqrt(13)*Vars.Sigma_ell*Cons.delta/(16144128*Cons.M**2) + 0.0*I
    Vars.hHat_spin_Symm_6_0_6 = Cons.nu*(-215*sqrt(13)*Vars.S_ell*Vars.Sigma_ell*Cons.delta/(6864*Cons.M**4) - sqrt(13)*(577920*Vars.S_ell**2 + 495875*Vars.Sigma_ell**2)/(18450432*Cons.M**4)) + 8005*sqrt(13)*Vars.S_ell*Vars.Sigma_ell*Cons.delta/(576576*Cons.M**4) + 215*sqrt(13)*Vars.Sigma_ell**2*Cons.nu**2/(6864*Cons.M**4) + sqrt(13)*(782720*Vars.S_ell**2 + 381315*Vars.Sigma_ell**2)/(73801728*Cons.M**4) + 0.0*I
    Vars.hHat_spin_Symm_6_0_7 = Cons.nu**2*(198335*sqrt(13)*Vars.S_ell/(247104*Cons.M**2) + 10165*sqrt(13)*Vars.Sigma_ell*Cons.delta/(20592*Cons.M**2)) + Cons.nu*(-10654115*sqrt(13)*Vars.S_ell/(15567552*Cons.M**2) - 7026595*sqrt(13)*Vars.Sigma_ell*Cons.delta/(13837824*Cons.M**2)) + 1986903733*sqrt(13)*Vars.S_ell/(15692092416*Cons.M**2) + 106983697*sqrt(13)*Vars.Sigma_ell*Cons.delta/(951035904*Cons.M**2) + 0.0*I
    Vars.hHat_spin_Symm_8_0_7 = Cons.nu**2*(349685*sqrt(17)*Vars.S_ell/(3150576*Cons.M**2) + 99085*sqrt(17)*Vars.Sigma_ell*Cons.delta/(1400256*Cons.M**2)) + Cons.nu*(-2586845*sqrt(17)*Vars.S_ell/(37806912*Cons.M**2) - 2640475*sqrt(17)*Vars.Sigma_ell*Cons.delta/(50409216*Cons.M**2)) + 7786903*sqrt(17)*Vars.S_ell/(680524416*Cons.M**2) + 9147661*sqrt(17)*Vars.Sigma_ell*Cons.delta/(907365888*Cons.M**2) + 0.0*I
    Vars.hHat_spin_Asymm_2_1_3 = -0.611297497215587*Cons.delta*(-I*Vars.Sigma_lambda + Vars.Sigma_n)/Cons.M**2 - sqrt(30)*(-61*I*Vars.S_lambda + 61*Vars.S_n)/(420*Cons.M**2)
    Vars.hHat_spin_Asymm_2_1_4 = -0.391230398217976*Vars.Sigma_ell*Cons.nu*(-I*Vars.Sigma_lambda + Vars.Sigma_n)/Cons.M**4 + 0.195615199108988*Cons.delta*(Vars.S_ell*(-I*Vars.Sigma_lambda + Vars.Sigma_n) + Vars.Sigma_ell*(-I*Vars.S_lambda + Vars.S_n))/Cons.M**4 + (2.73861278752583*Vars.S_ell*(-I*Vars.S_lambda + Vars.S_n) + 0.912870929175277*Vars.Sigma_ell*(-I*Vars.Sigma_lambda + Vars.Sigma_n))/(7*Cons.M**4)
    Vars.hHat_spin_Asymm_4_1_3 = -17*sqrt(5)*Cons.delta*(-I*Vars.Sigma_lambda + Vars.Sigma_n)/(420*Cons.M**2) - sqrt(5)*(-13*I*Vars.S_lambda + 13*Vars.S_n)/(280*Cons.M**2)
    Vars.hHat_spin_Asymm_4_1_4 = -sqrt(5)*Vars.Sigma_ell*Cons.nu*(-I*Vars.Sigma_lambda + Vars.Sigma_n)/(84*Cons.M**4) + sqrt(5)*Cons.delta*(Vars.S_ell*(-I*Vars.Sigma_lambda + Vars.Sigma_n) + Vars.Sigma_ell*(-I*Vars.S_lambda + Vars.S_n))/(168*Cons.M**4) + sqrt(5)*(3*Vars.S_ell*(-I*Vars.S_lambda + Vars.S_n) + Vars.Sigma_ell*(-I*Vars.Sigma_lambda + Vars.Sigma_n))/(252*Cons.M**4)

@njit(cache=True)
def Modes_4p0(modes,Cons,Vars):
    # (ell, m) = (0, +/- 0)
    Symm = Vars.rhOverM_coeff*(Cons.hHat_0_0_0 + Vars.v**2*(Cons.hHat_0_0_2 + Vars.v*(Vars.hHat_spin_Symm_0_0_3 + Vars.v*(Cons.hHat_0_0_4 + Cons.hHat_0_0_6*Vars.v**2 + Vars.hHat_spin_Symm_0_0_4))))
    Asymm = 0
    modes[0] = (Symm + Asymm)[0]
    # (ell, m) = (1, +/- 1)
    Symm = Cons.hHat_1_1_6*Vars.rhOverM_coeff*Vars.v**6
    Asymm = 0
    modes[3] = (Symm + Asymm)[0]
    modes[1] = (-np.conjugate(Symm - Asymm))[0]
    # (ell, m) = (2, +/- 0)
    Symm = Vars.rhOverM_coeff*(Cons.hHat_2_0_0 + Vars.v**2*(Cons.hHat_2_0_2 + Vars.v*(Vars.hHat_spin_Symm_2_0_3 + Vars.v*(Cons.hHat_2_0_4 + Vars.hHat_spin_Symm_2_0_4 + Vars.v*(Cons.hHat_2_0_5 + Vars.hHat_spin_Symm_2_0_5 + Vars.v*(Cons.hHat_2_0_6 + Vars.hHat_spin_Symm_2_0_6 + Vars.hHat_spin_Symm_2_0_7*Vars.v))))))
    Asymm = 0
    modes[6] = (Symm + Asymm)[0]
    # (ell, m) = (2, +/- 1)
    Symm = 0
    Asymm = Vars.rhOverM_coeff*Vars.v**3*(Vars.hHat_spin_Asymm_2_1_3 + Vars.hHat_spin_Asymm_2_1_4*Vars.v)
    modes[7] = (Symm + Asymm)[0]
    modes[5] = (np.conjugate(Symm - Asymm))[0]
    # (ell, m) = (3, +/- 1)
    Symm = Vars.rhOverM_coeff*Vars.v**6*(Cons.hHat_3_1_6 + Vars.hHat_spin_Symm_3_1_7*Vars.v)
    Asymm = 0
    modes[13] = (Symm + Asymm)[0]
    modes[11] = (-np.conjugate(Symm - Asymm))[0]
    # (ell, m) = (3, +/- 3)
    Symm = Vars.rhOverM_coeff*Vars.v**6*(Cons.hHat_3_3_6 + Vars.hHat_spin_Symm_3_3_7*Vars.v)
    Asymm = 0
    modes[15] = (Symm + Asymm)[0]
    modes[9] = (-np.conjugate(Symm - Asymm))[0]
    # (ell, m) = (4, +/- 0)
    Symm = Vars.rhOverM_coeff*(Cons.hHat_4_0_0 + Vars.v**2*(Cons.hHat_4_0_2 + Vars.v*(Vars.hHat_spin_Symm_4_0_3 + Vars.v*(Cons.hHat_4_0_4 + Vars.hHat_spin_Symm_4_0_4 + Vars.v*(Cons.hHat_4_0_5 + Vars.hHat_spin_Symm_4_0_5 + Vars.v*(Cons.hHat_4_0_6 + Vars.hHat_spin_Symm_4_0_6 + Vars.hHat_spin_Symm_4_0_7*Vars.v))))))
    Asymm = 0
    modes[20] = (Symm + Asymm)[0]
    # (ell, m) = (4, +/- 1)
    Symm = 0
    Asymm = Vars.rhOverM_coeff*Vars.v**3*(Vars.hHat_spin_Asymm_4_1_3 + Vars.hHat_spin_Asymm_4_1_4*Vars.v)
    modes[21] = (Symm + Asymm)[0]
    modes[19] = (np.conjugate(Symm - Asymm))[0]
    # (ell, m) = (4, +/- 4)
    Symm = Cons.hHat_4_4_5*Vars.rhOverM_coeff*Vars.v**5
    Asymm = 0
    modes[24] = (Symm + Asymm)[0]
    modes[16] = (np.conjugate(Symm - Asymm))[0]
    # (ell, m) = (5, +/- 1)
    Symm = Cons.hHat_5_1_6*Vars.rhOverM_coeff*Vars.v**6
    Asymm = 0
    modes[31] = (Symm + Asymm)[0]
    modes[29] = (-np.conjugate(Symm - Asymm))[0]
    # (ell, m) = (5, +/- 3)
    Symm = Cons.hHat_5_3_6*Vars.rhOverM_coeff*Vars.v**6
    Asymm = 0
    modes[33] = (Symm + Asymm)[0]
    modes[27] = (-np.conjugate(Symm - Asymm))[0]
    # (ell, m) = (5, +/- 5)
    Symm = Cons.hHat_5_5_6*Vars.rhOverM_coeff*Vars.v**6
    Asymm = 0
    modes[35] = (Symm + Asymm)[0]
    modes[25] = (-np.conjugate(Symm - Asymm))[0]
    # (ell, m) = (6, +/- 0)
    Symm = Vars.rhOverM_coeff*(Cons.hHat_6_0_0 + Vars.v**2*(Cons.hHat_6_0_2 + Vars.v**2*(Cons.hHat_6_0_4 + Vars.v*(Cons.hHat_6_0_5 + Vars.hHat_spin_Symm_6_0_5 + Vars.v*(Cons.hHat_6_0_6 + Vars.hHat_spin_Symm_6_0_6 + Vars.hHat_spin_Symm_6_0_7*Vars.v)))))
    Asymm = 0
    modes[42] = (Symm + Asymm)[0]
    # (ell, m) = (8, +/- 0)
    Symm = Vars.rhOverM_coeff*(Cons.hHat_8_0_0 + Vars.v**2*(Cons.hHat_8_0_2 + Vars.v**2*(Cons.hHat_8_0_4 + Vars.v**2*(Cons.hHat_8_0_6 + Vars.hHat_spin_Symm_8_0_7*Vars.v))))
    Asymm = 0
    modes[72] = (Symm + Asymm)[0]

def Modes(wHat_i, xHat_i, yHat_i, zHat_i, M1_i, M2_i, v_i, S_chi1_i, S_chi2_i, frame, y, PNWaveformModeOrder=3.5) :
    Recalculate={        0:Recalculate_0,
        1:Recalculate_0p50,
        2:Recalculate_1p0,
        3:Recalculate_1p5,
        4:Recalculate_2p0,
        5:Recalculate_2p5,
        6:Recalculate_3p0,
        7:Recalculate_3p5,
        8:Recalculate_4p0}
    WaveformModes={        0:Modes_0,
        1:Modes_0p50,
        2:Modes_1p0,
        3:Modes_1p5,
        4:Modes_2p0,
        5:Modes_2p5,
        6:Modes_3p0,
        7:Modes_3p5,
        8:Modes_4p0}
    ModeData=np.zeros((len(y[0]),81), dtype=complex)
    zz=np.array([0.0])
    z=np.array([0.0*1j])
    Cons=Constants(zz,zz,zz,zz,zz,zz,zz,zz,zz,zz,zz,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z)
    Vars=Variables(zz,zz,zz,zz,zz,zz,zz,zz,zz,zz,zz,zz,zz,zz,zz,zz,zz,zz,zz,zz,zz,zz,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z)
    Initialization(Cons, wHat_i, xHat_i, yHat_i, zHat_i, M1_i, M2_i, v_i, S_chi1_i, S_chi2_i)
    for i in range(len(y[0])):
        Recalculate.get(2*PNWaveformModeOrder)(Cons,Vars,y[:,i])
        WaveformModes.get(2*PNWaveformModeOrder)(ModeData[i,:],Cons,Vars)
    return ModeData, [0,8]
