import scri
from . import PNEvolution
from . import PNWaveformModes
from ..waveforms import WaveformModes
import numpy as np
import quaternionic

def PNWaveform(q,M,omega_0,chi1_0,chi2_0,frame_0,t_0=0.0, t_PNStart=False, t_PNEnd=False, return_chi=False, PNEvolutionOrder=4.0, PNWaveformModeOrder=4.0, TaylorTn=1, StepsPerOrbit=32, ForwardInTime=True, tol=1e-10, MinStep=1e-7):
    """
    q = m1/m2, float number,
    M = m1+m2, float number,
    omega_0: magnititude of angular velocity at t_0, float number,
    chi1_0 and chi2_0: spin vectors at t_0, 3-d vectors,
    frame_0: the frame quaternion at t_0, quaternionic_array object,
    t_0: the corresponding time of the above given initial values, float number,
    t_PNStart: the start time of PN relative to t_0: t_PNStart=t_real_start-t_0, float number. If false, default is t_0-t_real_start=3(t_merger-t_0),
    t_PNEnd: the end time of PN relative to t_0: t_PNEnd=t_real_end-t_0, float number. If false, default is merger time,
    return_chi: whether to return chi as quaternionic array of time, bool number,
    PNEvolutionOrder: float number in [0,0.5,1,1.5,2,2.5,3,3.5,4], default is 3.5,
    PNWaveformModeOrder: float number in [0,0.5,1,1.5,2,2.5,3,3.5,4], default is 3.5,
    TaylorTn: now only TaylorT1 is working, so its int number in [1], default is 1,
    StepsPerOrbit: float number,
    ForwardInTime: whether to evlove PN forward in time, bool number, default if True,
    tol: tolerance of the integrator, float number,
    MinStep: minimal time interval for the PN waveform, float number.
    """

    if not PNEvolutionOrder in [0,0.5,1,1.5,2,2.5,3,3.5,4]:
        message=("PNEvolutionOrder must be a float number in [0,0.5,1,1.5,2,2.5,3,3.5,4].")
        raise ValueError(message)
    if not PNWaveformModeOrder in [0,0.5,1,1.5,2,2.5,3,3.5,4]:
        message=("PNWaveformModeOrder must be a float number in [0,0.5,1,1.5,2,2.5,3,3.5,4].")
        raise ValueError(message)
    if not TaylorTn in [1]:
        message=("TaylorTn must be an int number in [1].")
        raise ValueError(message)          

    xHat=quaternionic.x
    yHat=quaternionic.y
    zHat=quaternionic.z
    m1=q/(1+q)*M
    m2=1/(1+q)*M
    v_0=omega_0**(1/3)
    chi1Mag=quaternionic.array([0,chi1_0[0],chi1_0[1],chi1_0[2]]).abs
    chi2Mag=quaternionic.array([0,chi2_0[0],chi2_0[1],chi2_0[2]]).abs

    # Quaternions that rotate z-axis to spin vectors
    S_chi1_0=0.0*quaternionic.one 
    S_chi2_0=0.0*quaternionic.one
    if chi1Mag>1e-12:
        S_chi1_0=np.sqrt(chi1Mag)*np.sqrt(
            -quaternionic.array([0,chi1_0[0],chi1_0[1],chi1_0[2]]).normalized*zHat).normalized
    if chi2Mag>1e-12:
        S_chi2_0=np.sqrt(chi2Mag)*np.sqrt(
            -quaternionic.array([0,chi2_0[0],chi2_0[1],chi2_0[2]]).normalized*zHat).normalized
        
    rfrak_frame_0=np.log(frame_0).vector # logarithm of frame quaternion  
    PN=PNEvolution.PNEv.Evolution(xHat, yHat, zHat, m1, m2, v_0,S_chi1_0, S_chi2_0, rfrak_frame_0, t_PNStart, t_PNEnd,
        PNEvolutionOrder, TaylorTn, StepsPerOrbit, ForwardInTime, tol, MinStep)# Evolve PN parameters, PN.t is PN time, PN.y=[v, chi1_x, chi1_y
        # chi2_x, chi2_y, rfrak_frame_x, rfrak_frame_y, rfrak_frame_z]

    data, [ellmin,ellmax] = PNWaveformModes.Modes(xHat, yHat, zHat, m1, m2, v_0,S_chi1_0, S_chi2_0, rfrak_frame_0, PN.y, PNWaveformModeOrder)
    W_PN_corot=WaveformModes(data,
        time=PN.t+t_0,
        frame=np.exp(quaternionic.array(np.column_stack((0.0*PN.t,PN.y[5],PN.y[6],PN.y[7])))),
        modes_axis=1,
        ell_min=ellmin,
        ell_max=ellmax,
        frame_type="corotating",
        data_type="h",
        spin_weight=-2)

    if return_chi:
        R_S1=np.exp(PN.y[1]*xHat + PN.y[2]*yHat)
        R_S2=np.exp(PN.y[3]*xHat + PN.y[4]*yHat)
        chiVec1=S_chi1_0*R_S1*zHat*R_S1.conjugate()*S_chi1_0.conjugate()
        chiVec2=S_chi2_0*R_S2*zHat*R_S2.conjugate()*S_chi2_0.conjugate()
        return W_PN_corot, chiVec1, chiVec2
    else:
        return W_PN_corot
