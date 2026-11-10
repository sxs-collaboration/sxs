from gwtemplate.waveform import CompactBinaryCoalescenceGenerator
from gwtemplate.mass_time import MassTime
from gwtemplate import constants
import astropy.units as u
from astropy.constants import G, c ,M_sun
import numpy as np
from .. import load, TimeSeries, WaveformModes
from scipy.interpolate import CubicSpline

MTSUN_SI = constants.get_constant("MTSUN_SI")

class SXS(CompactBinaryCoalescenceGenerator):

    def __init__(self, **kwargs):
        super(SXS, self).__init__()

    @property
    def epoch(self):
        return self._epoch
    @epoch.setter
    def epoch(self, value):
        self._epoch = value

    @property
    def duration(self):
        return self._duration
    @duration.setter
    def duration(self, value):
        self._duration = value
        
    
        
    def initialize_model(self, **parameters):
        """Loads an SXS simulation, reads the parameters at the reference frequency, and 
        converts to the appropriate conventions for LISA. 
        Parameters should include only extrinsic parameters and properties of the waveform
        Intrinsic parameters will be filled in from the NR metadata."""

        self.simulation = load(parameters['SXS_ID'])
        metadata = self.simulation.metadata
        t_scaling = (parameters['total_mass'] / u.solMass * MTSUN_SI)
        times,_, dynamics_dict = self.simulation.to_lvk(
            f_ref = (t_scaling * parameters['f22_ref']).value / 2, f_low = None)
        parameters['mass_ratio'] = metadata['reference_mass_ratio']

        ###The vector between black holes is reversed in the LISA conventions
        parameters['spin1x'] = -dynamics_dict['chi1_ref'][0][0]
        parameters['spin1y'] = -dynamics_dict['chi1_ref'][0][1]
        parameters['spin1z'] =  dynamics_dict['chi1_ref'][0][2]
        
        parameters['spin2x'] = -dynamics_dict['chi2_ref'][0][0]
        parameters['spin2y'] = -dynamics_dict['chi2_ref'][0][1]
        parameters['spin2z'] =  dynamics_dict['chi2_ref'][0][2]

        self.params = parameters

        self.epoch = dynamics_dict['t_low'] * t_scaling
        self.duration = (times[-1] - times[0]) * t_scaling

    def generate_time_domain_modes(self, times = None, deltaT = None):
        """
        Produce time domain modes from the loaded SXS waveform

        Parameters
        ----------

            Either an array of times, or a time step deltaT

        Returns
        -------

            sxs.WaveformModes containing the modes, interpolated to the desired time array

        Notes
        -----

            At least one of the ``times`` or ``deltaT`` must be provided.

        """
        if times is None and deltaT is None:
            raise ValueError("At least one of times and deltaT must be provided")
        t_scaling = (self.params['total_mass'] / u.solMass * MTSUN_SI)
        h_scaling = G * (self.params['total_mass']).to('kg') / (c**2 * self.params['distance'].to('m'))
        t, hlm, _ = self.simulation.to_lvk(
                f_ref = (t_scaling * self.params['f22_ref']).value / 2, f_low = None, 
                dt = ( deltaT / t_scaling).value if deltaT is not None else None)
        t *= t_scaling.value
        h_array = []
        for l in range(2, self.simulation.h.ell_max + 1):
            for m in range(-l, l + 1):
                h_array.append(hlm[(l,m)])

        hlm = WaveformModes(time = t,
                          input_array = h_scaling * np.array(h_array),
                          modes_axis = 0,
                          ell_min = 2,
                          ell_max = self.simulation.h.ell_max)
        if times is not None:
            hlm = hlm.interpolate(times)
            self.duration = (times[-1] - times[0])
        return hlm
    
    def generate_time_domain_polarizations(self, times = None, deltaT = None):
        """
        Generate time domain polarizations from the loaded SXS waveform.

        Parameters
        ----------

            Either an array of times, or a time step deltaT

        Returns
        -------

            sxs.TimeSeries containing the plus and cross polarizations

        Notes
        -----

            At least one of the ``times`` or ``deltaT`` must be provided.

        """
        if times is None and deltaT is None:
            raise ValueError("At least one of times and deltaT must be provided")
        t_scaling = (self.params['total_mass'] / u.solMass * MTSUN_SI)
        h_scaling = G * (self.params['total_mass']).to('kg') / (c**2 * self.params['distance'].to('m'))
        t, hp, hc, _ = self.simulation.to_lvk(
                f_ref = (t_scaling * self.params['f22_ref']).value / 2, f_low = None, 
                dt = (deltaT / t_scaling).value if deltaT is not None else None,
                inclination = self.params['inclination'].value, phi_ref = (np.pi / 2 * u.rad - self.params['phi_ref'].to(u.rad)).value)
        ##Convert to SI units
        t *= t_scaling.value
        hp *= h_scaling.value
        hc *= h_scaling.value
        hp = TimeSeries(input_array = hp, time = t)
        hc = TimeSeries(input_array = hc, time = t)
        if times is not None:
            hp = hp.interpolate(times)
            hc = hc.interpolate(times)
            self.duration = (times[-1] - times[0])

        return hp, hc
    
    def generate_frequency_domain_polarizations(self, **kwargs):
        """
        Generate fourier domain polarizations from the loaded SXS waveform.

    
        Returns
        -------

            sxs.TimeSeries containing the frequency domain plus and cross polarizations
            the frequency array is in the time attribute of the returned TimeSeries
        """
        t_scaling = (self.params['total_mass'] / u.solMass * MTSUN_SI)
        t, hp, hc, _ = self.simulation.to_lvk(
                f_ref = (t_scaling * self.params['f22_ref']).value / 2, f_low = None, 
                inclination = self.params['inclination'].value, phi_ref = (np.pi / 2 * u.rad - self.params['phi_ref'].to(u.rad)).value)
        ##Preprocess and fourier transform. This also rescales to SI units
        hp = TimeSeries(time = t, input_array = hp)
        hp = hp.preprocess()
        hp = hp.fourier_transform(total_mass = self.params['total_mass'], luminosity_distance=self.params['distance'].to(u.Mpc).value)
        
        hc = TimeSeries(time = t, input_array = hc)
        hc = hc.preprocess()
        hc = hc.fourier_transform(total_mass = self.params['total_mass'], luminosity_distance=self.params['distance'].to(u.Mpc).value)
        return hp, hc
        

    def generate_frequency_domain_modes(self, **kwargs):
        """Generate fourier domain modes from the loaded SXS waveform.
        Returns
        -------

            sxs.WaveformModes containing the frequency domain plus and cross polarizations
            the frequency array is in the time attribute of the returned WaveformModes
        """
        t_scaling = (self.params['total_mass'] / u.solMass * MTSUN_SI)
        t, hlm, _ = self.simulation.to_lvk(
                f_ref = (t_scaling * self.params['f22_ref']).value / 2, f_low = None, )
        h_array = []
        for l in range(2, self.simulation.h.ell_max + 1):
            for m in range(-l, l + 1):
                h_array.append(hlm[(l,m)])
        ##Fourier transform
        h = WaveformModes(time = t,
                          input_array = h_array,
                          modes_axis = 1,
                          ell_min = 2,
                          ell_max = self.simulation.h.ell_max)
        h = h.preprocess()
        h = h.fourier_transform(total_mass =self.params['total_mass'], luminosity_distance=self.params['distance'].to(u.Mpc).value)
        return h

    def uniform_time_array(self, deltaT):
        time_array = np.arange(0.0, self.duration.value, deltaT.value) + self.epoch.value
        time_array *= (self.params['total_mass'].value * MTSUN_SI)**(-1)
        return MassTime(time_array, self.params['total_mass'])
        



