#
# Copyright (C) 2012  Leo Singer
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
from __future__ import division
"""
Miscellaneous routines related to sky localization.
"""
__author__ = "Leo Singer <leo.singer@ligo.org>"


import numpy as np
import lal, lalsimulation


def mchirp(m1=1.4, m2=1.4):
    """Find chirp mass from component masses."""
    return (m1 * m2) ** 0.6 * (m1 + m2) ** -0.2


def get_noise_psd_func(ifo):
    """Find a function that describes the given interferometer's noise PSD."""
    if ifo in ("H1", "H2", "L1", "I1"):
        func = lalsimulation.SimNoisePSDaLIGOZeroDetHighPower
    elif ifo == "V1":
        func = lalsimulation.SimNoisePSDAdvVirgo
    elif ifo == "K1":
        func = lalsimulation.SimNoisePSDKAGRA
    else:
        raise ValueError("Unknown interferometer: %s", ifo)
    return func


def get_noise_integral(ifo, power, mass1, mass2, order=7, f_low=10, norm=False):
    """Calculate the integral of 4 |h(f)|^2 * f^power / S(f) from f=f_low
    to f_LSO, where h(f) is a post-Newtonian signal of the given amplitude
    order for an optimally oriented inspiral source at a distance of 1 Mpc,
    and S(f) is the noise power spectral density for the interferometer whose
    two-letter site name is given by ifo. If norm=True, then divide the
    result by the integral of 4 |h(f)|^2 / S(f) before returning."""

    # Integration step in Hz.
    df = 1.

    # Frequency-domain post-Newtonian inspiral waveform.
    h = lalsimulation.SimInspiralTaylorF2(0, df,
        mass1 * lal.LAL_MSUN_SI, mass2 * lal.LAL_MSUN_SI,
        f_low, 1e6 * lal.LAL_PC_SI, 0, order)

    # Find indices of first and last nonzero samples.
    nonzero = np.nonzero(h.data.data)[0]
    first_nonzero = nonzero[0]
    last_nonzero = nonzero[-1]

    # Frequency sample points
    f = h.f0 + h.deltaF * np.arange(first_nonzero, last_nonzero + 1)

    # Throw away leading and trailing zeros.
    h = h.data.data[first_nonzero:last_nonzero + 1]

    # Noise PSD function.
    S = get_noise_psd_func(ifo)
    S = [S(ff) for ff in f]

    denom_integrand = 4 * (np.square(h.real) + np.square(h.imag)) / S
    num_integrand = denom_integrand * f ** power

    ret = np.trapz(num_integrand, dx=df)

    if norm:
        ret /= np.trapz(denom_integrand, dx=df)
    return ret


def get_noise_moment(ifo, power, mass1, mass2, order=7, f_low=10):
    return get_noise_integral(ifo, power, mass1, mass2, order, f_low, True)


def get_horizon_distance(ifo, mass1, mass2, order=7, f_low=10, snr_thresh=1):
    """Compute the distance at which a source would produce a maximum SNR of
    snr_thresh in the given interferometer."""

    return np.sqrt(get_noise_integral(ifo, 0, mass1, mass2, order, f_low)) / snr_thresh


def get_effective_bandwidth(ifo, mass1, mass2, order=7, f_low=10):
    return np.sqrt(get_noise_moment(ifo, 2, mass1, mass2, order, f_low) - np.square(get_noise_moment(ifo, 1, mass1, mass2, order, f_low)))


def get_f_lso(mass1, mass2):
    """Calculate the GW frequency during the last stable orbit of a compact binary."""
    return 1 / (6 ** 1.5 * np.pi * (mass1 + mass2) * lal.LAL_MTSUN_SI)


def toa_uncertainty(snr, bandwidth):
    """Estimate uncertainty in TOA measurements using Cramer-Rao bound."""
    return 1 / (2 * np.pi * abs(snr) * bandwidth)
