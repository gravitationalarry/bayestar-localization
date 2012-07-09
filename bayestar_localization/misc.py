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
import swiglal, swiglalsimulation


def mchirp(m1=1.4, m2=1.4):
    """Find chirp mass from component masses."""
    return (m1 * m2) ** 0.6 * (m1 + m2) ** -0.2


def get_noise_psd_func(ifo):
    """Find a function that describes the given interferometer's noise PSD."""
    if ifo in ("H1", "H2", "L1", "I1"):
        func = swiglalsimulation.XLALSimNoisePSDaLIGOZeroDetHighPower
    elif ifo == "V1":
        func = swiglalsimulation.XLALSimNoisePSDAdvVirgo
    elif ifo == "K1":
        func = swiglalsimulation.XLALSimNoisePSDKAGRA
    else:
        raise ValueError("Unknown interferometer: %s", ifo)
    return func


def get_horizon_distance(ifo, mchirp=mchirp(1.4, 1.4), f_low=10., f_high=1560., snr_thresh=1):
    """Compute the distance at which a source would produce a maximum SNR of
    snr_thresh in the given interferometer."""

    # Chirp mass in seconds.
    tchirp = swiglal.LAL_MTSUN_SI * mchirp

    # Integration step in Hz.
    df = 1.

    # Evaluation frequencies.
    f = np.linspace(f_low, f_high, (f_high - f_low) / df)

    # Noise PSD function.
    S = get_noise_psd_func(ifo)
    S = [S(ff) for ff in f]

    Dhor = np.sum(f ** (-7/3) / S * df)
    Dhor = np.sqrt(5/6 * Dhor)
    Dhor *= tchirp ** (5/6) * np.pi ** (-2/3) * swiglal.LAL_C_SI / snr_thresh
    Dhor /= 1e6 * swiglal.LAL_PC_SI
    return Dhor


def get_effective_bandwidth(ifo, f_low=10., f_high=1560.):
    # Integration step in Hz.
    df = 1.

    # Evaluation frequencies.
    f = np.linspace(f_low, f_high, (f_high - f_low) / df)

    # Noise PSD function.
    S = get_noise_psd_func(ifo)
    S = [S(ff) for ff in f]

    return np.sqrt(np.sum(f ** (-1/3) / S) / np.sum(f ** (-7/3) / S))


def get_f_lso(mass1, mass2):
    """Calculate the GW frequency during the last stable orbit of a compact binary."""
    return 1 / (6 ** 1.5 * np.pi * (mass1 + mass2) * swiglal.LAL_MTSUN_SI)


def toa_uncertainty(snr, bandwidth):
    """Estimate uncertainty in TOA measurements using Cramer-Rao bound."""
    return 1 / (2 * np.pi * abs(snr) * bandwidth)
