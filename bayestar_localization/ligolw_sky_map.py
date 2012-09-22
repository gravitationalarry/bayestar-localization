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
"""
Convenience function to produce a sky map from LIGO-LW rows.
"""
from __future__ import division
__author__ = "Leo Singer <leo.singer@ligo.org>"


import time
import numpy as np
import bayestar_localization
from bayestar_localization import timing
from pylal import date
import lal, lalsimulation


def complex_from_polar(r, phi):
    return r * (np.cos(phi) + np.sin(phi) * 1j)


def ligolw_sky_map(sngl_inspirals, order, f_low, min_distance=None, max_distance=None, prior=None, method="toa_snr", reference_frequency=None, fd=False, nside=-1):
    """Convenience function to produce a sky map from LIGO-LW rows. Note that
    min_distance and max_distance should be in Mpc."""

    if method == "toa_snr" and (min_distance is None or max_distance is None or prior is None):
        raise ValueError("For method='toa_snr', the arguments min_distance, max_distance, and prior are required.")

    ifos = [sngl_inspiral.ifo for sngl_inspiral in sngl_inspirals]

    # Extract masses from the table.
    mass1s = np.asarray([sngl_inspiral.mass1 for sngl_inspiral in sngl_inspirals])
    mass2s = np.asarray([sngl_inspiral.mass2 for sngl_inspiral in sngl_inspirals])

    # Extract SNRs from table.
    snrs = np.asarray([complex_from_polar(sngl_inspiral.snr, sngl_inspiral.coa_phase)
        for sngl_inspiral in sngl_inspirals])

    # Extract TOAs from table.
    toas_ns = np.asarray([sngl_inspiral.get_end().ns()
        for sngl_inspiral in sngl_inspirals])

    # Optionally apply reference frequency shift.
    if reference_frequency is not None:
        toas_ns -= [int(round(1e9 * lalsimulation.
            SimInspiralTaylorF2ReducedSpinChirpTime(
            reference_frequency,
            m1 * lal.LAL_MSUN_SI,
            m2 * lal.LAL_MSUN_SI,
            0, 4))) for m1, m2 in zip(mass1s, mass2s)]

    # Convert TOAs from nanoseconds to seconds.
    toas = 1e-9 * toas_ns

    # Find average Greenwich mean sidereal time of event.
    mean_toa_ns = sum(toas_ns) // len(toas_ns)
    epoch = date.XLALINT8NSToGPS(mean_toa_ns)
    gmst = date.XLALGreenwichMeanSiderealTime(epoch)

    # Signal models for each detector.
    signal_models = [timing.SignalModel(mass1, mass2, timing.get_noise_psd_func(ifo), order=order, fd=fd)
        for mass1, mass2, ifo in zip(mass1s, mass2s, ifos)]

    # Get SNR=1 horizon distances for each detector.
    horizons = [signal_model.get_horizon_distance()
        for signal_model in signal_models]

    # Estimate TOA uncertainty (squared) using CRB or BRB evaluated at MEASURED
    # values of the SNRs.
    s2_toas = [np.square(signal_model.get_toa_uncert(np.abs(snr)))
        for signal_model, snr in zip(signal_models, snrs)]

    # Look up physical parameters for detector.
    detectors = [lalsimulation.InstrumentNameToLALDetector(str(ifo))
        for ifo in ifos]
    responses = [det.response for det in detectors]
    locations = [det.location for det in detectors]

    # Time and run sky localization.
    start_time = time.time()
    if method == "toa":
        sky_map = bayestar_localization.sky_map_tdoa(gmst, toas, s2_toas, locations, nside=nside)
    elif method == "toa_snr":
        sky_map = bayestar_localization.sky_map_tdoa_snr(gmst, toas, snrs, s2_toas, responses, locations, horizons, min_distance, max_distance, prior, nside=nside)
    else:
        raise ValueError("Unrecognized method: %s" % method)
    end_time = time.time()

    # Find elapsed run time.
    elapsed_time = end_time - start_time

    # Done!
    return sky_map, epoch, elapsed_time
