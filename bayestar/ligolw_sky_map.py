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
from . import timing
from . import sky_map
from pylal import date
import lal, lalsimulation


# Copied and adapted from gstlal.reference_psd.read_psd_xmldoc.
# FIXME: remove this when read_psd_xmldoc is available in releases of gstlal
# that are installed on all clusters.
def read_psd_xmldoc(xmldoc):
    """
    Parse a dictionary of PSD frequency series objects from an XML
    document.  See also make_psd_xmldoc() for the construction of XML documents
    from a dictionary of PSDs.  Interprets an empty freuency series for an
    instrument as None.
    """
    from glue.ligolw import ligolw
    from glue.ligolw import param
    from pylal import datatypes as laltypes
    from pylal import series as lalseries

    out = dict((param.get_pyvalue(elem, u"instrument"), lalseries.parse_REAL8FrequencySeries(elem)) for elem in xmldoc.getElementsByTagName(ligolw.LIGO_LW.tagName) if elem.hasAttribute(u"Name") and elem.getAttribute(u"Name") == u"REAL8FrequencySeries")
    # Interpret empty frequency series as None
    for k in out:
        if len(out[k].data) == 0:
            out[k] = None
    return out
# End section copied and adapted from gstlal.reference_psd.read_psd_xmldoc.


def complex_from_polar(r, phi):
    return r * (np.cos(phi) + np.sin(phi) * 1j)


def ligolw_sky_map(sngl_inspirals, approximant, amplitude_order, phase_order, f_low, min_distance=None, max_distance=None, prior=None, method="toa_snr", reference_frequency=None, psds=None, nside=-1):
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

    # Power spectra for each detector.
    if psds is None:
        psds = [timing.get_noise_psd_func(ifo) for ifo in ifos]

    # Signal models for each detector.
    signal_models = [timing.SignalModel(mass1, mass2, psd, f_low, approximant, amplitude_order, phase_order)
        for mass1, mass2, psd in zip(mass1s, mass2s, psds)]

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
        prob = sky_map.tdoa(gmst, toas, s2_toas, locations, nside=nside)
    elif method == "toa_snr":
        prob = sky_map.tdoa_snr(gmst, toas, snrs, s2_toas, responses, locations, horizons, min_distance, max_distance, prior, nside=nside)
    else:
        raise ValueError("Unrecognized method: %s" % method)
    end_time = time.time()

    # Find elapsed run time.
    elapsed_time = end_time - start_time

    # Done!
    return prob, epoch, elapsed_time


def gracedb_sky_map(coinc_file, psd_file, waveform, f_low, min_distance, max_distance, prior, reference_frequency=None, nside=-1):
    # LIGO-LW XML imports.
    from glue.ligolw import table as ligolw_table
    from glue.ligolw import utils as ligolw_utils
    from glue.ligolw import lsctables

    # gstlal imports
    from gstlal import reference_psd

    # BAYESTAR imports.
    from bayestar import io

    # Determine approximant, amplitude order, and phase order from command line arguments.
    approximant, amplitude_order, phase_order = timing.get_approximant_and_orders_from_string(waveform)

    # Read input file.
    xmldoc, _ = ligolw_utils.load_fileobj(coinc_file)

    # Locate the tables that we need.
    coinc_inspiral_table = ligolw_table.get_table(xmldoc,
        lsctables.CoincInspiralTable.tableName)
    coinc_map_table = ligolw_table.get_table(xmldoc,
        lsctables.CoincMapTable.tableName)
    sngl_inspiral_table = ligolw_table.get_table(xmldoc,
        lsctables.SnglInspiralTable.tableName)

    # Locate the sngl_inspiral rows that we need.
    coinc_inspiral = coinc_inspiral_table[0]
    coinc_event_id = coinc_inspiral.coinc_event_id
    event_ids = [coinc_map.event_id for coinc_map in coinc_map_table
        if coinc_map.coinc_event_id == coinc_event_id]
    sngl_inspirals = [(sngl_inspiral for sngl_inspiral in sngl_inspiral_table
        if sngl_inspiral.event_id == event_id).next() for event_id in event_ids]

    # Read PSDs.
    xmldoc, _ = ligolw_utils.load_fileobj(psd_file)
    psds = read_psd_xmldoc(xmldoc)

    # Rearrange PSDs into the same order as the sngl_inspirals.
    psds = [psds[sngl_inspiral.ifo] for sngl_inspiral in sngl_inspirals]

    # Interpolate PSDs.
    psds = [timing.interpolate_psd(psd.f0 + np.arange(len(psd.data)) * psd.deltaF, psd.data) for psd in psds]

    # TOA+SNR sky localization
    return ligolw_sky_map(sngl_inspirals, approximant, amplitude_order, phase_order, f_low,
        min_distance, max_distance, prior,
        reference_frequency=reference_frequency, nside=nside, psds=psds)
