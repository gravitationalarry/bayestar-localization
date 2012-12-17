#!/usr/bin/env python
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
Test generation of colored noise.
"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import mlab
import lal
import lalsimulation
from bayestar import filter
from bayestar import timing

# Sample rate in Hz
sample_rate = 4096

# Data duration in seconds
data_duration = 100

# Number of samples
data_length = data_duration * sample_rate

# Start time of data
epoch = lal.LIGOTimeGPS()

# Discretely sampled power spectrum
psd = lal.CreateREAL8FrequencySeries(None, epoch, 0, 1 / data_duration, filter.unitInverseHertz, data_length // 2 + 1)
S = timing.get_noise_psd_func("H1")
f = filter.abscissa(psd)
psd.data.data = [S(ff) for ff in f]
psd.data.data[(f < 10) | (f > 2048)] = 0

# Generate colored Gaussian noise
data = filter.colored_noise(lal.LIGOTimeGPS(), data_duration, sample_rate, psd)

# Plot measured power spectrum of data
[power, freqs] = plt.psd(data.data.data, NFFT=sample_rate, Fs=sample_rate, detrend=mlab.detrend_mean, noverlap=sample_rate//2)
plt.clf()
plt.loglog(freqs, power)

# Plot goal PSD
plt.loglog(f, psd.data.data, label='goal')

plt.xlabel('frequency (Hz)')
plt.ylabel(r'power spectral density (Hz$^{-1}$)')
plt.grid()
plt.legend()

plt.savefig('psd.png')
