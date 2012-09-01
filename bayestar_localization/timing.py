# -*- coding: utf-8
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
Functions for predicting timing accuracy of matched filters.
"""
__author__ = "Leo Singer <leo.singer@ligo.org>"


import numpy as np
import lal, lalsimulation
from scipy import interpolate
from scipy import linalg
from scipy import special


def mchirp(m1=1.4, m2=1.4):
    """Find chirp mass from component masses."""
    return (m1 * m2) ** 0.6 * (m1 + m2) ** -0.2


def get_f_lso(mass1, mass2):
    """Calculate the GW frequency during the last stable orbit of a compact binary."""
    return 1 / (6 ** 1.5 * np.pi * (mass1 + mass2) * lal.LAL_MTSUN_SI)


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


def interpolate_psd(f, S):
    """Create a (linear in log-log) interpolating function for a discretely
    sampled power spectrum S(f)."""
    func = interpolate.interp1d(np.log(f), np.log(S))
    return (lambda f: np.exp(func(np.log(f))))


class SignalModel(object):
    """Class to speed up computation of signal/noise-weighted integrals and
    Barankin and Cramer-Rao lower bounds on time and phase estimation."""

    def __init__(self, mass1, mass2, S, order=7, f_low=10):
        """Create a TaylorF2 signal model with the given masses, PSD function
        S(f), PN amplitude order, and low-frequency cutoff."""

        # Integration step in Hz.
        df = 1
        self.dw = 2 * np.pi * df

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
        self.w = 2 * np.pi * f

        # Throw away leading and trailing zeros.
        h = h.data.data[first_nonzero:last_nonzero + 1]

        # Noise PSD function.
        S = [S(ff) for ff in f]

        self.denom_integrand = 4 / (2 * np.pi) * (np.square(h.real) + np.square(h.imag)) / S
        self.den = np.trapz(self.denom_integrand, dx=self.dw)

    def get_horizon_distance(self, snr_thresh=1):
        return np.sqrt(self.den) / snr_thresh

    def get_sn_average(self, func):
        """Get the average of a function of angular frequency, weighted by the
        signal to noise per unit angular frequency."""
        return np.trapz(func(self.w) * self.denom_integrand, dx=self.dw) / self.den

    def get_sn_moment(self, power):
        """Get the average of angular frequency to the given power, weighted by
        the signal to noise per unit frequency."""
        return self.get_sn_average(lambda w: w**power)

    def get_crb(self, snr):
        """Get the Cramer-Rao bound, or inverse Fisher information matrix,
        describing the phase and time estimation covariance."""
        w1 = self.get_sn_moment(1)
        w2 = self.get_sn_moment(2)
        I = np.asarray(((1, -w1), (-w1, w2)))
        return linalg.inv(I) / np.square(snr)

    def get_brb(self, snr):
        """Get the Barankin bound on the phase and time estimation covariance."""

        # Mean-square angular frequency.
        w1 = self.get_sn_moment(1)
        w2 = self.get_sn_moment(2)

        # FIXME: stray factor of pi from somewhere?
        snr = snr / np.pi

        # Fisher information (factor of snr**2 kept separate)
        snr2 = np.square(snr)
        Lambda = np.asmatrix([[1, -w1], [-w1, w2]])

        # Test values in S/N, phase, and time.
        dphi = np.pi / 9
        dtau = np.pi / 4 / np.sqrt(w2)
        n = 8
        i = np.arange(1, 2 * n + 1)

        sinw_moments = np.concatenate(([0], [self.get_sn_average(lambda w: np.sin(w * ii * dtau)) for ii in i]))
        cosw_moments = np.concatenate(([1], [self.get_sn_average(lambda w: np.cos(w * ii * dtau)) for ii in i]))
        wsinw_moments = np.concatenate(([0], [self.get_sn_average(lambda w: w * np.sin(w * ii * dtau)) for ii in i]))
        wcosw_moments = np.concatenate(([w1], [self.get_sn_average(lambda w: w * np.cos(w * ii * dtau)) for ii in i]))

        i = np.concatenate((np.arange(-n, 0), np.arange(1, n + 1)))
        iphi, itau = np.meshgrid(i, i)
        iphi = np.atleast_2d(iphi.flatten())
        itau = np.atleast_2d(itau.flatten())

        A1 = (np.sin(iphi * dphi) * cosw_moments[np.abs(itau)]
            - np.cos(iphi * dphi) * np.sign(itau) * sinw_moments[np.abs(itau)])
        A2 = (-np.sin(iphi * dphi) * wcosw_moments[np.abs(itau)]
            + np.cos(iphi * dphi) * np.sign(itau) * wsinw_moments[np.abs(itau)])
        A = np.asmatrix(np.vstack((A1, A2)))

        # Create test points at all possible combinations of test values.
        Phi = np.asmatrix(np.vstack((iphi * dphi, itau * dtau)))

        iphi, kphi = np.meshgrid(iphi.flatten(), iphi.flatten())
        itau, ktau = np.meshgrid(itau.flatten(), itau.flatten())

        B = np.asmatrix(np.exp(snr2 * (1
            + np.cos((iphi - kphi) * dphi) * cosw_moments[np.abs(itau - ktau)]
            + np.sin((iphi - kphi) * dphi) * np.sign(itau - ktau) * sinw_moments[np.abs(itau - ktau)]
            - np.cos(iphi * dphi) * cosw_moments[np.abs(itau)]
            - np.sin(iphi * dphi) * np.sign(itau) * sinw_moments[np.abs(itau)]
            - np.cos(kphi * dphi) * cosw_moments[np.abs(ktau)]
            - np.sin(kphi * dphi) * np.sign(ktau) * sinw_moments[np.abs(ktau)])))

        LambdaInvA = np.asmatrix(linalg.solve(Lambda, A, sym_pos=True))
        Y = Phi - LambdaInvA
        Delta = B - snr2 * A.T * LambdaInvA
        x = np.asmatrix(linalg.solve(Delta, Y.T))
        BRB = Lambda.I / snr2 + Y * x

        # FIXME: stray factor of pi from somewhere?
        BRB /= np.square(np.pi)

        return BRB

    def get_cov(self, snr):
        """Get the Barankin bound if snr < 20, or the Cramer-Rao bound otherwise."""
        if snr < 20:
            return self.get_brb(snr)
        else:
            return self.get_crb(snr)

    # FIXME: np.vectorize doesn't work on unbound instance methods. The excluded
    # keyword, added in Numpy 1.7, could be used here to exclude the zeroth
    # argument, self.
    def __get_brb_toa_uncert(self, snr):
        return np.sqrt(self.get_brb(snr)[1, 1])
    def get_brb_toa_uncert(self, snr):
        return np.vectorize(self.__get_brb_toa_uncert)(snr)

    # FIXME: np.vectorize doesn't work on unbound instance methods. The excluded
    # keyword, added in Numpy 1.7, could be used here to exclude the zeroth
    # argument, self.
    def __get_crb_toa_uncert(self, snr):
        return np.sqrt(self.get_crb(snr)[1, 1])
    def get_crb_toa_uncert(self, snr):
        return np.vectorize(self.__get_crb_toa_uncert)(snr)

    def get_toa_uncert(self, snr):
        """Get timing uncertainty for a given signal to noise ratio. Use the
        Barankin bound for snr < 20, and the Cramer-Rao bound for snr >= 20."""
        return np.piecewise(snr, [snr < 20],
            [self.get_brb_toa_uncert, self.get_crb_toa_uncert])


if __name__ == '__main__':
    # Demo; execute with "python -m bayestar_localization.timing"
    from matplotlib import pyplot as plt
    mass_pairs = ((1.4, 1.4), (1.4, 10.), (10., 10.))
    colors = 'rgb'
    plt.figure(figsize=(6, 6))
    ax = plt.subplot(111, aspect=1)
    for color, mass_pair in zip(colors, mass_pairs):
        signal_model = SignalModel(*mass_pair, S=get_noise_psd_func("H1"), f_low=30)
        snr = np.logspace(0, 2, 20)
        plt.loglog(snr, signal_model.get_toa_uncert(snr), '-o', mew=2, lw=2, mfc='none', mec=color, color=color)
        plt.loglog(snr, signal_model.get_crb_toa_uncert(snr), '--', color=color)
        plt.annotate(ur"%g—%g $M_\odot$" % mass_pair, (15, signal_model.get_crb_toa_uncert(15)), rotation=-45)
    plt.xlabel("S/N")
    plt.ylabel("TOA uncertainty (s)")
    plt.grid(which='minor')
    plt.grid(which='major', lw=1)
    plt.ylim(1e-5, 1e-3)
    plt.savefig('toa_uncert.pdf')
