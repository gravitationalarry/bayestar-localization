from __future__ import division
"""
"""
__author__ = "Leo Singer <leo.singer@ligo.org>"


from itertools import izip
import numpy as np
from bayestar_localization._sky_map import *


def create_plan(nodes):
    if nodes % 2 != 1:
        raise ValueError("Nodes should be odd")
    plan = np.empty((nodes ** 3, 2, 2))
    angles = 2 * np.pi * np.arange(nodes) / nodes
    cosines = np.cos(angles)
    sines = np.sin(angles)
    ifac2s = np.arange(-(nodes // 2), nodes // 2 + 1) / (nodes // 2)
    ifac1s = 0.5 * (1 + ifac2s * ifac2s)

    for i, (cosphi, sinphi) in enumerate(izip(cosines, sines)):
        mA = np.asmatrix(((cosphi, sinphi), (-sinphi, cosphi)))
        for j, (ifac1, ifac2) in enumerate(izip(ifac1s, ifac2s)):
            mB = mA * np.diag((ifac1, ifac2))
            for k, (cospsi, sinpsi) in enumerate(izip(cosines, sines)):
                mC = mB * ((cospsi, sinpsi), (-sinpsi, cospsi))
                plan[i * nodes ** 2 + j * nodes + k, :, :] = mC

    return plan
