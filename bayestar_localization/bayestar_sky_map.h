/*
 * Copyright (C) 2012  Leo Singer
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <lal/LALDatatypes.h>

#ifndef BAYESTAR_SKY_MAP_H
#define BAYESTAR_SKY_MAP_H

/* Perform sky localization based on TDOAs alone. */
int bayestar_sky_map_tdoa(
    long npix, /* Input: number of HEALPix pixels. */
    double *restrict P, /* Output: pre-allocated array of length npix to store posterior map. */
    double gmst, /* Greenwich mean sidereal time in radians. */
    int nifos, /* Input: number of detectors. */
    const double **restrict locs, /* Input: array of detector positions. */
    const double *restrict toas, /* Input: array of times of arrival. */
    const double *restrict s2_toas /* Input: uncertainties in times of arrival. */
);

/* Perform sky localization based on TDOAs and amplitude. */
int bayestar_sky_map(
    long npix, /* Input: number of HEALPix pixels. */
    double *restrict P, /* Output: pre-allocated array of length npix to store posterior map. */
    long nplan, /* Number of entries in plan; plan should be an array of 4 * nplan doubles. */
    const double *restrict plan, /* Pre-evaluated transformation matrices. */
    double gmst, /* Greenwich mean sidereal time in radians. */
    int nifos, /* Input: number of detectors. */
    const float **restrict responses, /* Pointers to detector responses. */
    const double **restrict locations, /* Pointers to locations of detectors in Cartesian geographic coordinates. */
    const double *restrict toas, /* Input: array of times of arrival with arbitrary relative offset. (Make toas[0] == 0.) */
    const double complex *restrict snrs, /* Input: array of SNRs. */
    const double *restrict s2_toas, /* Measurement variance of TOAs. */
    const double *restrict horizons /* Distances at which a source would produce an SNR of 1 in each detector. */
);

#endif /* BAYESTAR_SKY_MAP_H */
