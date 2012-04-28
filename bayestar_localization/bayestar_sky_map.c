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

#include "bayestar_sky_map.h"

#include <complex.h>
#include <float.h>
#include <math.h>
#include <string.h>

#include <lal/DetResponse.h>
#include <lal/LALSimulation.h>

#include <chealpix.h>

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sort_vector_double.h>
#include <gsl/gsl_vector.h>


/* Copied from lal's TimeDelay.c */
/* scalar product of two 3-vectors */
static double dotprod(const double vec1[3], const double vec2[3])
{
	return vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2];
}

/* Copied from lal's TimeDelay.c, but takes gmst as an argument instead of GPS time. */
static double
_XLALArrivalTimeDiff(
	const double detector1_earthfixed_xyz_metres[3],
	const double detector2_earthfixed_xyz_metres[3],
	const double source_right_ascension_radians,
	const double source_declination_radians,
	const double gmst
)
{
	double delta_xyz[3];
	double ehat_src[3];
	const double greenwich_hour_angle = gmst - source_right_ascension_radians;

	if(XLAL_IS_REAL8_FAIL_NAN(greenwich_hour_angle))
		XLAL_ERROR_REAL8(XLAL_EFUNC);

	/*
	 * compute the unit vector pointing from the geocenter to the
	 * source
	 */

	ehat_src[0] = cos(source_declination_radians) * cos(greenwich_hour_angle);
	ehat_src[1] = cos(source_declination_radians) * -sin(greenwich_hour_angle);
	ehat_src[2] = sin(source_declination_radians);

	/*
	 * position of detector 2 with respect to detector 1
	 */

	delta_xyz[0] = detector2_earthfixed_xyz_metres[0] - detector1_earthfixed_xyz_metres[0];
	delta_xyz[1] = detector2_earthfixed_xyz_metres[1] - detector1_earthfixed_xyz_metres[1];
	delta_xyz[2] = detector2_earthfixed_xyz_metres[2] - detector1_earthfixed_xyz_metres[2];

	/*
	 * Arrival time at detector 1 - arrival time at detector 2.  This
	 * is positive when the wavefront arrives at detector 1 after
	 * detector 2 (and so t at detector 1 is greater than t at detector
	 * 2).
	 */

	return dotprod(ehat_src, delta_xyz) / LAL_C_SI;
}


static double square(double a)
{
    return a * a;
}


/* Perform sky localization based on TDOAs alone. Not normalized. */
static int bayestar_sky_map_tdoa_not_normalized(
    long npix, /* Input: number of HEALPix pixels. */
    double *restrict P, /* Output: pre-allocated array of length npix to store posterior map. */
    double gmst, /* Greenwich mean sidereal time in radians. */
    int nifos, /* Input: number of detectors. */
    const double **restrict locs, /* Input: array of detector positions. */
    const double *restrict toas, /* Input: array of times of arrival. */
    const double *restrict s2_toas /* Input: uncertainties in times of arrival. */
) {
    double tdoas[nifos - 1], w_tdoas[nifos - 1];
    long nside;
    int ret = -1;
    long i;

    /* Check that none of the inputs are NULL. */
    if (!P || !locs || !toas || !s2_toas)
    {
        errno = EFAULT;
        goto fail;
    }

    /* Determine the lateral HEALPix resolution. */
    nside = npix2nside(npix);
    if (nside < 0)
        goto fail;

    /* Compute time delays on arrival (TDOAS) relative to detector 0,
     * and weights for their measurement uncertainties. */
    for (i = 1; i < nifos; i ++)
    {
        tdoas[i - 1] = toas[i] - toas[0];
        w_tdoas[i - 1] = 1 / (s2_toas[i] + s2_toas[0]);
    }

    /* Loop over pixels. */
    for (i = 0; i < npix; i ++)
    {
        long j;
        double theta, phi, accum;

        /* Determine spherical polar coordinates of this pixel. */
        pix2ang_ring(nside, i, &theta, &phi);

        /* Loop over detector pairs. */
        for (accum = 0, j = 1; j < nifos; j ++)
        {
            /* Compute expected TDOA for this sky location and detector pair. */
            const double expected_tdoa = _XLALArrivalTimeDiff(locs[j], locs[0],
                phi, M_PI_2 - theta, gmst);

            /* Accumulate the measurement error. */
            accum += square(tdoas[j - 1] - expected_tdoa) * w_tdoas[j - 1];
        }

        /* Evaluate the (un-normalized) Gaussian PDF for the total measurement error. */
        P[i] = exp(-0.5 * accum);
    }

    /* Done! */
    ret = 0;
fail:
    return ret;
}


/* Perform sky localization based on TDOAs alone. */
int bayestar_sky_map_tdoa(
    long npix, /* Input: number of HEALPix pixels. */
    double *restrict P, /* Output: pre-allocated array of length npix to store posterior map. */
    double gmst, /* Greenwich mean sidereal time in radians. */
    int nifos, /* Input: number of detectors. */
    const double **restrict locs, /* Input: array of detector positions. */
    const double *restrict toas, /* Input: array of times of arrival. */
    const double *restrict s2_toas /* Input: uncertainties in times of arrival. */
) {
    int ret = bayestar_sky_map_tdoa_not_normalized(npix, P, gmst, nifos, locs, toas, s2_toas);
    if (ret == 0)
    {
        /* Sort pixel indices by ascending significance. */
        double accum;
        long i;
        gsl_vector_view P_vector = gsl_vector_view_array(P, npix);
        gsl_permutation *pix_perm = gsl_permutation_alloc(npix);
        if (!pix_perm)
            return -1;
        gsl_sort_vector_index(pix_perm, &P_vector.vector);

        for (accum = 0, i = 0; i < npix; i ++)
            accum += P[gsl_permutation_get(pix_perm, npix - i - 1)];
        for (i = 0; i < npix; i ++)
            P[i] /= accum;
        gsl_permutation_free(pix_perm);
    }
    return ret;
}


#define INTEGRAND_COUNT_NODES 15
#define INTEGRAND_COUNT_SAMPLES (INTEGRAND_COUNT_NODES * INTEGRAND_COUNT_NODES)


typedef struct {
    double args[INTEGRAND_COUNT_SAMPLES][2];
} inner_integrand_params;


static double inner_integrand(double x, void *params)
{
    const inner_integrand_params *my_params = (inner_integrand_params *) params;

    /* We'll need to divide by 1/x and also 1/x^2. We can save one division by
     * computing 1/x and then squaring it. */
    const double onebyx = 1 / x;
    const double onebyx2 = onebyx * onebyx;

    double ret;
    int i;

    for (ret = 0, i = 0; i < INTEGRAND_COUNT_SAMPLES; i ++)
    {
        const double I0_arg = onebyx * my_params->args[i][1];
        ret += exp(I0_arg + onebyx2 * my_params->args[i][0])
            * gsl_sf_bessel_I0_scaled(I0_arg);
    }

    return ret * x * x;
}


int bayestar_sky_map(
    long npix, /* Input: number of HEALPix pixels. */
    double *restrict P, /* Output: pre-allocated array of length npix to store posterior map. */
    double gmst, /* Greenwich mean sidereal time in radians. */
    int nifos, /* Input: number of detectors. */
    const float **restrict responses, /* Pointers to detector responses. */
    const double **restrict locations, /* Pointers to locations of detectors in Cartesian geographic coordinates. */
    const double *restrict toas, /* Input: array of times of arrival with arbitrary relative offset. (Make toas[0] == 0.) */
    const double complex *restrict snrs, /* Input: array of SNRs. */
    const double *restrict s2_toas, /* Measurement variance of TOAs. */
    const double *restrict horizons /* Distances at which a source would produce an SNR of 1 in each detector. */)
{
    long nside;
    long maxpix;
    int ret = -1;
    long i;
    double d1[nifos];
    gsl_permutation *pix_perm = NULL;
    gsl_integration_glfixed_table *glfixed_table = NULL;

    /* Precalculate trigonometric that occur in the integrand. */
    double u4_6u2_1[INTEGRAND_COUNT_NODES];
    double u4_2u2_1[INTEGRAND_COUNT_NODES];
    double u3_u[INTEGRAND_COUNT_NODES];
    double cosines[INTEGRAND_COUNT_NODES];
    double sines[INTEGRAND_COUNT_NODES];
    for (i = -(INTEGRAND_COUNT_NODES / 2); i <= (INTEGRAND_COUNT_NODES / 2); i ++)
    {
        const double u = (double) i / (INTEGRAND_COUNT_NODES / 2);
        const double u2 = u * u;
        const double u3 = u2 * u;
        const double u4 = u3 * u;
        const double angle = M_PI * u;
        u4_6u2_1[i + (INTEGRAND_COUNT_NODES / 2)] = u4 + 6 * u2 + 1;
        u4_2u2_1[i + (INTEGRAND_COUNT_NODES / 2)] = u4 - 2 * u2 + 1;
        u3_u[i + (INTEGRAND_COUNT_NODES / 2)] = u3 + u;
        cosines[i + (INTEGRAND_COUNT_NODES / 2)] = cos(angle);
        sines[i + (INTEGRAND_COUNT_NODES / 2)] = sin(angle);
    }

    /* Determine the lateral HEALPix resolution. */
    nside = npix2nside(npix);
    if (nside < 0)
        goto fail;

    /* Check that none of the inputs are NULL. */
    if (!P|| !responses || !locations || !toas || !snrs || !horizons || !s2_toas)
    {
        errno = EFAULT;
        goto fail;
    }

    /* Allocate temporary spaces. */
    pix_perm = gsl_permutation_alloc(npix);
    if (!pix_perm)
        goto fail;

    /* Precalculate weights and nodes for Gaussian quadrature. */
    glfixed_table = gsl_integration_glfixed_table_alloc(32);
    if (!glfixed_table)
        goto fail;

    /* Rescale distances so that furthest horizon distance is 1. */
    {
        double d1max;
        memcpy(d1, horizons, sizeof(d1));
        for (d1max = d1[0], i = 1; i < nifos; i ++)
            if (d1[i] > d1max)
                d1max = d1[i];
        for (i = 0; i < nifos; i ++)
            d1[i] /= d1max;
    }

    /* Evaluate posterior term only first. */
    if (bayestar_sky_map_tdoa_not_normalized(npix, P, gmst, nifos, locations, toas, s2_toas))
        goto fail;

    /* Sort pixel indices by ascending significance. */
    {
        gsl_vector_view P_vector = gsl_vector_view_array(P, npix);
        gsl_sort_vector_index(pix_perm, &P_vector.vector);
    }

    /* Find the number of pixels needed to account for 99.99% of the posterior
     * conditioned on TDOAs. */
    {
        double accum, Ptotal;
        for (Ptotal = 0, i = 0; i < npix; i ++)
            Ptotal += P[gsl_permutation_get(pix_perm, npix - i - 1)];
        for (accum = 0, maxpix = 0; maxpix < npix && accum <= 0.9999 * Ptotal; maxpix ++)
            accum += P[gsl_permutation_get(pix_perm, npix - maxpix - 1)];
    }

    /* Compute posterior factor for amplitude consistency. */
    #pragma omp parallel for
    for (i = 0; i < maxpix; i ++)
    {
        int j;
        long ipix = gsl_permutation_get(pix_perm, npix - i - 1);
        inner_integrand_params params;

        /* Evaluate integrand on a regular lattice in psi and cos(i). */
        {
            double e, f, g;
            double c1, c2, c3, c4;
            {
                double theta, phi;
                double a, b, c, d;
                pix2ang_ring(nside, ipix, &theta, &phi);
                for (j = 0, a = 0, b = 0, c = 0, d = 0, e = 0, f = 0, g = 0; j < nifos; j ++)
                {
                    double Fp, Fx;
                    XLALComputeDetAMResponse(&Fp, &Fx, responses[j], phi, M_PI_2 - theta, 0, gmst);
                    Fp *= d1[j];
                    Fx *= d1[j];
                    a += Fp * creal(snrs[j]);
                    b += Fx * cimag(snrs[j]);
                    c += Fx * creal(snrs[j]);
                    d += Fp * cimag(snrs[j]);
                    g += Fp * Fx;
                    Fp *= Fp;
                    Fx *= Fx;
                    e += Fp + Fx;
                    f += Fp - Fx;
                }
                c1 = a * b - c * d;
                c4 = (a * c + b * d) / 4;
                a *= a;
                b *= b;
                c *= c;
                d *= d;
                c2 = (a + b + c + d) / 8;
                c3 = (a - b - c + d) / 8;
            }
            e /= -16;
            f /= -16;
            g /= -16;

            for (j = 0; j < INTEGRAND_COUNT_NODES; j ++)
            {
                const double args0 = u4_6u2_1[j] * e;
                const double args1 = u3_u[j] * c1 + u4_6u2_1[j] * c2;
                int k;

                for (k = 0; k < INTEGRAND_COUNT_NODES; k ++)
                {
                    int l = j * INTEGRAND_COUNT_NODES + k;
                    params.args[l][0] = args0 + u4_2u2_1[j] * (f * cosines[k] + g * sines[k]);
                    params.args[l][1] = sqrt(args1 + u4_2u2_1[j] * (c3 * cosines[k] + c4 * sines[k]));
                }
            }
        }

        /* Perform Gaussian quadrature over luminosity distance. */
        {
            const gsl_function func = {inner_integrand, &params};
            P[ipix] *= gsl_integration_glfixed(&func, 0.05, 4, glfixed_table);
        }
    }

    /* Finish up by zeroing pixels that didn't meet the TDOA cut. */
    for (i = maxpix; i < npix; i ++)
    {
        long ipix = gsl_permutation_get(pix_perm, npix - i - 1);
        P[ipix] = 0;
    }

    /* Normalize posterior. */
    {
        double accum;
        for (accum = 0, i = 0; i < npix; i ++)
            accum += P[gsl_permutation_get(pix_perm, npix - i - 1)];
        for (i = 0; i < npix; i ++)
            P[i] /= accum;
    }

    ret = 0;
fail:
    gsl_permutation_free(pix_perm);
    gsl_integration_glfixed_table_free(glfixed_table);
    return ret;
}
