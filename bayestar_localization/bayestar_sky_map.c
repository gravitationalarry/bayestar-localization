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

#include <gsl/gsl_errno.h>
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


static double cabs2(const double complex x)
{
    return square(creal(x)) + square(cimag(x));
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


/* Custom error handler. */
static void my_gsl_error(const char * reason, const char * file, int line, int gsl_errno)
{
    /* For 'maximum number of subdivisions reached' errors, print a warning and move on. */
    if (errno == GSL_EMAXITER)
        gsl_stream_printf("WARNING", file, line, reason);
    else
    /* For all other errors, use the default handler: print an error message and abort. */
        gsl_error(reason, file, line, errno);
}


struct inner_integrand_params {
    int nifos;
    double complex *response;
    const double complex *snrs;
};


static double inner_integrand(double x, void *params)
{
    int i;
    double ret;
    int nifos = ((struct inner_integrand_params *) params)->nifos;
    const double complex *response = ((struct inner_integrand_params *) params)->response;
    const double complex *snrs = ((struct inner_integrand_params *) params)->snrs;

    for (ret = 0, i = 0; i < nifos; i ++)
        ret += cabs2(response[i] / x - snrs[i]);
    ret = exp(-0.5 * ret) * x * x;
    return ret;
}


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
    const double *restrict horizons /* Distances at which a source would produce an SNR of 1 in each detector. */)
{
    long nside;
    long maxpix;
    int ret = -1;
    long i;
    double d1[nifos];
    gsl_permutation *pix_perm;
    gsl_error_handler_t *old_handler = gsl_set_error_handler(my_gsl_error);

    /* Determine the lateral HEALPix resolution. */
    nside = npix2nside(npix);
    if (nside < 0)
        goto fail;

    /* Check that none of the inputs are NULL. */
    if (!P|| !plan || !responses || !locations || !toas || !snrs || !horizons || !s2_toas)
    {
        errno = EFAULT;
        goto fail;
    }

    /* Allocate temporary spaces. */
    pix_perm = gsl_permutation_alloc(npix);
    if (!pix_perm)
        goto fail;

    /* Rescale distances so that furthest horizon distance is 1 */
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
        long j;
        double theta, phi, accum;
        double Fp[nifos], Fx[nifos];
        long ipix = gsl_permutation_get(pix_perm, npix - i - 1);
        gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(20);
        pix2ang_ring(nside, ipix, &theta, &phi);

        for (j = 0; j < nifos; j ++)
        {
            XLALComputeDetAMResponse(&Fp[j], &Fx[j], responses[j],
                phi, M_PI_2 - theta, 0, gmst);
            Fp[j] *= d1[j];
            Fx[j] *= d1[j];
        }
        for (accum = 0, j = 0; j < nplan; j ++)
        {
            long k;
            const double *M = &plan[4 * j];
            double complex response[nifos];
            double result, abserr;
            struct inner_integrand_params params = {nifos, response, snrs};
            gsl_function func = {inner_integrand, &params};
            for (k = 0; k < nifos; k ++)
                response[k] = (M[0] * Fp[k] + M[1] * Fx[k]) + (M[2] * Fp[k] + M[3] * Fx[k]) * 1i;
            gsl_integration_qag(&func, 0, 4, DBL_MIN, 0.05, 20, GSL_INTEG_GAUSS15, workspace, &result, &abserr);
            accum += result;
        }
        gsl_integration_workspace_free(workspace);
        P[ipix] *= accum;
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
    gsl_set_error_handler(old_handler);
    return ret;
}
