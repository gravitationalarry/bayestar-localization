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

#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
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
    double *P, /* Output: pre-allocated array of length npix to store posterior map. */
    double gmst, /* Greenwich mean sidereal time in radians. */
    int nifos, /* Input: number of detectors. */
    const double **locs, /* Input: array of detector positions. */
    const double *toas, /* Input: array of times of arrival. */
    const double *s2_toas /* Input: uncertainties in times of arrival. */
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
    double *P, /* Output: pre-allocated array of length npix to store posterior map. */
    double gmst, /* Greenwich mean sidereal time in radians. */
    int nifos, /* Input: number of detectors. */
    const double **locs, /* Input: array of detector positions. */
    const double *toas, /* Input: array of times of arrival. */
    const double *s2_toas /* Input: uncertainties in times of arrival. */
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


#define INTEGRAND_COUNT_NODES 31
#define INTEGRAND_COUNT_SAMPLES (INTEGRAND_COUNT_NODES * INTEGRAND_COUNT_NODES)


typedef struct {
    double a;
    double log_offset;
} inner_integrand_params;


/* Inner (radial) integrand. */
static double inner_integrand(double log_x, void *params)
{
    const inner_integrand_params *integrand_params = (const inner_integrand_params *) params;

    /* We'll need to divide by 1/x and also 1/x^2. We can save one division by
     * computing 1/x and then squaring it. */
    const double onebyx = exp(-log_x);
    const double onebyx2 = onebyx * onebyx;
    const double I0_arg = onebyx / integrand_params->a;
    return exp(I0_arg - 0.5 * onebyx2 - integrand_params->log_offset) * gsl_sf_bessel_I0_scaled(I0_arg);
}


int bayestar_sky_map_tdoa_snr(
    long npix, /* Input: number of HEALPix pixels. */
    double *P, /* Output: pre-allocated array of length npix to store posterior map. */
    double gmst, /* Greenwich mean sidereal time in radians. */
    int nifos, /* Input: number of detectors. */
    const float **responses, /* Pointers to detector responses. */
    const double **locations, /* Pointers to locations of detectors in Cartesian geographic coordinates. */
    const double *toas, /* Input: array of times of arrival with arbitrary relative offset. (Make toas[0] == 0.) */
    const double complex *snrs, /* Input: array of SNRs. */
    const double *s2_toas, /* Measurement variance of TOAs. */
    const double *horizons, /* Distances at which a source would produce an SNR of 1 in each detector. */
    double min_distance,
    double max_distance)
{
    long nside;
    long maxpix;
    int ret = -1;
    long i;
    double d1[nifos];
    gsl_permutation *pix_perm = NULL;

    /* Maximum number of subdivisions for adaptive integration. */
    static const size_t subdivision_limit = 64;

    /* Constants that determine how much of the peak in the likelihood to enclose in integration break points. */
    static const double y1 = 0.01, y2 = 0.005;
    const double upper_breakpoint_default = (log(y2) - sqrt(log(y1) * log(y2))) / (sqrt(-2 * log(y2)) * (log(y2) - log(y1)));

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

    /* Rescale distances so that furthest horizon distance is 1. */
    {
        double d1max;
        memcpy(d1, horizons, sizeof(d1));
        for (d1max = d1[0], i = 1; i < nifos; i ++)
            if (d1[i] > d1max)
                d1max = d1[i];
        for (i = 0; i < nifos; i ++)
            d1[i] /= d1max;
        min_distance /= d1max;
        max_distance /= d1max;
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
        long ipix = gsl_permutation_get(pix_perm, npix - i - 1);

        /* Pre-compute some coefficients in the integrand that are determined by antenna factors and SNR. */
        double e, f, g;
        double c1, c2, c3, c4;
        {
            int j;
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
        e /= 8;
        f /= 8;
        g /= 8;

        /* Evaluate integral on a regular lattice in psi and cos(i) and using
         * adaptive quadrature over log(distance). */
        {
            /* Prepare workspace for adaptive integrator. */
            gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(subdivision_limit);

            int j;
            double accum;

            /* Loop over cos(i). */
            for (accum = -INFINITY, j = 0; j < INTEGRAND_COUNT_NODES; j ++)
            {
                /* Coefficients in integrand that are determined by antenna factors, SNR, and cos(i), but not psi. */
                const double args0 = u4_6u2_1[j] * e;
                const double args1 = u3_u[j] * c1 + u4_6u2_1[j] * c2;

                int k;

                /* Loop over psi. */
                for (k = 0; k < INTEGRAND_COUNT_NODES; k ++)
                {
                    /* Variables to store output from integrator. */
                    double result, abserr;

                    /* Coefficients in integrand that are determined by antenna factors, SNR, cos(i), and psi. */
                    const double num = args1 + u4_2u2_1[j] * (c3 * cosines[k] + c4 * sines[k]);
                    const double den = args0 + u4_2u2_1[j] * (f * cosines[k] + g * sines[k]);
                    const double a2 = den / num;
                    const double a = sqrt(a2);
                    const double sqrt_den = sqrt(den);

                    /* Create data structures for integrand callback. */
                    const inner_integrand_params integrand_params = {a, 0.5 / a2};
                    const gsl_function func = {inner_integrand, &integrand_params};

                    /* Limits of integration. */
                    const double x1 = min_distance / sqrt_den;
                    const double x2 = max_distance / sqrt_den;

                    /* Find break points in integration interval that enclose the peak in the likelihood function. */
                    const double lower_breakpoint = (a - a2 * sqrt(-2 * log(y1))) / (1 + 2 * a2 * log(y1));
                    const double upper_breakpoint = (a < 1 / sqrt(-2 * log(y2)))
                        ? ((a + a2 * sqrt(-2 * log(y1))) / (1 + 2 * a2 * log(y1)))
                        : upper_breakpoint_default;

                    /* Create list of integration break points. */
                    double breakpoints[4];
                    int num_breakpoints = 0;
                    /* Always start with lower limit of integration. */
                    breakpoints[num_breakpoints++] = log(x1);
                    /* If integration interval contains lower break point, add it. */
                    if (lower_breakpoint > x1 && lower_breakpoint < x2)
                        breakpoints[num_breakpoints++] = log(lower_breakpoint);
                    /* If integration interval contains upper break point, add it too. */
                    if (upper_breakpoint > x1 && upper_breakpoint < x2)
                        breakpoints[num_breakpoints++] = log(upper_breakpoint);
                    /* Always end with upper limit of integration. */
                    breakpoints[num_breakpoints++] = log(x2);

                    /* Perform adaptive integration. Stop when a relative accuracy of 0.01 has been reached. */
                    gsl_integration_qagp(&func, &breakpoints[0], num_breakpoints, DBL_MIN, 0.01, subdivision_limit, workspace, &result, &abserr);

                    /* Accumulate the (log) posterior for this cos(i) and psi. */
                    {
                        double max_log_p;
                        result = log(result) + integrand_params.log_offset;
                        max_log_p = fmax(result, accum);
                        accum = log(exp(result - max_log_p) + exp(accum - max_log_p)) + max_log_p;
                    }
                }
            }

            /* Discard workspace for adaptive integrator. */
            gsl_integration_workspace_free(workspace);

            /* Accumulate (log) posterior terms for SNR and TDOA. */
            P[ipix] = log(P[ipix]) + log(accum);
        }
    }

    /* Finish up by zeroing pixels that didn't meet the TDOA cut. */
    for (i = maxpix; i < npix; i ++)
    {
        long ipix = gsl_permutation_get(pix_perm, npix - i - 1);
        P[ipix] = -INFINITY;
    }

    /* Normalize posterior. */
    {
        /* Find maximum of log posterior. */
        double accum;
        for (accum = P[0], i = 1; i < npix; i ++)
        {
            double new_log_p = P[i];
            if (new_log_p > accum)
                accum = new_log_p;
        }

        /* Subtract off maximum of log posterior, and then exponentiate to get the posterior itself. */
        for (i = 0; i < npix; i ++)
            P[i] = exp(P[i] - accum);

        /* Sum posterior to get normalization. */
        for (accum = 0, i = 0; i < npix; i ++)
            accum += P[gsl_permutation_get(pix_perm, npix - i - 1)];

        /* Normalize posterior. */
        for (i = 0; i < npix; i ++)
            P[i] /= accum;
    }

    /* Done !*/
    ret = 0;
fail:
    gsl_permutation_free(pix_perm);
    return ret;
}
