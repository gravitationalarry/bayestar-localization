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

#include <Python.h>
#include <numpy/arrayobject.h>
#include <chealpix.h>
#include <gsl/gsl_errno.h>
#include "bayestar_sky_map.h"


static void
my_gsl_error (const char * reason, const char * file, int line, int gsl_errno)
{
    PyObject *exception_type;
    switch (gsl_errno)
    {
        case GSL_EINVAL:
            exception_type = PyExc_ValueError;
            break;
        case GSL_ENOMEM:
            exception_type = PyExc_MemoryError;
            break;
        default:
            exception_type = PyExc_RuntimeError;
            break;
    }
    PyErr_Format(exception_type, "%s:%d: %s\n", file, line, reason);
}


static PyObject *AsArrayObj(PyObject *obj, int typenum, int depth)
{
    return PyArray_CheckFromAny(obj,
        PyArray_DescrFromType(typenum), depth, depth,
        NPY_NOTSWAPPED | NPY_ELEMENTSTRIDES | NPY_CONTIGUOUS | NPY_ALIGNED,
        NULL);
}


static PyObject *sky_map_tdoa(PyObject *module, PyObject *args, PyObject *kwargs)
{
    long i;
    Py_ssize_t n;
    long nside;
    long npix;
    long nifos;
    double gmst;
    PyObject *toas_obj, *toa_variances_obj, *locations_obj;

    PyObject *toas_npy = NULL, *toa_variances_npy = NULL, **locations_npy = NULL;

    double *toas;
    double *toa_variances;
    double **locations = NULL;

    npy_intp dims[1];
    PyObject *out = NULL, *ret = NULL;
    double *P;
    int result;
    gsl_error_handler_t *old_handler;

    /* Names of arguments */
    static char *keywords[] = {"nside", "gmst", "toas",
        "toa_variances", "locations", NULL};

    /* Parse arguments */
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "ldOOO|", keywords,
        &nside, &gmst, &toas_obj, &toa_variances_obj, &locations_obj))
        goto fail;

    npix = nside2npix(nside);
    if (npix < 0)
    {
        PyErr_SetString(PyExc_ValueError, "nside must be a power of 2");
        goto fail;
    }
    dims[0] = npix;
    out = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (!out)
        goto fail;
    P = PyArray_DATA(out);

    toas_npy = AsArrayObj(toas_obj, NPY_DOUBLE, 1);
    if (!toas_npy) goto fail;
    nifos = PyArray_DIM(toas_npy, 0);
    toas = PyArray_DATA(toas_npy);

    toa_variances_npy = AsArrayObj(toa_variances_obj, NPY_DOUBLE, 1);
    if (!toa_variances_npy) goto fail;
    if (PyArray_DIM(toa_variances_npy, 0) != nifos)
    {
        PyErr_SetString(PyExc_ValueError, "toas and toa_variances must have the same length");
        goto fail;
    }
    toa_variances = PyArray_DATA(toa_variances_npy);

    locations_npy = malloc(nifos * sizeof(PyObject *));
    if (!locations_npy)
    {
        PyErr_SetNone(PyExc_MemoryError);
        goto fail;
    }
    for (i = 0; i < nifos; i ++)
        locations_npy[i] = NULL;
    locations = malloc(nifos * sizeof(double *));
    if (!locations)
    {
        PyErr_SetNone(PyExc_MemoryError);
        goto fail;
    }

    n = PySequence_Length(locations_obj);
    if (n < 0) goto fail;
    if (n != nifos)
    {
        PyErr_SetString(PyExc_ValueError, "toas and locations must have the same length");
        goto fail;
    }
    for (i = 0; i < nifos; i ++)
    {
        PyObject *obj = PySequence_GetItem(locations_obj, i);
        if (!obj) goto fail;
        locations_npy[i] = AsArrayObj(obj, NPY_DOUBLE, 1);
        Py_XDECREF(obj);
        if (!locations_npy[i]) goto fail;
        if (PyArray_DIM(locations_npy[i], 0) != 3)
        {
            PyErr_SetString(PyExc_ValueError, "expected every element of locations to be a vector of length 3");
            goto fail;
        }
        locations[i] = PyArray_DATA(locations_npy[i]);
    }

    old_handler = gsl_set_error_handler(my_gsl_error);
    result = bayestar_sky_map_tdoa(npix, P, gmst, nifos, locations, toas, toa_variances);
    gsl_set_error_handler(old_handler);
    if (result != GSL_SUCCESS)
        goto fail;

    ret = out;
    out = NULL;
fail:
    Py_XDECREF(toas_npy);
    Py_XDECREF(toa_variances_npy);
    if (locations_npy)
        for (i = 0; i < nifos; i ++)
            Py_XDECREF(locations_npy[i]);
    free(locations_npy);
    free(locations);
    Py_XDECREF(out);
    return ret;
};


static PyObject *sky_map_tdoa_snr(PyObject *module, PyObject *args, PyObject *kwargs)
{
    long i;
    Py_ssize_t n;
    long nside;
    long npix;
    long nifos;
    double gmst;
    PyObject *toas_obj, *snrs_obj, *toa_variances_obj, *responses_obj,
        *locations_obj, *horizons_obj;

    PyObject *toas_npy = NULL, *snrs_npy = NULL, *toa_variances_npy = NULL, **responses_npy = NULL, **locations_npy = NULL, *horizons_npy = NULL;
    char *prior_str = NULL;

    double *toas;
    double complex *snrs;
    double *toa_variances;
    float **responses = NULL;
    double **locations = NULL;
    double *horizons;

    double min_distance, max_distance;
    bayestar_prior_t prior = -1;

    npy_intp dims[1];
    PyObject *out = NULL, *ret = NULL;
    double *P;
    int result;
    gsl_error_handler_t *old_handler;

    /* Names of arguments */
    static char *keywords[] = {"nside", "gmst", "toas", "snrs",
        "toa_variances", "responses", "locations", "horizons",
        "min_distance", "max_distance", "prior", NULL};

    /* Parse arguments */
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "ldOOOOOOdds|", keywords,
        &nside, &gmst, &toas_obj, &snrs_obj, &toa_variances_obj,
        &responses_obj, &locations_obj, &horizons_obj,
        &min_distance, &max_distance, &prior_str)) goto fail;

    npix = nside2npix(nside);
    if (npix < 0)
    {
        PyErr_SetString(PyExc_ValueError, "nside must be a power of 2");
        goto fail;
    }
    dims[0] = npix;
    out = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (!out)
        goto fail;
    P = PyArray_DATA(out);

    toas_npy = AsArrayObj(toas_obj, NPY_DOUBLE, 1);
    if (!toas_npy) goto fail;
    nifos = PyArray_DIM(toas_npy, 0);
    toas = PyArray_DATA(toas_npy);

    snrs_npy = AsArrayObj(snrs_obj, NPY_CDOUBLE, 1);
    if (!snrs_npy) goto fail;
    if (PyArray_DIM(snrs_npy, 0) != nifos)
    {
        PyErr_SetString(PyExc_ValueError, "toas and snrs must have the same length");
        goto fail;
    }
    snrs = PyArray_DATA(snrs_npy);

    toa_variances_npy = AsArrayObj(toa_variances_obj, NPY_DOUBLE, 1);
    if (!toa_variances_npy) goto fail;
    if (PyArray_DIM(toa_variances_npy, 0) != nifos)
    {
        PyErr_SetString(PyExc_ValueError, "toas and toa_variances must have the same length");
        goto fail;
    }
    toa_variances = PyArray_DATA(toa_variances_npy);

    responses_npy = malloc(nifos * sizeof(PyObject *));
    if (!responses_npy)
    {
        PyErr_SetNone(PyExc_MemoryError);
        goto fail;
    }
    for (i = 0; i < nifos; i ++)
        responses_npy[i] = NULL;
    responses = malloc(nifos * sizeof(float *));
    if (!responses)
    {
        PyErr_SetNone(PyExc_MemoryError);
        goto fail;
    }

    n = PySequence_Length(responses_obj);
    if (n < 0) goto fail;
    if (n != nifos)
    {
        PyErr_SetString(PyExc_ValueError, "toas and responses must have the same length");
        goto fail;
    }
    for (i = 0; i < nifos; i ++)
    {
        PyObject *obj = PySequence_GetItem(responses_obj, i);
        if (!obj) goto fail;
        responses_npy[i] = AsArrayObj(obj, NPY_FLOAT, 2);
        Py_XDECREF(obj);
        if (!responses_npy[i]) goto fail;
        if (PyArray_DIM(responses_npy[i], 0) != 3 || PyArray_DIM(responses_npy[i], 1) != 3)
        {
            PyErr_SetString(PyExc_ValueError, "expected every element of responses to be a 3x3 matrix");
            goto fail;
        }
        responses[i] = PyArray_DATA(responses_npy[i]);
    }

    locations_npy = malloc(nifos * sizeof(PyObject *));
    if (!locations_npy)
    {
        PyErr_SetNone(PyExc_MemoryError);
        goto fail;
    }
    for (i = 0; i < nifos; i ++)
        locations_npy[i] = NULL;
    locations = malloc(nifos * sizeof(double *));
    if (!locations)
    {
        PyErr_SetNone(PyExc_MemoryError);
        goto fail;
    }

    n = PySequence_Length(locations_obj);
    if (n < 0) goto fail;
    if (n != nifos)
    {
        PyErr_SetString(PyExc_ValueError, "toas and locations must have the same length");
        goto fail;
    }
    for (i = 0; i < nifos; i ++)
    {
        PyObject *obj = PySequence_GetItem(locations_obj, i);
        if (!obj) goto fail;
        locations_npy[i] = AsArrayObj(obj, NPY_DOUBLE, 1);
        Py_XDECREF(obj);
        if (!locations_npy[i]) goto fail;
        if (PyArray_DIM(locations_npy[i], 0) != 3)
        {
            PyErr_SetString(PyExc_ValueError, "expected every element of locations to be a vector of length 3");
            goto fail;
        }
        locations[i] = PyArray_DATA(locations_npy[i]);
    }

    horizons_npy = AsArrayObj(horizons_obj, NPY_DOUBLE, 1);
    if (!horizons_npy) goto fail;
    if (PyArray_DIM(horizons_npy, 0) != nifos)
    {
        PyErr_SetString(PyExc_ValueError, "toas and horizons must have the same length");
        goto fail;
    }
    horizons = PyArray_DATA(horizons_npy);

    if (prior_str)
    {
        if (strcmp(prior_str, "uniform in log distance"))
            prior = BAYESTAR_PRIOR_UNIFORM_IN_LOG_DISTANCE;
        else if (strcmp(prior_str, "uniform in volume"))
            prior = BAYESTAR_PRIOR_UNIFORM_IN_VOLUME;
    }

    old_handler = gsl_set_error_handler(my_gsl_error);
    result = bayestar_sky_map_tdoa_snr(npix, P, gmst, nifos, responses, locations, toas, snrs, toa_variances, horizons, min_distance, max_distance, prior);
    gsl_set_error_handler(old_handler);
    if (result != GSL_SUCCESS)
        goto fail;

    ret = out;
    out = NULL;
fail:
    Py_XDECREF(toas_npy);
    Py_XDECREF(snrs_npy);
    Py_XDECREF(toa_variances_npy);
    if (responses_npy)
        for (i = 0; i < nifos; i ++)
            Py_XDECREF(responses_npy[i]);
    free(responses_npy);
    free(responses);
    if (locations_npy)
        for (i = 0; i < nifos; i ++)
            Py_XDECREF(locations_npy[i]);
    free(locations_npy);
    free(locations);
    Py_XDECREF(horizons_npy);
    Py_XDECREF(out);
    return ret;
};


static PyMethodDef methods[] = {
    {"sky_map_tdoa", (PyCFunction)sky_map_tdoa, METH_VARARGS | METH_KEYWORDS, "fill me in"},
    {"sky_map_tdoa_snr", (PyCFunction)sky_map_tdoa_snr, METH_VARARGS | METH_KEYWORDS, "fill me in"},
    {NULL, NULL, 0, NULL}
};


PyMODINIT_FUNC
init_sky_map(void) {
    Py_InitModule("_sky_map", methods);
    import_array();
}
