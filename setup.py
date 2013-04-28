#
# Copyright (C) 2013  Leo Singer
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#


from setuptools import setup
from misc import *
from misc.setuptools_openmp import *
import numpy as np


setup(
    name='bayestar-localization',
    version='0.0.6',
    author='Leo Singer',
    author_email='leo.singer@ligo.org',
    description='Rapid Bayesian sky maps for gravitational wave inspirals',
    license='GNU GPLv2+',
    requires=['bayestar', 'healpy', 'numpy', 'glue', 'pylal', 'lal', 'lalsimulation'],
    namespace_packages=['bayestar'],
    packages=['bayestar'],
    ext_modules=[
        Extension('bayestar.sky_map', ['bayestar/sky_map.c', 'bayestar/bayestar_sky_map.c'],
            **copy_library_dirs_to_runtime_library_dirs(
            **pkgconfig('lal', 'lalsimulation', 'gsl', 'chealpix',
                include_dirs=[np.get_include()],
                extra_compile_args=['-std=c99'],
                define_macros=[('HAVE_INLINE', None)],
                openmp=True
            ))
        )
    ],
    scripts=[
        'bin/bayestar_aggregate_found_injections',
        'bin/bayestar_plot_found_injections',
        'bin/bayestar_lattice_tmpltbank',
        'bin/bayestar_localize_gracedb',
        'bin/bayestar_localize_lvalert',
        'bin/bayestar_cluster_coincs',
        'bin/bayestar_prune_neighborhood_tmpltbank',
        'bin/bayestar_realize_coincs',
        'bin/bayestar_localize_coincs',
        'bin/bayestar_sim_to_tmpltbank',
        'bin/ligolw_coire_to_coinc',
        'bin/bayestar_littlehope'
    ],
    cmdclass={'build_ext': build_ext}
)
