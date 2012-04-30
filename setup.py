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


from distutils.core import setup
from distutils.extension import Extension
import numpy as np
import os

healpixdir = os.getenv('HEALPIXDIR')
if healpixdir is None:
    healpix_include_dirs = []
    healpix_library_dirs = []
else:
    healpix_include_dirs = [os.path.join(healpixdir, 'include')]
    healpix_library_dirs = [os.path.join(healpixdir, 'lib')]


## {{{ http://code.activestate.com/recipes/502261/ (r1)
## But modified to throw an error if the package cannot be found.
def pkgconfig(*packages, **kw):
    import commands
    from distutils.errors import DistutilsExecError
    flag_map = {'-I': 'include_dirs', '-L': 'library_dirs', '-l': 'libraries'}
    status, output = commands.getstatusoutput("pkg-config --libs --cflags %s" % ' '.join(["'%s'" % s for s in packages]))
    if status != 0:
        raise DistutilsExecError("pkg-config: error %d\n%s" % (status, output))
    for token in output.split():
        kw.setdefault(flag_map.get(token[:2]), []).append(token[2:])
    return kw
## end of http://code.activestate.com/recipes/502261/ }}}


def copy_library_dirs_to_runtime_library_dirs(**kwargs):
    """Add every entry of library_dirs to runtime_library_dirs so that linker
    adds rpath entries."""
    kwargs.setdefault('runtime_library_dirs', []).extend(
        kwargs.get('library_dirs', []))
    return kwargs


setup(
    name='bayestar-localization',
    version='0.0.1',
    author='Leo Singer',
    author_email='leo.singer@ligo.org',
    license='GNU General Public License Version 3',
    packages=['bayestar_localization'],
    ext_modules=[
        Extension('bayestar_localization._sky_map', ['bayestar_localization/_sky_map.c', 'bayestar_localization/bayestar_sky_map.c'],
            **copy_library_dirs_to_runtime_library_dirs(
            **pkgconfig('lal', 'lalsimulation', 'gsl',
                include_dirs=[np.get_include()] + healpix_include_dirs,
                library_dirs=healpix_library_dirs,
                libraries=['cfitsio', 'chealpix', 'gomp'],
                extra_compile_args=['-fopenmp', '-std=c99']
            ))
        )
    ],
    scripts=['bin/bayestar_localization_fake_coincs'],
)
