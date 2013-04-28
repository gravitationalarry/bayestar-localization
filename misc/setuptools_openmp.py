#
# Copyright (C) 2013  Leo Singer
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
"""
Enable OpenMP compiler support in Python C extensions. By default, your
`setup.py` script will check if the compiler supports OpenMP by compiling a test
program. If you specify the command line option `--enable-openmp`, then if the
compiler does not support OpenMP, then `setup.py` will fail with an error
message. If you specify the command line option `--disable-openmp`, then the
OpenMP support will be disabled unconditionally.


Examples
--------

    # Enable OpenMP support the compiler supports it
    $ python setup.py build_ext

    # Demand OpenMP support, fail if the compiler does not support it
    $ python setup.py build_ext --enable-openmp

    # Disable OpenMP support. Don't even check if the compiler supports it.
    $ python setup.py build_ext --disable-openmp


Usage
-----

To use, add the following line to your `setup.py` script::

    from misc.setuptools_openmp import *

Next, for any C extension for which you want optional OpenMP support, add the
`openmp=True` keyword argument to your invocation of `Extension()`::

    Extension(
        ...,
        openmp=True
    )

And finally, add the `cmdclass` keyword to your invocation of `setup()`::

    setup(
        ...
        cmdclass={'build_ext': build_ext},
    )

FIXME: Unify this with distutils_openmp.py.
"""
__author__ = "Leo Singer <leo.singer@ligo.org>"


from distutils import log as _log
from setuptools.extension import Extension as _Extension
from setuptools.command.build_ext import build_ext as _build_ext
from distutils.dir_util import mkpath as _mkpath
from distutils.errors import DistutilsOptionError as _DistutilsOptionError
from distutils.errors import DistutilsPlatformError as _DistutilsPlatformError
from distutils.errors import CCompilerError as _CCompilerError
import os as _os
import tempfile as _tempfile


class build_ext(_build_ext):
    """Extend the standard distutils build_ext command to compile and link with OpenMP support."""

    # Add the --enable-openmp and --disable-openmp command line arguments.
    user_options = _build_ext.user_options + [
        ('enable-openmp', None, 'enable OpenMP acceleration, or fail if compiler does not support it'),
        ('disable-openmp', None, 'disable OpenMP acceleration')
        ]

    # Mark the --enable-openmp and --disable-openmp command line arguments as booleans.
    boolean_options = _build_ext.boolean_options + ['enable-openmp',
        'disable-openmp']

    def initialize_options(self):
        # Initialize some variables.
        self.enable_openmp = 0
        self.disable_openmp = 0
        # Chain to method in parent class.
        _build_ext.initialize_options(self)

    def finalize_options(self):
        # Validate the --enable-openmp and --disable-openmp command line options:
        # the user should provide at most one of them.
        if self.enable_openmp and self.disable_openmp:
            raise _DistutilsOptionError("--enable-openmp and --disable-openmp are mutually exclusive")
        # Chain to method in parent class.
        _build_ext.finalize_options(self)

    def build_extensions(self):
        if self.disable_openmp:
            _log.info("disabling OpenMP support")
            use_openmp = False
        else:
            # Compile a test program to determine if compiler supports OpenMP.
            _mkpath(self.build_temp)
            with _tempfile.NamedTemporaryFile(mode='w', dir=self.build_temp, prefix='openmptest', suffix='.c') as srcfile:
                _log.info("checking if compiler supports OpenMP")
                srcfile.write("""
                #include <omp.h>
                int testfunc() {
                    int i;
                    #pragma omp parallel for
                    for (i = 0; i < 10; i ++)
                        ;
                    return omp_get_num_threads();
                }
                """)
                srcfile.flush()
                try:
                    objects = self.compiler.compile([srcfile.name], extra_postargs=["-fopenmp"], output_dir="/")
                except _CCompilerError:
                    if self.enable_openmp:
                        raise _DistutilsOptionError("compiler does not support OpenMP")
                    else:
                        _log.info("compiler does not support OpenMP")
                    use_openmp = False
                else:
                    _log.info("enabling OpenMP support")
                    use_openmp = True
                    for o in objects:
                        _os.remove(o)

        if use_openmp:
            for ext in self.extensions:
                if getattr(ext, 'openmp', False):
                    if ext.extra_compile_args:
                        ext.extra_compile_args += ['-fopenmp']
                    else:
                        ext.extra_compile_args = ['-fopenmp']
                    if ext.extra_link_args:
                        ext.extra_link_args += ['-fopenmp']
                    else:
                        ext.extra_link_args = ['-fopenmp']

        # Chain to method in parent class.
        _build_ext.build_extensions(self)


class Extension(_Extension):
    """Extend the standard distutils Extension to take an option for enabling OpenMP support."""
    def __init__(self, names, sources, openmp=False, **kw):
        self.openmp = openmp
        _Extension.__init__(self, names, sources, **kw)