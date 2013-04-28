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
Miscellaneous Distutils utilities.
"""
__author__ = "Leo Singer <leo.singer@ligo.org>"


import commands as _commands
from distutils.errors import DistutilsExecError as _DistutilsExecError


## {{{ http://code.activestate.com/recipes/502261/ (r1)
## But modified to throw an error if the package cannot be found.
def pkgconfig(*packages, **kw):
    flag_map = {'-I': 'include_dirs', '-L': 'library_dirs', '-l': 'libraries'}
    status, output = _commands.getstatusoutput("pkg-config --libs --cflags %s" % ' '.join(["'%s'" % s for s in packages]))
    if status != 0:
        raise _DistutilsExecError("pkg-config: error %d\n%s" % (status, output))
    for token in output.split():
        if token[:2] in flag_map:
            kw.setdefault(flag_map[token[:2]], []).append(token[2:])
    return kw
## end of http://code.activestate.com/recipes/502261/ }}}


def copy_library_dirs_to_runtime_library_dirs(**kwargs):
    """Add every entry of library_dirs to runtime_library_dirs so that linker
    adds rpath entries."""
    kwargs.setdefault('runtime_library_dirs', []).extend(
        kwargs.get('library_dirs', []))
    return kwargs