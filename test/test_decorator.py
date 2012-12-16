#!/usr/bin/env python
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
"""
Test cases for is_in_fov function.
"""
__author__ = "Leo Singer <leo.singer@ligo.org>"

import unittest
from bayestar_localization.decorator import *


class TestMemoize(unittest.TestCase):

    def test(self):

        @memoized
        class foo:
            def __init__(self, bar, bat, baz=1, baf=2):
                self.bar = bar
                self.bat = bat
                self.baz = baz
                self.baf = baf

            @memoized
            def xyzzy(self, plugh):
                return self.baz + plugh

        class bar:
            def __init__(self, bar, bat, baz=1, baf=2):
                self.bar = bar
                self.bat = bat
                self.baz = baz
                self.baf = baf

            @memoized
            def xyzzy(self, plugh):
                return self.baz + plugh

        foo1 = foo("bar", "bat")
        foo2 = foo("bar", "bat")
        foo3 = foo("bar", "baz")
        foo4 = foo("bar", "bat")
        foo5 = foo("bat", "bar")
        foo6 = foo("bar", "bat", baz=3)
        foo7 = foo("bar", "bat", baz=3)
        foo8 = foo("bar", "bat", baz=4)

        bar1 = bar("bar", "bat")
        bar2 = bar("bar", "bat")
        bar3 = bar("bar", "baz")
        bar4 = bar("bar", "bat")
        bar5 = bar("bat", "bar")
        bar6 = bar("bar", "bat", baz=3)
        bar7 = bar("bar", "bat", baz=3)
        bar8 = bar("bar", "bat", baz=4)

        # Check identity of memoized instances.
        self.assertIs(foo1, foo2)
        self.assertIsNot(foo1, foo3)
        self.assertIs(foo1, foo4)
        self.assertIsNot(foo1, foo5)
        self.assertIsNot(foo1, foo6)
        self.assertIs(foo6, foo7)
        self.assertIsNot(foo6, foo8)

        # Check that memoized instance methods do not alias instances.
        self.assertEqual(foo1.xyzzy(1), 2)
        self.assertEqual(foo6.xyzzy(1), 4)
        self.assertEqual(foo8.xyzzy(1), 5)
        self.assertEqual(bar1.xyzzy(1), 2)
        self.assertEqual(bar6.xyzzy(1), 4)
        self.assertEqual(bar8.xyzzy(1), 5)

if __name__ == '__main__':
    unittest.main()
