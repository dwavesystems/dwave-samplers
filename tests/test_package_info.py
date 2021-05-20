# Copyright 2021 D-Wave Systems Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import unittest

import greedy


class TestSteepestDescentSampler(unittest.TestCase):

    def test_import(self):
        from greedy import package_info, __version__

        self.assertEqual(__version__, package_info.__version__)
        self.assertEqual(__version__, greedy.package_info.__version__)
        self.assertIn('__version__', dir(package_info))

    def test_attributes(self):
        pi = greedy.package_info

        self.assertIsNotNone(pi.__package_name__)
        self.assertIsNotNone(pi.__version__)
        self.assertIsNotNone(pi.__author__)
        self.assertIsNotNone(pi.__author_email__)
        self.assertIsNotNone(pi.__description__)
        self.assertIsNotNone(pi.__url__)
        self.assertIsNotNone(pi.__license__)
