# Copyright 2019 D-Wave Systems Inc.
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
from itertools import permutations


from dwave.samplers.stopwatch import (Stopwatch, NonmonotonicTimestampError,
                                      RepeatedTimestampError, MissingTimestampError)


class TestStopwatch(unittest.TestCase):
    # Validity of some tests depend on this order; do not change list order!
    TIMING_METHOD_NAMES = ["start_preprocessing", "start_sampling",
                           "start_postprocessing", "end_postprocessing"]

    def test_missing_timestamps(self):
        # Skipping any one timestamp should raise an error
        num_timestamps = 4
        for skip_index in range(num_timestamps):
            sw = Stopwatch()
            for method_index in range(num_timestamps):
                if method_index == skip_index:
                    continue
                timing_method = getattr(sw, self.TIMING_METHOD_NAMES[method_index])
                timing_method()
            self.assertRaises(MissingTimestampError, sw.report)

    def test_monotonicity_and_report(self):
        num_timestamps = 4

        # Any nonmonotonic call order should raise an error
        # otherwise, the check all keys are present and that the timings are nonnegative
        invocation_orders = permutations(range(num_timestamps))
        for invocation_order in invocation_orders:
            sw = Stopwatch()

            timing_methods = [getattr(sw, method) for method in self.TIMING_METHOD_NAMES]
            for method_index in invocation_order:
                timing_methods[method_index]()

            monotonic = invocation_order == (0, 1, 2, 3)
            if not monotonic:
                self.assertRaises(NonmonotonicTimestampError, sw.report)
            else:
                report = sw.report()
                self.assertIn("timing", report)

                for category, timing in report['timing'].items():
                    self.assertGreaterEqual(timing, 0)

                timings = report['timing']
                self.assertIn("preprocessing_ns", timings)
                self.assertIn("sampling_ns", timings)
                self.assertIn("postprocessing_ns", timings)

    def test_repeated_timestamps(self):
        sw = Stopwatch()

        # Check everything initializes to None (an assumption)
        num_timestamps = 4
        self.assertListEqual([None] * num_timestamps,
                             [sw.timestamp_preprocessing, sw.timestamp_sampling,
                              sw.timestamp_postprocessing, sw.timestamp_end])

        timing_methods = [getattr(sw, method) for method in self.TIMING_METHOD_NAMES]

        # Timestamp each process
        for timing_method in timing_methods:
            timing_method()

        # Check timestamps have been applied
        self.assertNotIn(None,
                         [sw.timestamp_preprocessing, sw.timestamp_sampling,
                          sw.timestamp_postprocessing, sw.timestamp_end])

        # Timestamping again should raise an error
        for timing_method in timing_methods:
            self.assertRaises(RepeatedTimestampError, timing_method)


if __name__ == "__main__":
    unittest.main()
