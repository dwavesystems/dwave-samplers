# Copyright 2023 D-Wave Systems Inc.
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

from time import perf_counter_ns as timer
from typing import Dict


class Stopwatch():
    """A stopwatch for timing samplers' preprocessing, sampling, and postprocessing time in units
    of nanoseconds. See documentation in each sampler for what each category entails.
    The order in which timing methods are invoked matter and should be:
    1. `start_preprocessing`; this marks the beginning of the preprocessing phase
    2. `start_sampling`; this marks the beginning of the sampling phase (end of preprocessing)
    3. `start_postprocessing`; this marks the beginning of the postprocessing phase (end of sampling)
    4. `stop_postprocessing`; this marks the end of the postprocessing phase
    """

    def __init__(self) -> None:
        self.timestamp_preprocessing = None
        self.timestamp_sampling = None
        self.timestamp_postprocessing = None
        self.timestamp_end = None

    def start_preprocessing(self) -> None:
        """Records the timestamp of the beginning of the preprocessing step.
        """
        if self.timestamp_preprocessing is not None:
            raise RepeatedTimestampError()
        self.timestamp_preprocessing = timer()

    def start_sampling(self) -> None:
        """Records the timestamp of the beginning of the sampling step.
        """
        if self.timestamp_sampling is not None:
            raise RepeatedTimestampError()
        self.timestamp_sampling = timer()

    def start_postprocessing(self) -> None:
        """Records the timestamp of the beginning of the postprocessing step.
        """
        if self.timestamp_postprocessing is not None:
            raise RepeatedTimestampError()
        self.timestamp_postprocessing = timer()

    def stop_postprocessing(self) -> None:
        """Records the final timestamp. This should be invoked last.
        """
        if self.timestamp_end is not None:
            raise RepeatedTimestampError()
        self.timestamp_end = timer()

    def report(self) -> Dict[str, float]:
        """Reports the duration of each process.

        - Preprocessing time is the duration from the beginning of preprocessing until the beginning
          of sampling.
        - Sampling time is the duration from the beginning of sampling until the beginning of
          postprocessing.
        - Postprocessing time is the duration from the beginning of postprocessing until the end of
          postprocessing.

        Returns:
            Dict[str: float]: timings of each category.
        """

        if None in {self.timestamp_preprocessing, self.timestamp_sampling,
                    self.timestamp_postprocessing, self.timestamp_end}:
            raise MissingTimestampError()

        monotonic = (self.timestamp_preprocessing
                     <= self.timestamp_sampling
                     <= self.timestamp_postprocessing
                     <= self.timestamp_end)

        if not monotonic:
            raise NonmonotonicTimestampError()

        return dict(timing=dict(
            preprocessing_ns=self.timestamp_sampling - self.timestamp_preprocessing,
            sampling_ns=self.timestamp_postprocessing - self.timestamp_sampling,
            postprocessing_ns=self.timestamp_end - self.timestamp_postprocessing
        ))


class RepeatedTimestampError(Exception):
    def __init__(self) -> None:
        super().__init__("Timestamp already exists")


class MissingTimestampError(Exception):
    def __init__(self) -> None:
        super().__init__(
            "Missing at least one timestamp; check all start and end times have been called."
        )


class NonmonotonicTimestampError(Exception):
    def __init__(self) -> None:
        super().__init__("Timings are non-monotonic; order of timestamps should be: "
                         + "preprocessing, sampling, postprocessing, then end.")
