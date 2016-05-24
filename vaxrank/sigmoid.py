# Copyright (c) 2016. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import numpy as np

class DecreasingSigmoid(object):
    def __init__(self, midpoint, width):
        assert width > 0
        self.midpoint = midpoint
        self.width = width

    def __call__(self, value):
        normalized = (float(value) - self.midpoint) / self.width
        # simplification of 1.0 - logistic(x) = logistic(-x)
        return 1.0 / (1.0 + np.exp(normalized))

class IncreasingSigmoid(object):
    def __init__(self, midpoint, width):
        assert width > 0
        self.midpoint = midpoint
        self.width = width

    def __call__(self, value):
        normalized = (float(value) - self.midpoint) / self.width
        return 1.0 / (1.0 + np.exp(-normalized))

