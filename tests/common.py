# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

def ok_(a, s=None):
    if s is None:
        assert a
    else:
        assert a, s

def eq_(a, b, s=None):
    if s is None:
        assert a == b
    else:
        assert a == b, s

def neq_(a, b, s=None):
    if s is None:
        assert a != b
    else:
        assert a != b, s

def gt_(a, b, s=None):
    if s is None:
        assert a > b
    else:
        assert a > b, s

def lt_(a, b, s=None):
    if s is None:
        assert a < b
    else:
        assert a < b, s

def gte_(a, b, s=None):
    if s is None:
        assert a >= b
    else:
        assert a >= b, s

def lte_(a, b, s=None):
    if s is None:
        assert a <= b
    else:
        assert a <= b, s

def almost_eq_(a, b, tol=1e-6, s=None):
    if s is None:
        assert abs(a - b) < tol
    else:
        assert abs(a - b) < tol, s