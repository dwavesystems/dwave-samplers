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

__package_name__ = 'dwave-greedy'
__version__ = '0.2.2'

from greedy.sampler import *
import greedy.sampler

from greedy.composite import *
import greedy.composite

# Note: provide access to `greedy.package_info` for backwards compatibility, but
# (lazy-)load values from setup.cfg. We load lazily to save time during `greedy`
# import, if `package_info` is not requested (package metadata has to be loaded
# from the filesystem).
#
# The workaround below can be replaced with module-level `__getattr__` (PEP-562)
# once we drop support for py36.
import sys
import types
import warnings

class _PackageInfoModule(types.ModuleType):

    # map module attribute name to metadata entry name
    _attr_map = {
        '__package_name__': 'name',
        '__version__': 'version',
        '__author__': 'author',
        '__author_email__': 'author-email',
        '__description__': 'summary',
        '__url__': 'home-page',
        '__license__': 'license',
    }

    class _LazyAttributeLoaderClassProperty:
        # roughly equivalent to ``classmethod(property(cached_attrs))``, but it
        # doesn't require chained decorators support available in py39+

        def __init__(self, attr_map):
            self._attr_map = attr_map
            self._attrs = None

        def _metadata_attributes(self, attr_map):
            try:
                from importlib import metadata
            except ImportError:
                import importlib_metadata as metadata

            m = metadata.metadata(__package_name__)

            return {attr: m.get(meta) for attr, meta in attr_map.items()}

        def __get__(self, obj, objtype=None):
            if self._attrs is None:
                self._attrs = self._metadata_attributes(self._attr_map)
            return self._attrs

    _attrs = _LazyAttributeLoaderClassProperty(_attr_map)

    def __dir__(self):
        return tuple(self._attr_map.keys())

    def __getattr__(self, name):
        try:
            val = self._attrs[name]
            warnings.warn(
                f"'{self.__name__}' module is deprecated. "
                "Please convert your code to use `importlib.metadata` instead.",
                DeprecationWarning, stacklevel=2)
            return val
        except KeyError:
            raise AttributeError(f"module '{self.__name__}' has no attribute '{name}'")

def _create_package_info_module():
    name = 'package_info'
    import_path = f'{__name__}.{name}'

    globals()[name] = sys.modules[import_path] = _PackageInfoModule(import_path)

# note: create both module and local attribute
_create_package_info_module()
