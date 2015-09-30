import os
import sys
from pkg_resources import Distribution, get_build_platform
import zipfile


PKG_INFO = (
    ('Metadata-Version', '1.0'),
    ('Name', 'orang'),
    ('Version', '1.0'),
    ('Summary', 'orang'),
    ('Home-page', 'http://www.dwavesys.com'),
    ('Author', 'D-Wave Systems Inc.'),
    ('Author-email', 'dw1support@dwavesys.com'),
    ('License', 'other'),
    ('Platform', 'POSIX, MacOS, Windows')
)

SWIG_MODULES = ['orang']

EGG_NAME = Distribution(
    project_name=dict(PKG_INFO)['Name'],
    version=dict(PKG_INFO)['Version'],
    platform=get_build_platform()
).egg_name() + '.egg'

EXT_EXT = '.pyd' if os.name == 'nt' else '.so'


def build_egg():
    egg = zipfile.ZipFile(EGG_NAME, mode='w', compression=zipfile.ZIP_DEFLATED)
    add_egg_info(egg)
    add_swig_modules(egg)


def add_egg_info(egg):
    # egg.writestr('EGG-INFO/dependency_links.txt', '')
    # egg.writestr('EGG-INFO/SOURCES.txt', '')
    egg.writestr('EGG-INFO/requires.txt', 'numpy\n')
    egg.writestr('EGG-INFO/not-zip-safe', '')
    egg.writestr('EGG-INFO/top_level.txt',
                 ''.join('{0}\n_{0}\n'.format(x) for x in SWIG_MODULES))
    egg.writestr('EGG-INFO/native_libs.txt',
                 ''.join('_{0}{1}'.format(x, EXT_EXT) for x in SWIG_MODULES))
    egg.writestr('EGG-INFO/PKG-INFO',
                 ''.join('{0}[0]: {0}[1]\n'.format(x) for x in PKG_INFO))


def add_swig_modules(egg):
    for module in SWIG_MODULES:
        egg.write('{0}.py'.format(module))
        egg.write('{0}.pyc'.format(module))
        egg.write('_{0}{1}'.format(module, EXT_EXT))


if __name__ == '__main__':
    if len(sys.argv) == 1:
        build_egg()
    elif sys.argv[1:] == ['--egg-name']:
        print EGG_NAME
    else:
        sys.stderr.write('Bad usage\n')
        sys.exit(1)

