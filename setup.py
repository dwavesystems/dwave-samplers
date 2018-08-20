from __future__ import absolute_import

import sys
import os

from setuptools import setup

# add __version__, __author__, __authoremail__, __description__ to this namespace
_PY2 = sys.version_info.major == 2
my_loc = os.path.dirname(os.path.abspath(__file__))
os.chdir(my_loc)
if _PY2:
    execfile(os.path.join(".", "savanna", "package_info.py"))
else:
    exec(open(os.path.join(".", "savanna", "package_info.py")).read())

install_requires = ['dimod>=0.6.9,<0.8.0',
                    'networkx>=2.0,<3.0',
                    'scipy>=1.1.0,<2.0.0',
                    ]

packages = ['savanna',
            'savanna.io',
            ]

classifiers = [
    'Operating System :: OS Independent',
    'Programming Language :: Python :: 2',
    'Programming Language :: Python :: 2.7',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.4',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',
    'Programming Language :: Python :: 3.7',
    ]

python_requires = '>=2.7,!=3.0.*,!=3.1.*,!=3.2.*,!=3.3.*'

setup(
    name='savanna',
    version=__version__,
    author=__author__,
    author_email=__authoremail__,
    description=__description__,
    long_description=open('README.rst').read(),
    url='https://github.com/dwavesystems/savanna',
    download_url='https://github.com/dwavesystems/savanna/releases',
    # license='',
    packages=packages,
    install_requires=install_requires,
    classifiers=classifiers,
    python_requires=python_requires
)
