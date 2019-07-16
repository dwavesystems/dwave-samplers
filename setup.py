import os
from io import open
from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext


# Load package info, without importing the package
basedir = os.path.dirname(os.path.abspath(__file__))
package_info_path = os.path.join(basedir, "greedy", "package_info.py")
package_info = {}
try:
    with open(package_info_path, encoding='utf-8') as f:
        exec(f.read(), package_info)
except SyntaxError:
    execfile(package_info_path, package_info)


packages = ['greedy']

# Package requirements, minimal pinning
install_requires = ['six>=1.10', 'numpy>=1.15.0,<1.16.0', 'dimod>=0.8.11,<0.9.0']

# Setup (extension build) requirements
setup_requires = ['numpy>=1.14.0,<1.16.0']

# Package extras requirements
extras_require = {
    'test': ['coverage', 'mock'],
}


# Custom build_ext subclass to add compile/link args and numpy include dirs
extra_compile_args = {
    'msvc': ['/std:c++14'],
    'unix': ['-std=c++11'],
}

extra_link_args = {
    'msvc': [],
    'unix': [],
}

class build_ext_compiler_check(build_ext):
    def build_extensions(self):
        compiler = self.compiler.compiler_type

        compile_args = extra_compile_args[compiler]
        for ext in self.extensions:
            ext.extra_compile_args = compile_args

        link_args = extra_compile_args[compiler]
        for ext in self.extensions:
            ext.extra_compile_args = link_args

        build_ext.build_extensions(self)

    def run(self):
        import numpy
        self.include_dirs.append(numpy.get_include())
        build_ext.run(self)


# We distribute cythonized source, so cython is not required
# for install (build) from package
if not os.path.exists(os.path.join(basedir, 'PKG-INFO')):
    try:
        from Cython.Build import cythonize
        USE_CYTHON = True
    except ImportError:
        USE_CYTHON = False
else:
    USE_CYTHON = False

ext = '.pyx' if USE_CYTHON else '.cpp'

extensions = [Extension(
    name='greedy.descent',
    sources=['./greedy/descent' + ext],
)]

if USE_CYTHON:
    extensions = cythonize(extensions)


classifiers = [
    'License :: OSI Approved :: Apache Software License',
    'Operating System :: OS Independent',
    'Development Status :: 3 - Alpha',
    'Programming Language :: Python :: 2',
    'Programming Language :: Python :: 2.7',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',
    'Programming Language :: Python :: 3.7',
]

setup(
    name=package_info['__packagename__'],
    version=package_info['__version__'],
    author=package_info['__author__'],
    author_email=package_info['__authoremail__'],
    description=package_info['__description__'],
    long_description=open('README.rst', encoding='utf-8').read(),
    url=package_info['__url__'],
    license=package_info['__license__'],
    packages=packages,
    install_requires=install_requires,
    setup_requires=setup_requires,
    extras_require=extras_require,
    ext_modules=extensions,
    cmdclass={'build_ext': build_ext_compiler_check},
    classifiers=classifiers,
    zip_safe=False,
)
