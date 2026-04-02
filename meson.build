project(
    'dwave-samplers',
    'c', 'cpp', 'cython',
    default_options: [
        'cpp_std=c++17',
        'buildtype=release',
        'b_ndebug=if-release',  # add -DNDEBUG when building for "release", which we do by default
        'optimization=3',
        'pkgconfig.relocatable=true',
        'warning_level=1',  # equivalent to -Wall. We should consider cranking it up to 3
        'python.allow_limited_api=false',  # Our build system can then opt-in for Python 3.12+
    ],
    meson_version: '>= 1.3.0',
    # As of March 2026, there is no nicer way that I know of to dynamically get the version
    version: run_command(
        ['dwave/samplers/_build/_version.py'],
        capture: true,
        check: true,
    ).stdout(),
)

add_project_arguments(['-g1'], language: ['cpp'])

# Get access to Python
py = import('python').find_installation(pure: false)

# We build the limited API when we have Python>=3.12, which is the lowest
# version that Cython makes use of vectorcall.
# See https://cython.readthedocs.io/en/stable/src/userguide/limited_api.html
# By default this will do nothing unless the user overrides the
# python.allow_limited_api=false we specify above, which we do when building
# with cibuildwheel
if py.version().version_compare('>=3.12')
    limited_api = '3.12'
else
    limited_api = ''  # i.e. no limited api
endif

# We need the include files from NumPy and/or dimod for some of our Cython extensions.
# So let's get them and save the information for later.
dimod_dep = declare_dependency(
    include_directories: run_command(
        py,
        ['-c', 'import os.path; import dimod; print(os.path.relpath(dimod.get_include()))'],
      ).stdout().strip(),
)
numpy_dep = declare_dependency(
    include_directories: run_command(
        py,
        ['-c', 'import os.path; import numpy as np; print(os.path.relpath(np.get_include()))'],
      ).stdout().strip(),
)

# Ok, now we can build each of our samplers one-by-one. They all have some idiosyncratic
# build needs so we do them individually.

# -- Greedy
greedy_dep = declare_dependency(
    include_directories: ['dwave/samplers/greedy/src/'],
    sources: ['dwave/samplers/greedy/src/descent.cpp'],
)
py.extension_module(
    'descent',
    'dwave/samplers/greedy/descent.pyx',
    dependencies: [greedy_dep, numpy_dep],
    install: true,
    override_options : ['cython_language=cpp'],
    cython_args : ['-Xfreethreading_compatible=True'],
    subdir: 'dwave/samplers/greedy',
    limited_api: limited_api,
)

# -- Random
py.extension_module(
    'cyrandom',
    'dwave/samplers/random/cyrandom.pyx',
    dependencies: [dimod_dep, numpy_dep],
    install: true,
    override_options : ['cython_language=cpp'],
    cython_args : ['-Xfreethreading_compatible=True'],
    subdir: 'dwave/samplers/random',
    limited_api: limited_api,
)

# -- Simulated Annealing
sa_dep = declare_dependency(
    include_directories: ['dwave/samplers/sa/src/'],
    sources: ['dwave/samplers/sa/src/cpu_sa.cpp'],
)
py.extension_module(
    'simulated_annealing',
    'dwave/samplers/sa/simulated_annealing.pyx',
    dependencies: [numpy_dep, sa_dep],
    install: true,
    override_options : ['cython_language=cpp'],
    cython_args : ['-Xfreethreading_compatible=True'],
    subdir: 'dwave/samplers/sa',
    limited_api: limited_api,
)

# -- Simulated Quantum Annealing
sqa_dep = declare_dependency(
    include_directories: ['dwave/samplers/sqa/src/'],
    sources: [
        'dwave/samplers/sqa/src/cpu_rotormc.cpp',
        'dwave/samplers/sqa/src/localPIMC.cpp',
    ],
)
py.extension_module(
    'pimc_annealing',
    'dwave/samplers/sqa/pimc_annealing.pyx',
    dependencies: [numpy_dep, sqa_dep],
    install: true,
    override_options : ['cython_language=cpp'],
    cython_args : ['-Xfreethreading_compatible=True'],
    subdir: 'dwave/samplers/sqa',
    limited_api: limited_api,
)
py.extension_module(
    'rotormc_annealing',
    'dwave/samplers/sqa/rotormc_annealing.pyx',
    dependencies: [numpy_dep, sqa_dep],
    install: true,
    override_options : ['cython_language=cpp'],
    cython_args : ['-Xfreethreading_compatible=True'],
    subdir: 'dwave/samplers/sqa',
    limited_api: limited_api,
)

# -- Tabu
tabu_dep = declare_dependency(
    include_directories: ['dwave/samplers/tabu/src/'],
    sources: [
        'dwave/samplers/tabu/src/bqp.cpp',
        'dwave/samplers/tabu/src/tabu_search.cpp',
        'dwave/samplers/tabu/src/tabu_utils.cpp',
    ],
)
py.extension_module(
    'tabu_search',
    'dwave/samplers/tabu/tabu_search.pyx',
    dependencies: [tabu_dep],
    install: true,
    override_options : ['cython_language=cpp'],
    cython_args : ['-Xfreethreading_compatible=True'],
    subdir: 'dwave/samplers/tabu',
    limited_api: limited_api,
)

# -- Tree
tree_dep = declare_dependency(
    include_directories: ['dwave/samplers/tree/src/include/'],
)
py.extension_module(
    'sample',
    'dwave/samplers/tree/sample.pyx',
    dependencies: [dimod_dep, numpy_dep, tree_dep],
    install: true,
    override_options : ['cython_language=cpp'],
    cython_args : ['-Xfreethreading_compatible=True'],
    subdir: 'dwave/samplers/tree',
    limited_api: limited_api,
)
py.extension_module(
    'solve',
    'dwave/samplers/tree/solve.pyx',
    dependencies: [dimod_dep, numpy_dep, tree_dep],
    install: true,
    override_options : ['cython_language=cpp'],
    cython_args : ['-Xfreethreading_compatible=True'],
    subdir: 'dwave/samplers/tree',
    limited_api: limited_api,
)
py.extension_module(
    'utilities',
    'dwave/samplers/tree/utilities.pyx',
    dependencies: [dimod_dep],
    install: true,
    override_options : ['cython_language=cpp'],
    cython_args : ['-Xfreethreading_compatible=True'],
    subdir: 'dwave/samplers/tree',
    limited_api: limited_api,
)

install_subdir('dwave/samplers', install_dir: py.get_install_dir(pure: false) / 'dwave')

# meson doesn't disable tests by default, so we let the user decide via the build-tests
# option. If set to 'auto', then we build the tests for non-release builds.
if get_option('build-tests').enabled() or (get_option('build-tests').auto() and get_option('buildtype') != 'release')
    subdir('tests/')
endif
