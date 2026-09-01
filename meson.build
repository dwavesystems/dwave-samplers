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
    meson_version: '>= 1.10.0',
    # As of March 2026, there is no nicer way that I know of to dynamically get the version
    version: run_command(
        ['dwave/samplers/_build/_version.py'],
        capture: true,
        check: true,
    ).stdout(),
)

# Configure and/or get access to our languages
py = import('python').find_installation(pure: false)

cpp = meson.get_compiler('cpp')
if cpp.has_argument('-g1')
    add_project_arguments('-g1', language: 'cpp')
endif

# Declare our C++ libraries as dependencies (for those samplers that have one).
# We want these each to be a separate dependency so that other packages can use
# dwave-samplers as a subproject to build against them.
greedy_dep = declare_dependency(
    include_directories: ['dwave/samplers/greedy/src/'],
    sources: ['dwave/samplers/greedy/src/descent.cpp'],
)
sa_dep = declare_dependency(
    include_directories: ['dwave/samplers/sa/src/'],
    sources: ['dwave/samplers/sa/src/cpu_sa.cpp'],
)
sbm_dep = declare_dependency(
    include_directories: ['dwave/samplers/sbm/src/'],
    sources: ['dwave/samplers/sbm/src/cpu_dsbm.cpp'],
)
sqa_dep = declare_dependency(
    include_directories: ['dwave/samplers/sqa/src/'],
    sources: [
        'dwave/samplers/sqa/src/cpu_rotormc.cpp',
        'dwave/samplers/sqa/src/localPIMC.cpp',
    ],
)
tabu_dep = declare_dependency(
    include_directories: ['dwave/samplers/tabu/src/'],
    sources: [
        'dwave/samplers/tabu/src/bqp.cpp',
        'dwave/samplers/tabu/src/tabu_search.cpp',
        'dwave/samplers/tabu/src/tabu_utils.cpp',
    ],
)
tree_dep = declare_dependency(
    include_directories: ['dwave/samplers/tree/src/include/'],
)

# Now the Cython modules that are used in the Python package.
# The build-python option toggles this off to save time for packages that
# use dwave-samplers as a subproject for the C++ dependencies.
if get_option('build-python')
    extension_kwargs = {
        'install': true,
        'override_options' : ['cython_language=cpp'],
        'cython_args' : ['-Xfreethreading_compatible=True'],
    }

    # We build the limited API when we have Python>=3.12, which is the lowest
    # version that Cython makes use of vectorcall.
    # See https://cython.readthedocs.io/en/stable/src/userguide/limited_api.html
    # By default this will do nothing unless the user overrides the
    # python.allow_limited_api=false we specify above, which we do when building
    # with cibuildwheel
    if py.version().version_compare('>=3.12')
        extension_kwargs += {'limited_api': '3.12'}
    endif

    # We need the include files from dimod for some of our Cython extensions.
    dimod_dep = subproject(
        'dimod',
        default_options: {
            'build-python': false,  # We only want the dependency
            'build-tests': 'disabled',  # Likewise we don't need to build the tests
        },
    ).get_variable('dimod_dep')

    # Ok, now we can build each of our samplers one-by-one. They all have some idiosyncratic
    # build needs so we do them individually.

    # -- Greedy
    py.extension_module(
        'descent',
        'dwave/samplers/greedy/descent.pyx',
        dependencies: [greedy_dep],
        subdir: 'dwave/samplers/greedy',
        kwargs: extension_kwargs,
    )

    # -- Random
    py.extension_module(
        'cyrandom',
        'dwave/samplers/random/cyrandom.pyx',
        dependencies: [dimod_dep],
        subdir: 'dwave/samplers/random',
        kwargs: extension_kwargs,
    )

    # -- Simulated Annealing
    py.extension_module(
        'simulated_annealing',
        'dwave/samplers/sa/simulated_annealing.pyx',
        dependencies: [sa_dep],
        subdir: 'dwave/samplers/sa',
        kwargs: extension_kwargs,
    )

    # -- Discrete Simulated Bifurcation Machine
    py.extension_module(
        'dsbm',
        'dwave/samplers/sbm/dsbm.pyx',
        dependencies: [sbm_dep],
        subdir: 'dwave/samplers/sbm',
        kwargs: extension_kwargs,
    )

    # -- Simulated Quantum Annealing
    py.extension_module(
        'pimc_annealing',
        'dwave/samplers/sqa/pimc_annealing.pyx',
        dependencies: [sqa_dep],
        subdir: 'dwave/samplers/sqa',
        kwargs: extension_kwargs,
    )
    py.extension_module(
        'rotormc_annealing',
        'dwave/samplers/sqa/rotormc_annealing.pyx',
        dependencies: [sqa_dep],
        subdir: 'dwave/samplers/sqa',
        kwargs: extension_kwargs,
    )

    # -- Tabu
    py.extension_module(
        'tabu_search',
        'dwave/samplers/tabu/tabu_search.pyx',
        dependencies: [tabu_dep],
        subdir: 'dwave/samplers/tabu',
        kwargs: extension_kwargs,
    )

    # -- Tree
    py.extension_module(
        'sample',
        'dwave/samplers/tree/sample.pyx',
        dependencies: [dimod_dep, tree_dep],
        subdir: 'dwave/samplers/tree',
        kwargs: extension_kwargs,
    )
    py.extension_module(
        'solve',
        'dwave/samplers/tree/solve.pyx',
        dependencies: [dimod_dep, tree_dep],
        subdir: 'dwave/samplers/tree',
        kwargs: extension_kwargs,
    )
    py.extension_module(
        'utilities',
        'dwave/samplers/tree/utilities.pyx',
        dependencies: [dimod_dep],
        subdir: 'dwave/samplers/tree',
        kwargs: extension_kwargs,
    )

    # If we're not building the extensions, we shouldn't install the package
    # because imports will fail.
    install_subdir('dwave/samplers', install_dir: py.get_install_dir(pure: false) / 'dwave')

endif  # get_option('build-python')

# meson doesn't disable tests by default, so we let the user decide via the build-tests
# option. If set to 'auto', then we build the tests for non-release builds.
if get_option('build-tests').enabled() or (get_option('build-tests').auto() and get_option('buildtype') != 'release')
    subdir('tests/')
endif
