catch2 = subproject(
    'catch2',
    default_options: {'tests': false},  # disable Catch2's self-tests
).get_variable('catch2_with_main_dep')

tests_greedy = executable(
    'tests-greedy',
    [
        'test_greedy.cpp',
    ],
    dependencies: [
        catch2,
        dep_greedy,
    ],
)
test('tests-greedy', tests_greedy, verbose: true)

tests_tabu = executable(
    'tests-tabu',
    [
        'test_tabu_bqp.cpp',
    ],
    dependencies: [
        catch2,
        dep_tabu,
    ],
)
test('tests-tabu', tests_tabu, verbose: true)
