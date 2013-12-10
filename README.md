# orang solver

Generic tree decomposition-based solver.

## MATLAB usage

### Building

Boost is required to build the MATLAB commands.  On Ubuntu, you can install
Boost like this:

    $ sudo apt-get install libboost-dev

Run the `buildOrang` command in the `matlab` subdirectory to build.

### Running

There are four main functions available:

* `orang_minsum` for optimization
* `orang_sample` for Boltzmann sampling
* `orang_mincount` for counting minimum-energy states
* `orang_greedyvarorder` for construction variable elimination orders

MATLAB documentation is available for all of them.

