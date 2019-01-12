
```
# built for MacOS
git clone --recursive https://github.com/tt6746690/sparse_solver.git
brew install libomp

# compile
mkdir build && cd build
cmake ..
make

# need to put lhs for `torso` and `tsopf` inside `data/`
# rename `data/TSOPF_RS_b678_c2` -> `data/s_torsoL.mtx` and similarly for tsopf

# unit tests
cd .. && ./build/tests
./build/tests "triangular_large"    # long running tests on torso and tospof


// ./build/main matrix_name solver_type number_of_threads
//
//     matrix_name: `torso` will find `./data/{torso}L.mtx` and `./data/{torso}b.mtx`
//     solver_type: one of {all, simple, eigen, reachset, eigen_par}
//     n_threads  : a digit or a range, i.e. `3` or `1-4`
// 
./build/main small all 1-4

# debugging
lldb build/main
process launch
```