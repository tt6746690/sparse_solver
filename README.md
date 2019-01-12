
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

# debugging
lldb build/main
process launch
```