

```
# built for MacOS
git clone --recursive https://github.com/tt6746690/sparse_solver.git
brew install libomp

mkdir build && cd build
cmake -DWITH_DEBUG=OFF ..
make

# unit tests
cd .. && ./build/tests
./build/tests "triangular_large"    # long running tests

# debugging
lldb build/main
process launch
```