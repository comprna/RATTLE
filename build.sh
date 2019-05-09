# build SPOA lib
cd spoa
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make

# build RATTLE
cd ../..
make