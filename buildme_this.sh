module purge
module load shared
module load cmake/3.11.3
module load openmpi/gcc/4.0.1
module load hdf5/gcc_par/1.10.5

which g++
which gcc

cmake \
    -DCMAKE_INSTALL_PREFIX=$HOME/PSC_Install \
    -DCMAKE_CXX_COMPILER=/cm/shared/apps/gcc/6.4/bin/g++ \
    -DCMAKE_C_COMPILER=/cm/shared/apps/gcc/6.4/bin/gcc \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DCMAKE_C_FLAGS_RELWITHDEBINFO="-g -O2" \
    -DCMAKE_CXX_FLAGS_RELWITHDEBINFO="-g -O2" \
    -DUSE_CUDA=OFF \
    -DUSE_VPIC=OFF \
    ..
