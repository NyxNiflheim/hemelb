# Build the main application: source setup_hemeLB.sh

## unload conflicting modules
module purge
module use /mnt/lustre/e1000/home/y07/shared/cirrus-modulefiles
module load epcc/setup-env

module unload intel-20.4/compilers  # 先卸载 Intel 编译器包
module unload intel-20.4/cc intel-20.4/fc  # 再卸载 Intel C 和 Fortran 编译器
module unload gcc  # unload old version GCC
module unload hdf5parallel

## load required modules
module load gcc/12.3.0-offload  # load GCC
module load cmake/3.25.2
module load intel-19.5/mpi
module use /work/z04/shared/rwn/modules
module load boost/1.84.0
module load ctemplate/2.4
module load parmetis/4.0.3
module load tinyxml/2.6.2
module load catch2/2.13.6
## NEW: Load HDF5 module (compatible with GCC 12.3.0) 
module load hdf5parallel/1.14.1-2_cuda118_ompi414

## paths
source=/work/m24oc/m24oc/s2450341/hemelb/Code
build=/work/m24oc/m24oc/s2450341/hemelb/build
install=/work/m24oc/m24oc/s2450341/hemelb/install

rm -rf $build
mkdir -p $build
mkdir -p $install

## run CMake with explicit HDF5 paths
cmake -S $source -B $build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=$install \
  -DHEMELB_BUILD_TESTS=OFF
  -DHDF5_INCLUDE_DIR=/work/y07/shared/cirrus-software/hdf5parallel/1.14.3-gcc10.2.0-ompi4.1.6/include \
  -DHDF5_LIBRARY=/work/y07/shared/cirrus-software/hdf5parallel/1.14.3-gcc10.2.0-ompi4.1.6/lib/libhdf5.so \
  -DHDF5_HL_LIBRARY=/work/y07/shared/cirrus-software/hdf5parallel/1.14.3-gcc10.2.0-ompi4.1.6/lib/libhdf5_hl.so
# cmake -S $source -B $build \
#   -DCMAKE_BUILD_TYPE=Release \
#   -DCMAKE_INSTALL_PREFIX=$install \
#   -DHEMELB_BUILD_TESTS=OFF \
  # -DHDF5_INCLUDE_DIR=/work/y07/shared/cirrus-software/hdf5parallel/1.14.1-2-cuda11.8-ompi4.1.4/include \
  # -DHDF5_LIBRARY=/work/y07/shared/cirrus-software/hdf5parallel/1.14.1-2-cuda11.8-ompi4.1.4/lib/libhdf5.so \
  # -DHDF5_HL_LIBRARY=/work/y07/shared/cirrus-software/hdf5parallel/1.14.1-2-cuda11.8-ompi4.1.4/lib/libhdf5_hl.so

# cmake -S $source -B $build \
#   -DCMAKE_BUILD_TYPE=Release \
#   -DCMAKE_INSTALL_PREFIX=$install \
#   -DHDF5_INCLUDE_DIR=/work/y07/shared/cirrus-software/hdf5parallel/1.14.1-2-cuda11.8-ompi4.1.4/include \
#   -DHDF5_LIBRARY=/work/y07/shared/cirrus-software/hdf5parallel/1.14.1-2-cuda11.8-ompi4.1.4/lib/libhdf5.so \
#   -DHDF5_HL_LIBRARY=/work/y07/shared/cirrus-software/hdf5parallel/1.14.1-2-cuda11.8-ompi4.1.4/lib/libhdf5_hl.so

## Configure
# cmake -S $source -B $build -DCMAKE_BUILD_TYPE=Release
## Build and test
cmake --build $build -j 36
$build/tests/hemelb-tests

## Install
cmake -S /work/m24oc/m24oc/s2450341/hemelb/Code -B /work/m24oc/m24oc/s2450341/hemelb/build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$HOME/.local
cmake --install $build