GMXPATH="/your/gmx/src!"
GMXBUILD="${GMXPATH}/build"
GMXLOC="/desired/path/gromacs/gromacs_5_1"

GMXFLAGS="-DGMX_MPI=ON -DGMX_GPU=OFF"

ORIGDIR=${PWD}

#modules
module load open_mpi 
module load cmake/3.9.2
module load gsl
module load gcc
module load fftw

#setup building
mkdir ${GMXBUILD}
cd ${GMXBUILD}

cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON -DCMAKE_INSTALL_PREFIX=${GMXLOC} ${GMXFLAGS} 1>${ORIGDIR}/cmake.log 2>${ORIGDIR}/cmake.err || $(eval "echo -e "\n\n FAILED in cmake!!!\n\n" && exit 1")
echo  "finished Configure"
sleep 2
make 1> ${ORIGDIR}/make.log 2>${ORIGDIR}/make.err || $(eval "echo -e "\n\n FAILED in make!!!\n\n" && exit 1")

make check 1> ${ORIGDIR}/check.log 2> ${ORIGDIR}/check.err || $(eval "echo -e "\n\n FAILED in make check!!!\n\n" && exit 1")

make install 1> ${ORIGDIR}/install.log 2> ${ORIGDIR}/install.err || $(eval "echo -e "\n\n FAILED in make install!!!\n\n" && exit 1")



cd ${ORIGDIR}
(END)

