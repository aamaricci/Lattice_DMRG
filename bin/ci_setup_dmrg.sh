#Building DMRG
#Errors
set -e

cd Lattice_DMRG
mkdir build
cd build

echo "cmake .."
cmake ..

echo "make"
make

echo "make install"
make install

echo "source ~/opt/Lattice_DMRG/gnu/*/bin/dmrg_config_user.sh" >> ~/.dmrg_config_user
echo -e "\e[32m DMRG installed and sourced \e[0m"

