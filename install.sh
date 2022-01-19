cd /usr/local/include
sudo ln -sf eigen3/Eigen Eigen
sudo ln -sf eigen3/unsupported unsupported
cd /usr/include
sudo ln -sf eigen3/Eigen Eigen
sudo ln -sf eigen3/unsupported unsupported
sudo apt update
sudo apt install libeigen3-dev
sudo apt-get install libboost-all-dev