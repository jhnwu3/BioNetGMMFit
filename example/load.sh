# loads necessary modules on NCH cluster, use source load.sh
git pull
module load Eigen GCC/9.3.0 OpenMPI/4.0.3 Boost/1.72.0 gnuplot
make
