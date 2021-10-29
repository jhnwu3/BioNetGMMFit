# **Snapshot**
Snapshot is a C++ software designed for rate constant estimation of CYTOF Snapshot Data. 
It takes into account both linear and nonlinear models of evolution for estimating rate constants.

## **3rd Party Libraries**

### *Compatibility*
All software currently has only been tested on Linux-based systems.

### *Eigen*
Snapshot uses the Eigen 3.3.9 linear algebra library for matrix computations. If on Ubuntu, you can do a quick install using:

    sudo apt install libeigen3-dev

otherwise you can see more detailed install instructions [here](https://eigen.tuxfamily.org/dox/GettingStarted.html)

### *Boost*
Snapshot uses the Boost 1.7.2 odeint C++ library for ODE estimations for nonlinear systems. To install the whole boost C++ library, you can try:

    sudo apt-get install libboost-all-dev

However, Snapshot only uses the C++ odeint library, so if storage space is an explicit concern, more
detailed install intructions can be found [here](https://www.boost.org/doc/libs/1_77_0/more/getting_started/unix-variants.html)

## **Getting Started** ##











