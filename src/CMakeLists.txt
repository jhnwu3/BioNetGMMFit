cmake_minimum_required(VERSION 3.20)

project(CyGMM)


set(CMAKE_CXX_STANDARD 17)


# could be command line variables (-DROADRUNNER_INSTALL_PREFIX=...)

set(ROADRUNNER_INSTALL_PREFIX "../../buildroadrunner/roadrunner/install-Release/")

set(LLVM_INSTALL_PREFIX "../../buildroadrunner/llvm13-ubuntu-gcc10-rel")

set(CMAKE_PREFIX_PATH

        "${ROADRUNNER_INSTALL_PREFIX}/cmake" # needed because this is where poco installs config files to

        "${ROADRUNNER_INSTALL_PREFIX}/lib/cmake" # all other normal packages install config files to here

        )

# for registering location of ImportRoadRunnerAndDependencies
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../) # move executable to main directory
set(CMAKE_MODULE_PATH "${ROADRUNNER_INSTALL_PREFIX}/lib/cmake")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17  -static-libgcc -static-libstdc++ -lstdc++fs -O3")
# Import roadrunner and all necessary dependencies for linking roadrunner.

# This command is a roadrunner convenience cmake script which does all of the

# calls to find_package(xxx, CONFIG) for us. If we didn't include this

# script, we would manually need to call find_package for roadrunner-static

# as well as roadrunner dependencies. Checkout The script

# under <rr_install_prefix>/lib/cmake.

include(ImportRoadrunnerAndDependencies)


# add an executable
find_package(OpenMP) # openMP for parallelization
add_executable(${PROJECT_NAME} main.cpp main.hpp calc.cpp calc.hpp fileIO.cpp fileIO.hpp linear.cpp linear.hpp nonlinear.cpp nonlinear.hpp system.hpp system.cpp sbml.cpp sbml.hpp param.hpp cli.hpp cli.cpp)
target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)

# link roadrunner-static. Note that we have configured roadrunner-static target (which is imported

# under the alias roadrunner-static::roadrunner-static) with PUBLIC visibility, which means

# all includes and transitive linked libraries are automatically handled by cmake.

# without this, you would have to manually link roadrunner dependencies and manually

# include roadrunner include directories with something like:

#   - target_include_directories(UseRoadRunnerFromCxx PRIVATE "${ROADRUNNER_INSTALL_PREFIX}/include")
target_link_libraries(${PROJECT_NAME} PRIVATE roadrunner-static::roadrunner-static stdc++fs)
if(OpenMP_CXX_FOUND)
    target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
endif()
# target_include_directories(UseRoadRunnerFromCxx PRIVATE "${ROADRUNNER_INSTALL_PREFIX}/include")