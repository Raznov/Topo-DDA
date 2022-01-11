# Install script for directory: E:/Topo-DDA-forWin64/eigen-eigen-94875feeeeb9/unsupported/Eigen

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "E:/Topo-DDA-forWin64/eigen-eigen-94875feeeeb9/out/install/x64-Debug (默认值)")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Debug")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xDevelx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/unsupported/Eigen" TYPE FILE FILES
    "E:/Topo-DDA-forWin64/eigen-eigen-94875feeeeb9/unsupported/Eigen/AdolcForward"
    "E:/Topo-DDA-forWin64/eigen-eigen-94875feeeeb9/unsupported/Eigen/AlignedVector3"
    "E:/Topo-DDA-forWin64/eigen-eigen-94875feeeeb9/unsupported/Eigen/ArpackSupport"
    "E:/Topo-DDA-forWin64/eigen-eigen-94875feeeeb9/unsupported/Eigen/AutoDiff"
    "E:/Topo-DDA-forWin64/eigen-eigen-94875feeeeb9/unsupported/Eigen/BVH"
    "E:/Topo-DDA-forWin64/eigen-eigen-94875feeeeb9/unsupported/Eigen/EulerAngles"
    "E:/Topo-DDA-forWin64/eigen-eigen-94875feeeeb9/unsupported/Eigen/FFT"
    "E:/Topo-DDA-forWin64/eigen-eigen-94875feeeeb9/unsupported/Eigen/IterativeSolvers"
    "E:/Topo-DDA-forWin64/eigen-eigen-94875feeeeb9/unsupported/Eigen/KroneckerProduct"
    "E:/Topo-DDA-forWin64/eigen-eigen-94875feeeeb9/unsupported/Eigen/LevenbergMarquardt"
    "E:/Topo-DDA-forWin64/eigen-eigen-94875feeeeb9/unsupported/Eigen/MatrixFunctions"
    "E:/Topo-DDA-forWin64/eigen-eigen-94875feeeeb9/unsupported/Eigen/MoreVectorization"
    "E:/Topo-DDA-forWin64/eigen-eigen-94875feeeeb9/unsupported/Eigen/MPRealSupport"
    "E:/Topo-DDA-forWin64/eigen-eigen-94875feeeeb9/unsupported/Eigen/NonLinearOptimization"
    "E:/Topo-DDA-forWin64/eigen-eigen-94875feeeeb9/unsupported/Eigen/NumericalDiff"
    "E:/Topo-DDA-forWin64/eigen-eigen-94875feeeeb9/unsupported/Eigen/OpenGLSupport"
    "E:/Topo-DDA-forWin64/eigen-eigen-94875feeeeb9/unsupported/Eigen/Polynomials"
    "E:/Topo-DDA-forWin64/eigen-eigen-94875feeeeb9/unsupported/Eigen/Skyline"
    "E:/Topo-DDA-forWin64/eigen-eigen-94875feeeeb9/unsupported/Eigen/SparseExtra"
    "E:/Topo-DDA-forWin64/eigen-eigen-94875feeeeb9/unsupported/Eigen/SpecialFunctions"
    "E:/Topo-DDA-forWin64/eigen-eigen-94875feeeeb9/unsupported/Eigen/Splines"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xDevelx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/unsupported/Eigen" TYPE DIRECTORY FILES "E:/Topo-DDA-forWin64/eigen-eigen-94875feeeeb9/unsupported/Eigen/src" FILES_MATCHING REGEX "/[^/]*\\.h$")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("E:/Topo-DDA-forWin64/eigen-eigen-94875feeeeb9/out/build/x64-Debug (默认值)/unsupported/Eigen/CXX11/cmake_install.cmake")

endif()

