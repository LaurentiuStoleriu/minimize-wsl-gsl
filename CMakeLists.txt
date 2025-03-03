cmake_minimum_required(VERSION 3.0.0)
project(minimize-wsl-gsl VERSION 0.1.0)

include(CTest)
enable_testing()

add_executable(minimize-wsl-gsl main.cpp)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

# Packages
find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
target_link_libraries(minimize-wsl-gsl ${GSL_LIBRARIES})