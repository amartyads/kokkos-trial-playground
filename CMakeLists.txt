cmake_minimum_required(VERSION 3.16)
project(MyProject)

set(Kokkos_COMMON_SOURCE_DIR /home/amartya/software/kokkos-4.6.02/install)

find_package(Kokkos CONFIG)

add_executable(linkedcell_experiments linkedcell_experiments.cpp)
target_link_libraries(linkedcell_experiments Kokkos::kokkos)

