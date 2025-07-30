cmake_minimum_required(VERSION 3.16)
project(MyProject)

set(Kokkos_COMMON_SOURCE_DIR /home/amartya/software/kokkos-4.6.02/install)

find_package(Kokkos CONFIG)

include_directories(/home/amartya/Work/kokkos-trial)
add_executable(linkedcell_experiments linkedcell_experiments.cpp)
add_executable(linkedcell_experiments_parallel linkedcell_experiments_parallel.cpp)
add_executable(linkedcell_benchmarks linkedcell_benchmarks.cpp)
add_executable(linkedcell_tester linkedcell_tester.cpp)
add_executable(subview_playground subview_playground.cpp)
target_link_libraries(linkedcell_experiments Kokkos::kokkos)
target_link_libraries(linkedcell_experiments_parallel Kokkos::kokkos)
target_link_libraries(linkedcell_benchmarks Kokkos::kokkos)
target_link_libraries(linkedcell_tester Kokkos::kokkos)
target_link_libraries(subview_playground Kokkos::kokkos)
