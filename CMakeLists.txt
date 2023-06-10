cmake_minimum_required(VERSION 3.10)
project("Lab 0")

find_package(deal.II 9.4.2 REQUIRED)
deal_ii_initialize_cached_variables()

find_package(GTest)
enable_testing()

add_subdirectory(source)