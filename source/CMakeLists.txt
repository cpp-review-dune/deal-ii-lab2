# set(TARGET "step-1")
# SET(TARGET_SRC step-1.cc)

# add_executable(step-1 step-1.cc)
# deal_ii_setup_target(step-1)

# add_executable(step-2 step-2.cc)
# deal_ii_setup_target(step-2)

add_executable(gtest gtest.cc)
deal_ii_setup_target(gtest)
target_link_libraries(gtest GTest::gtest_main)

include(GoogleTest)
gtest_discover_tests(gtest)
# deal_ii_invoke_autopilot()