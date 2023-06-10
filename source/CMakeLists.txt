set(TARGET "step-1")

SET(TARGET_SRC step-1.cc)

# add_executable(step-1 step-1.cc)
# deal_ii_setup_target(step-1)
deal_ii_invoke_autopilot()
# target_link_libraries(pythagoras GTest::gtest_main)

# include(GoogleTest)
# gtest_discover_tests(pythagoras)