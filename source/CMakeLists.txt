set(TARGET "step-1")

SET(TARGET_SRC step-1.cc gtest.cc)

deal_ii_invoke_autopilot()

# add_executable(gtest gtest.cc)
target_link_libraries(gtest GTest::gtest_main)

include(GoogleTest)
gtest_discover_tests(gtest)