cmake_minimum_required(VERSION 3.16)

add_executable(test_circuit
  test_circuit.cpp
)
target_link_libraries(test_circuit PRIVATE ahkab GTest::gtest GTest::gtest_main)
add_test(NAME test_circuit COMMAND test_circuit)

add_executable(test_op
  test_op.cpp
)
target_link_libraries(test_op PRIVATE ahkab GTest::gtest GTest::gtest_main)
add_test(NAME test_op COMMAND test_op)

add_executable(test_netlist
  test_netlist.cpp
)
target_link_libraries(test_netlist PRIVATE ahkab GTest::gtest GTest::gtest_main)
add_test(NAME test_netlist COMMAND test_netlist)

add_executable(test_integration
  test_integration.cpp
)
target_link_libraries(test_integration PRIVATE ahkab GTest::gtest GTest::gtest_main)
add_test(NAME test_integration COMMAND test_integration)
