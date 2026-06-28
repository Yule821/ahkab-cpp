add_library(ahkab STATIC
  core/circuit.cpp
  devices/linear_devices.cpp
  devices/nonlinear_devices.cpp
  mna/mna_builder.cpp
  mna/nr_solver.cpp
  mna/linear_solver.cpp
  analysis/op_analysis.cpp
  analysis/dc_analysis.cpp
  analysis/transient_analysis.cpp
  analysis/ac_analysis.cpp
  analysis/pz_analysis.cpp
  parser/lexer.cpp
  parser/parser.cpp
  output/csv_writer.cpp
  output/table_writer.cpp
)

target_include_directories(ahkab PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(ahkab PUBLIC
  Eigen3::Eigen
)

target_compile_options(ahkab PRIVATE
  $<$<CXX_COMPILER_ID:GNU>:-Wall -Wextra -Wpedantic>
  $<$<CXX_COMPILER_ID:Clang>:-Wall -Wextra -Wpedantic>
)
