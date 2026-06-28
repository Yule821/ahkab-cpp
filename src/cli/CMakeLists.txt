set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(ahkab_cli main.cpp)
target_link_libraries(ahkab_cli PRIVATE ahkab)
