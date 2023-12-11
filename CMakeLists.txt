cmake_minimum_required(VERSION 3.18)
project(lu)
 
set(CMAKE_CXX_STANDARD 11)
 
add_executable(main main.cpp)
add_executable(main1 main1.cpp)
add_executable(lu lu.cpp)
add_executable(cg1 cg1.cpp)
add_executable(test test.cpp)
 
# matplotlibcppで使うPythonとNumpyを探す
find_package(Python3 COMPONENTS Development NumPy)
target_include_directories(main PRIVATE ${Python3_INCLUDE_DIRS} ${Python3_NumPy_INCLUDE_DIRS})
target_link_libraries(main Python3::Python Python3::NumPy)
target_include_directories(main1 PRIVATE ${Python3_INCLUDE_DIRS} ${Python3_NumPy_INCLUDE_DIRS})
target_link_libraries(main1 Python3::Python Python3::NumPy)
target_include_directories(lu PRIVATE ${Python3_INCLUDE_DIRS} ${Python3_NumPy_INCLUDE_DIRS})
target_link_libraries(lu Python3::Python Python3::NumPy)
target_include_directories(cg1 PRIVATE ${Python3_INCLUDE_DIRS} ${Python3_NumPy_INCLUDE_DIRS})
target_link_libraries(cg1 Python3::Python Python3::NumPy)
target_include_directories(test PRIVATE ${Python3_INCLUDE_DIRS} ${Python3_NumPy_INCLUDE_DIRS})
target_link_libraries(test Python3::Python Python3::NumPy)
