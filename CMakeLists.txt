project(UnscentedKF)
enable_testing()
include (CTest)

cmake_minimum_required (VERSION 3.5)

add_definitions(-std=c++11)

#setting env variables
set(CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS "${CXX_FLAGS}")

#I assume setting the source files from which all others will be built
#but what of the header files? Does it assume the header files are in the same directory?
set(sources src/ukf.cpp src/main.cpp src/tools.cpp)
set(test_sources src/ukf.cpp src/tools.cpp src/tests/ukf_test.cpp )
include_directories(src)

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin") 
#the next two seem to add directories to search for files under
include_directories(/usr/local/include)
include_directories(/usr/local/opt/openssl/include)
link_directories(/usr/local/lib)
link_directories(/usr/local/opt/openssl/lib)
link_directories(/usr/local/Cellar/libuv/1*/lib)

endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin") 

#names the destination executable of the make file
add_executable(UnscentedKF ${sources})
add_executable(TEST-1 ${test_sources})

target_link_libraries(UnscentedKF z ssl uv uWS)
target_link_libraries(TEST-1 z ssl uv uWS)

add_test(TEST-1 TEST-1)

#Questions: what directory does it start it's search from? Probably the one
#containing cmakelists.txt


