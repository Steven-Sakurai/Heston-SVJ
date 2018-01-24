project(Heston_MC)
set(CMAKE_CXX_STANDARD 11)
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})
include_directories(.)
add_executable(a.out main.cpp)
target_link_libraries(a.out ${ARMADILLO_LIBRARIES})

