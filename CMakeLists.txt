cmake_minimum_required(VERSION 3.1)

project(Snapshot)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Release)

find_package(CGAL)
find_package(Matlab)
find_package(Threads)

find_library(LIB_MATLAB_DATA_ARRAY 
	NAMES MatlabDataArray
	PATHS ${Matlab_ROOT_DIR}/extern/bin/glnxa64
	NO_DEFAULT_PATH)

find_library(LIB_MATLAB_ENGINE
	NAMES MatlabEngine
	PATHS ${Matlab_ROOT_DIR}/extern/bin/glnxa64
        NO_DEFAULT_PATH)

include(${CGAL_USE_FILE})
include_directories(${Matlab_INCLUDE_DIRS})

include_directories(.)
include_directories(eigen)

file(GLOB Sources
	src/main.cpp
    src/utility.hpp
    src/alpha_shape.hpp
    src/alpha_shape.cpp
    src/rips_shape.hpp
    src/rips_shape.cpp
    src/snapshot.hpp
    src/snapshot.cpp)

add_executable(Snapshot ${Sources})

target_link_libraries( Snapshot 
	${CGAL_LIBRARIES} 
	${CGAL_3RD_PARTY_LIBRARIES} 
	${LIB_MATLAB_DATA_ARRAY}
	${LIB_MATLAB_ENGINE}
	Threads::Threads 	
)	
