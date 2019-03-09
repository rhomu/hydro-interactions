# This file is a part of my soul, Copyright (C) 2018, Romain Mueller
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

cmake_minimum_required(VERSION 2.8)
project(hydro)
set(CMAKE_CXX_STANDARD 11)

# default is release
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

# set warning level
if(MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Wall")
elseif(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long")
endif()

################################################################################
# directories
################################################################################

# directories
include_directories(src)
add_subdirectory(src)

################################################################################
# target definitions
################################################################################

# use all files in src/
file(GLOB_RECURSE sources src/*.cpp src/*.hpp)

# main executable
add_executable(hydro ${sources})

################################################################################
# dependencies
################################################################################

# check for boost libraries
find_package(Boost 1.36.0 COMPONENTS program_options REQUIRED)
target_link_libraries(hydro PUBLIC ${Boost_LIBRARIES})

# eigen we need!
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries (hydro PUBLIC Eigen3::Eigen)

# add hydra support
option(HYDRA "Make it work on HYDRA (as of 2018)" OFF)
if(HYDRA)
  message("-- Setting HYDRA paths")
  set_target_properties(hydro PROPERTIES LINK_FLAGS "${CMAKE_LINK_FLAGS} -L/usr/local/shared/boost/1.64.0-gcc5.4.0/lib/")
  target_include_directories(hydro PRIVATE "/usr/local/shared/eigen3/include/eigen3/")
endif()

################################################################################
# testing
################################################################################

################################################################################
# packaging
################################################################################
