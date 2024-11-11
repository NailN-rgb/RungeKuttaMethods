cmake_minimum_required(VERSION 3.0)

project(Runge-Cuttas)

include_directories(${PROJECT_SOURCE_DIR}/..)

find_package(PythonLibs REQUIRED)
find_package(Python2 REQUIRED)

include_directories(${PYTHON_INCLUDE_DIRS})

add_executable(Runge-Cuttas main.cpp)

target_link_libraries(Runge-Cuttas PRIVATE -lpython3.5m)
