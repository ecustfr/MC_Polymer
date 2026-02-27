cmake_minimum_required(VERSION 3.16)
project(TestBulkDensity LANGUAGES CXX)

# 设置C++标准
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 包含目录
include_directories(${CMAKE_SOURCE_DIR}/../include)

# 源文件
set(SOURCES
    ${CMAKE_SOURCE_DIR}/../src/Bulk_RingPolymer.cpp
    ${CMAKE_SOURCE_DIR}/../src/Utils.cpp
    ${CMAKE_SOURCE_DIR}/../test_bulk_density_final.cpp
)

# 添加可执行文件
add_executable(test_bulk_density ${SOURCES})

# MSVC特定设置
if(MSVC)
    target_compile_options(test_bulk_density PRIVATE
        /W3
        /utf-8
    )
endif()