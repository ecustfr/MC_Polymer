cmake_minimum_required(VERSION 3.16)

project(MuVTPolymer VERSION 1.0 LANGUAGES CXX)


# 设置C++标准
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Python 3 COMPONENTS Interpreter Development REQUIRED)

set(pybind11_DIR "D:/myanaconda3/Lib/site-packages/pybind11/share/cmake/pybind11")
find_package(pybind11 CONFIG REQUIRED)

# 设置包含目录
include_directories(${PROJECT_SOURCE_DIR}/include)

# 源文件列表（只包含必要的文件）
set(SOURCES
    ${PROJECT_SOURCE_DIR}/src/MuVT_MC_LinearPolymer.cpp
    ${PROJECT_SOURCE_DIR}/src/MuVT_MC_RingPolymer.cpp
    ${PROJECT_SOURCE_DIR}/src/Utils.cpp
)

add_library(MuVTPolymer src/MuVT_MC_LinearPolymer.cpp src/Utils.cpp src/MuVT_MC_RingPolymer)
# === 添加 Python 模块 ===
pybind11_add_module(pymcpolymer ${PROJECT_SOURCE_DIR}/src/pybind11_module.cpp)
target_link_libraries(pymcpolymer PRIVATE MuVTPolymer)

# === 关键修改 ===
# MSVC 不需要链接 m 库，只有非 MSVC (如 Linux/MinGW) 才需要
if(NOT MSVC)
    target_link_libraries(MuVTPolymer PRIVATE m)
    target_link_libraries(pymcpolymer PRIVATE m)
endif()

# === 针对 Visual Studio 的优化建议 ===
if(MSVC)
    target_compile_options(MuVTPolymer PRIVATE 
        /MP   # 启用多核编译 (加快编译速度)
        /W3   # 开启常用警告
        /utf-8      # 强制源文件和执行字符集为 UTF-8，解决乱码
    )
    target_compile_options(pymcpolymer PRIVATE 
        /MP   # 启用多核编译 (加快编译速度)
        /W3   # 开启常用警告
        /utf-8      # 强制源文件和执行字符集为 UTF-8，解决乱码
    )
endif()

add_custom_command(TARGET pymcpolymer POST_BUILD  COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:pymcpolymer> ${PROJECT_SOURCE_DIR})

# 打印配置信息
message(STATUS "=========================================")
message(STATUS "Project: ${PROJECT_NAME}")
message(STATUS "Compiler: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "pybind11 found: ${pybind11_FOUND}")
message(STATUS "Source files included:")
foreach(source ${SOURCES})
    message(STATUS "  - ${source}")
endforeach()
message(STATUS "=========================================")
