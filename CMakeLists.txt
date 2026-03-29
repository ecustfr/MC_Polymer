cmake_minimum_required(VERSION 3.16)

project(MuVTPolymer VERSION 1.0 LANGUAGES CXX)



# 关键设置：指定源代码编码为 UTF-8


# 对于 pybind11 项目特别重要
add_compile_definitions(_SILENCE_CXX17_CODECVT_HEADER_DEPRECATION_WARNING)
# 直接指定 conda 环境中的 Python
set(CONDA_ENV_PATH "D:/myanaconda3")
#set(CONDA_ENV_PATH "T:/anaconda/envs/mcpolymer")
set(Python3_EXECUTABLE "${CONDA_ENV_PATH}/python.exe")
message(STATUS "Using Python from: ${Python3_EXECUTABLE}")

find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
# 设置C++标准
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#set(python_DIR "T:/anaconda/envs/mcpolymer")

set(pybind11_DIR "D:/myanaconda3/Lib/site-packages/pybind11/share/cmake/pybind11")
# set(pybind11_DIR "T:/anaconda/envs/mcpolymer/Lib/site-packages/pybind11/share/cmake/pybind11")
find_package(pybind11 CONFIG REQUIRED)

# 设置包含目录
include_directories(${PROJECT_SOURCE_DIR}/include)

# 源文件列表（只包含必要的文件）
set(SOURCES
    ${PROJECT_SOURCE_DIR}/src/MuVT_MC_LinearPolymer.cpp
    ${PROJECT_SOURCE_DIR}/src/MuVT_MC_RingPolymer.cpp
    ${PROJECT_SOURCE_DIR}/src/Utils.cpp
)


add_library(Utils STATIC src/Utils.cpp)


add_library(MuVTPolymer STATIC
 src/MuVT_MC_LinearPolymer.cpp 
 src/MuVT_MC_RingPolymer.cpp)

target_link_libraries(MuVTPolymer PRIVATE Utils)


# === 添加 Python 模块 ===
pybind11_add_module(pymcpolymer ${PROJECT_SOURCE_DIR}/src/pybind11_module.cpp)
target_link_libraries(pymcpolymer PRIVATE MuVTPolymer Python3::Python)


# === 关键修改 ===
# MSVC 不需要链接 m 库，只有非 MSVC (如 Linux/MinGW) 才需要
if(NOT MSVC)
    target_link_libraries(MuVTPolymer PRIVATE m)
    target_link_libraries(pymcpolymer PRIVATE m)
    target_link_libraries(Utils PRIVATE m)
    target_link_libraries(BulkRingPolymer PRIVATE m)
    target_link_libraries(BulkLinearPolymer PRIVATE m)
endif()

# === 针对 Visual Studio 的优化建议 ===
if(MSVC)
    # 启用多核编译 (加快编译速度)，开启常用警告，强制源文件和执行字符集为 UTF-8，解决乱码
    add_compile_options(/O2 /Oi /Ot /GL)
    # 链接器优化：启用引用消除和中间代码生成
    add_link_options(/LTCG /OPT:REF /OPT:ICF)
    target_compile_options(MuVTPolymer PRIVATE
        /MP
        /W3
        /utf-8
    )
    target_compile_options(pymcpolymer PRIVATE
        /MP
        /W3
        /utf-8
    )
    target_compile_options(Utils PRIVATE
        /MP
        /W3
        /utf-8
    )
    add_compile_options(/utf-8)
    # 设置执行字符集为 UTF-8（可选，/utf-8 已包含此功能）
    add_compile_definitions(_UTF8_SOURCE)
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
