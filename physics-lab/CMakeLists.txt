cmake_minimum_required(VERSION 3.16)
project(Graviton
    VERSION 1.0.0
    DESCRIPTION "High-Performance Symplectic Physics Engine & Spatial Hash Grid"
    LANGUAGES CXX
)

# -----------------------------------------------------------------------------
# C++ Standard & Compiler Configurations
# -----------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Enable high compiler warning levels & aggressive optimization
if(MSVC)
    add_compile_options(/W4 /O2 /permissive-)
else()
    add_compile_options(-Wall -Wextra -Wpedantic -O3)
endif()

# -----------------------------------------------------------------------------
# Library & Executable Targets
# -----------------------------------------------------------------------------
add_executable(graviton_benchmark
    src/main.cpp
)

target_include_directories(graviton_benchmark PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# Set output directory to binary root
set_target_properties(graviton_benchmark PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)
