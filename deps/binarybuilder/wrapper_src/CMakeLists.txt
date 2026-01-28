cmake_minimum_required(VERSION 3.15)
project(PROPOSAL_jl LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find dependencies
find_package(JlCxx REQUIRED)
find_package(PROPOSAL REQUIRED)

# Wrapper sources
set(WRAPPER_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/PROPOSAL_jl.cxx
)

# Create the wrapper library
add_library(PROPOSAL_jl SHARED ${WRAPPER_SOURCES})

target_link_libraries(PROPOSAL_jl
    PRIVATE
    JlCxx::cxxwrap_julia
    PROPOSAL::PROPOSAL
)

# Install
install(TARGETS PROPOSAL_jl
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
)
