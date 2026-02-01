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
add_library(PROPOSAL_cxxwrap SHARED ${WRAPPER_SOURCES})

target_link_libraries(PROPOSAL_cxxwrap
    PRIVATE
    JlCxx::cxxwrap_julia
    PROPOSAL::PROPOSAL
)

# On macOS, avoid directly linking libjulia to prevent dlopen deadlock.
# Julia symbols will be resolved at runtime since the host process already has them.
if(APPLE)
    target_link_options(PROPOSAL_cxxwrap PRIVATE "LINKER:-undefined,dynamic_lookup")
    # Replace libjulia in JlCxx interface with nothing - symbols resolve at runtime
    set_target_properties(JlCxx::cxxwrap_julia PROPERTIES
        INTERFACE_LINK_LIBRARIES "-framework CoreFoundation"
    )
endif()

# Install
install(TARGETS PROPOSAL_cxxwrap
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
)
