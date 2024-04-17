# Define the ParticleProperties library
add_library(DynamicBoundaryMover_lib
        DynamicBoundaryMover.cpp
        SpringOscillator.cpp
        ConstantMove.cpp
)

# Link with ParticleProperties_lib to use PropertyTypeID
target_link_libraries(DynamicBoundaryMover_lib
    Particle_lib
    ParticleProperties_lib
    BoundaryCondition_lib
)
# Specify the include directory for the ParticleProperties library
target_include_directories(DynamicBoundaryMover_lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)
