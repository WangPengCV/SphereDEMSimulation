# Define the ParticleProperties library
add_library(BoundaryCondition_lib
    BoundaryCondition.cpp
    PlaneWall.cpp
    #RectangularContainer.cpp
    CylinderContainer.cpp
    ConnectedGeometry.cpp
)

# Link with ParticleProperties_lib to use PropertyTypeID
target_link_libraries(BoundaryCondition_lib
    Particle_lib
    ParticleProperties_lib
)
# Specify the include directory for the ParticleProperties library
target_include_directories(BoundaryCondition_lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)
