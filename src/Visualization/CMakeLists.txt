# Define the ParticleProperties library
add_library(Visualization_lib
   Visualization.cpp
   LineSphereSource.cpp
   LineCylinderSource.cpp
)

# Link with ParticleProperties_lib to use PropertyTypeID
target_link_libraries(Visualization_lib
    Particle_lib
    ParticleProperties_lib
    BoundaryCondition_lib
    DEM_lib
    ${VTK_LIBRARIES}
)
# Specify the include directory for the ParticleProperties library
target_include_directories(Visualization_lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)
vtk_module_autoinit(
    TARGETS Visualization_lib
    MODULES ${VTK_LIBRARIES}
)