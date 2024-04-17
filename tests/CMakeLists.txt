add_executable(SphereParticleTest SphereParticleTest.cpp)
target_link_libraries(SphereParticleTest Particle_lib DEM_lib)

add_executable( ParticlePropertiesTest ParticlePropertiesTest.cpp)
target_link_libraries(ParticlePropertiesTest ParticleProperties_lib DEM_lib)



add_executable(VisualizeGridBasedContactDetection VisualizeGridBasedContactDetection.cpp)
target_link_libraries(VisualizeGridBasedContactDetection DEM_lib Visualization_lib ${VTK_LIBRARIES})

vtk_module_autoinit(
    TARGETS VisualizeGridBasedContactDetection
    MODULES ${VTK_LIBRARIES}
  )

add_executable(DEM DEM.cpp)
target_link_libraries(DEM DEM_lib Visualization_lib ${VTK_LIBRARIES})



vtk_module_autoinit(
    TARGETS DEM
    MODULES ${VTK_LIBRARIES}
  )

add_executable(PostProcess PostProcess.cpp)
target_link_libraries(PostProcess DEM_lib Visualization_lib ${VTK_LIBRARIES})
  

add_executable(ParallelTest ParallelTest.cpp)
target_link_libraries(ParallelTest DEM_lib Visualization_lib ${VTK_LIBRARIES})



vtk_module_autoinit(
    TARGETS ParallelTest
    MODULES ${VTK_LIBRARIES}
  )
add_executable(DamainDecomposition DamainDecomposition.cpp)
target_link_libraries(DamainDecomposition DEM_lib Visualization_lib ${VTK_LIBRARIES})



vtk_module_autoinit(
    TARGETS DamainDecomposition
    MODULES ${VTK_LIBRARIES}
   )