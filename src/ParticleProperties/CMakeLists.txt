# Define the ParticleProperties library
add_library(ParticleProperties_lib
    ParticleProperties.cpp
    SphereProperties.cpp
    PlanewallProperties.cpp
    ParticlePropertyManager.cpp
    FiberProperties.cpp
    CylinderwallProperties.cpp
)

# Specify the include directory for the ParticleProperties library
target_include_directories(ParticleProperties_lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)
