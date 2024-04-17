#include "ParticlePropertyManager.h"

void ParticlePropertyManager::addSphereProperties(const PropertyTypeID& id, std::shared_ptr<SphereProperties>& properties)
{
    // Check if the property already exists
    if (propertiesMap.find(id) != propertiesMap.end()) {
        throw std::runtime_error("Property with the same ID already exists.");
    }    
    // Add the property if no conflict
    propertiesMap[id] = properties;
}

void ParticlePropertyManager::addPlanewallProperties(const PropertyTypeID& id, std::shared_ptr<PlanewallProperties>& properties)
{
     // Check if the property already exists
    if (propertiesMap.find(id) != propertiesMap.end()) {
        throw std::runtime_error("Property with the same ID already exists.");
    } 
    propertiesMap[id] = properties;
}

void ParticlePropertyManager::addFiberProperties(const PropertyTypeID &id, std::shared_ptr<FiberProperties>& properties)
{

     // Check if the property already exists
    if (propertiesMap.find(id) != propertiesMap.end()) {
        throw std::runtime_error("Property with the same ID already exists.");
    } 
    propertiesMap[id] = properties;
}

void ParticlePropertyManager::addCylinerwallProperties(const PropertyTypeID &id, std::shared_ptr<CylinderwallProperties> &properties)
{
    if (propertiesMap.find(id) != propertiesMap.end()) {
        throw std::runtime_error("Property with the same ID already exists.");
    } 
    propertiesMap[id] = properties;
}