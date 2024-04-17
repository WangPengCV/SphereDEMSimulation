#pragma once
#include "SphereProperties.h"
#include "PlanewallProperties.h"
#include "FiberProperties.h"
#include "CylinderwallProperties.h"
#include <set>
#include <tuple>
#include <memory>
#include <map>

class PropertyTypeID
{
public:
    PropertyTypeID() {}
    PropertyTypeID(int category, int subType) : category(category), subType(subType) {}

    bool operator<(const PropertyTypeID &other) const
    {
        return std::tie(category, subType) < std::tie(other.category, other.subType);
    }
    int getCategory() const { return category; }
    int getSubType() const { return subType; }

private:
    int category; // For example, 1 for spheres, 2 for cylinders, etc.
    int subType;  // Subtype to distinguish different properties within the same category
};

/**
 * @brief Class for managing particle properties and types.
 */
class ParticlePropertyManager
{
public:
    // Method to add SphereProperties
    void addSphereProperties(const PropertyTypeID &id, std::shared_ptr<SphereProperties> &properties);
    // Method to add PlanewallProperties
    void addPlanewallProperties(const PropertyTypeID &id, std::shared_ptr<PlanewallProperties> &properties);
    // Method to add FiberProperties
    void addFiberProperties(const PropertyTypeID &id, std::shared_ptr<FiberProperties> &properties);
    // Method to add FiberProperties
    void addCylinerwallProperties(const PropertyTypeID &id, std::shared_ptr<CylinderwallProperties> &properties);
    // Method to get SphereProperties
    const std::shared_ptr<SphereProperties> getSphereProperties(const PropertyTypeID &id) const
    {
        auto it = propertiesMap.find(id);
        if (it != propertiesMap.end())
        {
            // Use dynamic_cast to ensure the correct type is returned
            const auto sphereProperties = std::static_pointer_cast<SphereProperties>(it->second);
            if (sphereProperties)
            {
                return sphereProperties;
            }
            else
            {
                throw std::runtime_error("Property with ID is not of type SphereProperties.");
            }
        }
        // Property not found
        throw std::runtime_error("Property with ID not found.");
    }

    const std::shared_ptr<PlanewallProperties> getPlanewallProperties(const PropertyTypeID &id) const
    {
        auto it = propertiesMap.find(id);
        if (it != propertiesMap.end())
        {
            // Use dynamic_cast to ensure the correct type is returned
            const auto planewallProperties = std::static_pointer_cast<PlanewallProperties>(it->second);
            if (planewallProperties)
            {
                return planewallProperties;
            }
            else
            {
                throw std::runtime_error("Property with ID is not of type PlanewallProperties.");
            }
        }
        // Property not found
        throw std::runtime_error("Property with ID not found.");
    }

    const std::shared_ptr<FiberProperties> getFiberProperties(const PropertyTypeID &id) const
    {
        auto it = propertiesMap.find(id);
        if (it != propertiesMap.end())
        {
            // Use dynamic_cast to ensure the correct type is returned
            const auto fiberProperties = std::static_pointer_cast<FiberProperties>(it->second);
            if (fiberProperties)
            {
                return fiberProperties;
            }
            else
            {
                throw std::runtime_error("Property with ID is not of type FiberProperties.");
            }
        }
        // Property not found
        throw std::runtime_error("Property with ID not found.");
    }

    const std::shared_ptr<CylinderwallProperties> getCylinderwallProperties(const PropertyTypeID &id) const
    {
        auto it = propertiesMap.find(id);
        if (it != propertiesMap.end())
        {
            // Use dynamic_cast to ensure the correct type is returned
            const auto cylinderwallProperties = std::static_pointer_cast<CylinderwallProperties>(it->second);
            if (cylinderwallProperties)
            {
                return cylinderwallProperties;
            }
            else
            {
                throw std::runtime_error("Property with ID is not of type CylinderwallProperties.");
            }
        }
        // Property not found
        throw std::runtime_error("Property with ID not found.");
    }

    const std::map<PropertyTypeID, std::shared_ptr<ParticleProperties>> &getParticleProperties() const
    {
        return propertiesMap;
    }

private:
    std::map<PropertyTypeID, std::shared_ptr<ParticleProperties>> propertiesMap;
};
