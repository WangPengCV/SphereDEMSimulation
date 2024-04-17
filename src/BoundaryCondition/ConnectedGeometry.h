#pragma once
#include "ParticlePropertyManager.h"
#include <vector>
#include <tuple>
#include <Eigen/Dense>
#include <bitset>

class ConnectedGeometry
{
public:
    enum class DegreeOfFreedom
    {
        X = 0,
        Y,
        Z
    };
    ConnectedGeometry() = default;
    ConnectedGeometry(int id, Eigen::Vector3d &center, Eigen::Vector3d &velocity, std::bitset<3> dofs = std::bitset<3>(7));
    
    void setDegreeOfFreedom(DegreeOfFreedom dof, bool value);

    bool getDegreeOfFreedom(DegreeOfFreedom dof) const
    {
        return dofs[static_cast<size_t>(dof)];
    }
    void addGeometry(const std::string &geoType, int geometryId);

    void setGeometry(const std::vector<std::tuple<std::string, int>> &geometrys);

    void setVelocity(const Eigen::Vector3d &velocity);

    const Eigen::Vector3d &getVelocity() const { return velocity; }

    void setMass(double mass);

    double getMass() const { return mass; }

    void setCenterOfMass(const Eigen::Vector3d &centerOfMass);

    void addForce(const Eigen::Vector3d &force);

    int GeometrySize() const;

    void resetForce();

    std::string save_tostring() const;

    const std::vector<std::tuple<std::string, int>> &getGeometrys() const { return geometrys; }

    const Eigen::Vector3d &getCenterOfMass() const { return centerOfMass; }

private:
    int id;
    double mass;
    std::vector<std::tuple<std::string, int>> geometrys;
    Eigen::Vector3d centerOfMass;
    Eigen::Vector3d velocity;
    Eigen::Vector3d force;
    std::bitset<3> dofs; 
};
