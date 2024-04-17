#pragma once
#include "BoundaryCondition.h"
#include "ParticlePropertyManager.h"
#include <unordered_map>
#include <Eigen/Dense>

class CylinderContainer : public BoundaryCondition
{
public:
    CylinderContainer(){};
    // Constructor
    CylinderContainer(int id, const PropertyTypeID &type, int state,
                      double radius,
                      const Eigen::Vector3d &startpoint,
                      const Eigen::Vector3d &endpoint,
                      const Eigen::Vector3d &velocity = Eigen::Vector3d::Zero(),
                      const Eigen::Vector3d &force = Eigen::Vector3d::Zero());

   
    void setMass(std::shared_ptr<ParticlePropertyManager> manager);

    void setVelocity(const Eigen::Vector3d &Velocity);

    void setStartpoint(const Eigen::Vector3d &startpoint);

    void setEndpoint(const Eigen::Vector3d &endpoint);

    void addForce(const Eigen::Vector3d &Force);

    void resetForce();


    // Override the save_tostring method from the BoundaryCondition class
    virtual std::string save_tostring() const override;

   
    double getMass() const { return mass; }
    double getRadius() const { return radius;}
    const Eigen::Vector3d &getVelocity() const { return velocity; }
    const Eigen::Vector3d &getForce() const { return force; }
    const Eigen::Vector3d &getStartpoint() const { return startpoint; }
    const Eigen::Vector3d &getEndpoint() const { return endpoint;}
    void  generateMesh(double meshResolution);
    const std::vector<Eigen::Vector3d> &getMeshVertices() const{ return meshVertices;};

private:
    Eigen::Vector3d startpoint;
    Eigen::Vector3d endpoint;
    double radius;
    double mass;
    Eigen::Vector3d force;
    Eigen::Vector3d velocity;
    std::vector<Eigen::Vector3d> meshVertices; // discreted points for contact detection

};
