#pragma once
#include <memory>
#include <Eigen/Dense>
#include "SphereParticle.h"
#include "SphereCylinderBond.h"
#include "PlanewallProperties.h"
#include "CylinderwallProperties.h"
#include "PlaneWall.h"
#include "CylinderContainer.h"
#include "ParticlePropertyManager.h"
#include <unordered_map>

class ContactInformation
{
public:
    int particleId1;                            // ID of the first particle in contact
    int particleId2;                            // ID of the second particle in contact
    Eigen::Vector3d tangentialDisplacement;     // Cumulative tangential displacement
    Eigen::Vector3d previousTangentialVelocity; // Tangential velocity in the previous timestep
    Eigen::Vector3d normalForce;                // Computed normal force during the last interaction
    Eigen::Vector3d tangentialForce;            // Computed tangential force during the last interaction

    ContactInformation(int id1, int id2) : particleId1(id1), particleId2(id2),
                                           tangentialDisplacement(Eigen::Vector3d::Zero()),
                                           previousTangentialVelocity(Eigen::Vector3d::Zero()),
                                           normalForce(Eigen::Vector3d::Zero()),
                                           tangentialForce(Eigen::Vector3d::Zero()) {}
    ContactInformation() : particleId1(-1), particleId2(-1), // -1 or another invalid ID
                           tangentialDisplacement(Eigen::Vector3d::Zero()),
                           previousTangentialVelocity(Eigen::Vector3d::Zero()),
                           normalForce(Eigen::Vector3d::Zero()),
                           tangentialForce(Eigen::Vector3d::Zero())
    {
    }
};

class ContactForce
{
public:
    ContactForce();
    void addParticleProperties(const std::shared_ptr<ParticlePropertyManager> &manager);

    void computerSphereSphereEffective(std::shared_ptr<SphereProperties> &sphereproperties1, PropertyTypeID &type1,
                                       std::shared_ptr<SphereProperties> &sphereproperties2, PropertyTypeID &type2);

    void computerSpherePlanewallEffective(std::shared_ptr<PlanewallProperties> &planewallproperties, PropertyTypeID &type1,
                                          std::shared_ptr<SphereProperties> &sphereproperties, PropertyTypeID &type2);

    void computerSphereCylinderwallEffective(std::shared_ptr<CylinderwallProperties> &cylinderwallproperties, PropertyTypeID &type1,
                                             std::shared_ptr<SphereProperties> &sphereproperties, PropertyTypeID &type2);

    void computerFiberFiberEffective(std::shared_ptr<FiberProperties> &fibereproperties1, PropertyTypeID &type1,
                                     std::shared_ptr<FiberProperties> &fiberproperties2, PropertyTypeID &type2);

    void computerSphereFiberEffective(std::shared_ptr<SphereProperties> &sphereproperties1, PropertyTypeID &type1,
                                      std::shared_ptr<FiberProperties> &fiberproperties2, PropertyTypeID &type2);

    void computerFiberPlanewallEffective(std::shared_ptr<PlanewallProperties> &planewallproperties, PropertyTypeID &type1,
                                         std::shared_ptr<FiberProperties> &fiberproperties2, PropertyTypeID &type2);

    void computerFiberCylinderwallEffective(std::shared_ptr<CylinderwallProperties> &cylinderwallproperties, PropertyTypeID &type1,
                                            std::shared_ptr<FiberProperties> &fiberproperties2, PropertyTypeID &type2);
    // Method to compute force between two spherical particles
    void computeSphereSphereForce(const std::unique_ptr<SphereParticle> &sphere1, const std::unique_ptr<SphereParticle> &sphere2, double timeStep);
    // Method to compute force between a sphere particle and a plane wall
    void computePlaneWallSphereForce(const std::unique_ptr<PlaneWall> &planeWall, const std::unique_ptr<SphereParticle> &sphere, double timeStep);

    void computeCylinderWallSphereForce(const std::unique_ptr<CylinderContainer> &cylinderWall, const std::unique_ptr<SphereParticle> &sphere, double timeStep);

    // Method to compute force between two spherical particles
    void computeSphereFiberForce(const std::unique_ptr<SphereParticle> &sphere, const std::unique_ptr<SphereCylinderBond> &fiberbond,
                                 const std::unique_ptr<SphereParticle> &node1, const std::unique_ptr<SphereParticle> &node2, double timeStep);
    // Method to compute force between a fiber particle and a plane wall
    void computePlaneWallFiberForce(const std::unique_ptr<PlaneWall> &planeWall, const std::unique_ptr<SphereCylinderBond> &fiberbond,
                                    const std::unique_ptr<SphereParticle> &node1, const std::unique_ptr<SphereParticle> &node2, double timeStep);

    void computeCylinderWallFiberForce(const std::unique_ptr<CylinderContainer> &cylinderWall, const std::unique_ptr<SphereCylinderBond> &fiberbond,
                                       const std::unique_ptr<SphereParticle> &node1, const std::unique_ptr<SphereParticle> &node2, double timeStep);
    // Method to compute force between two spherical particles
    void computeFiberFiberForce(const std::unique_ptr<SphereCylinderBond> &fiberbond1, const std::unique_ptr<SphereParticle> &bondonenode1, const std::unique_ptr<SphereParticle> &bondonenode2,
                                const std::unique_ptr<SphereCylinderBond> &fiberbond2, const std::unique_ptr<SphereParticle> &bondtwonode1, const std::unique_ptr<SphereParticle> &bondtwonode2,
                                double timeStep);
    void updateContactInformation();

    double getEnergydissipation() const { return energyDissipation;}
    void setEnergydissipation(double energyDissipation);


    std::unordered_map<int, std::unordered_map<int, ContactInformation>> &getSphereSphereContactInformationList();
    std::unordered_map<int, std::unordered_map<int, ContactInformation>> &getSphereFiberContactInformationList();

    std::unordered_map<int, std::unordered_map<int, ContactInformation>> &getplanewallSphereContactInformationList();
    std::unordered_map<int, std::unordered_map<int, ContactInformation>> &getplanewallFiberContactInformationList();

    std::unordered_map<int, std::unordered_map<int, ContactInformation>> &getFiberFiberContactInformationList();

    std::unordered_map<int, std::unordered_map<int, ContactInformation>> &getcylinderwallSphereContactInformationList();
    std::unordered_map<int, std::unordered_map<int, ContactInformation>> &getcylinderwallFiberContactInformationList();

private:
    std::unordered_map<int, std::unordered_map<int, ContactInformation>> spherespherecontactinformationlist;
    std::unordered_map<int, std::unordered_map<int, ContactInformation>> nextspherespherecontactinformationlist;

    std::unordered_map<int, std::unordered_map<int, ContactInformation>> spherefibercontactinformationlist;
    std::unordered_map<int, std::unordered_map<int, ContactInformation>> nextspherefibercontactinformationlist;

    std::unordered_map<int, std::unordered_map<int, ContactInformation>> fiberfibercontactinformationlist;
    std::unordered_map<int, std::unordered_map<int, ContactInformation>> nextfiberfibercontactinformationlist;

    std::unordered_map<int, std::unordered_map<int, ContactInformation>> planewallfibercontactinformationlist;
    std::unordered_map<int, std::unordered_map<int, ContactInformation>> nextplanewallfibercontactinformationlist;

    std::unordered_map<int, std::unordered_map<int, ContactInformation>> planewallspherecontactinformationlist;
    std::unordered_map<int, std::unordered_map<int, ContactInformation>> nextplanewallspherecontactinformationlist;

    std::unordered_map<int, std::unordered_map<int, ContactInformation>> cylinderwallspherecontactinformationlist;
    std::unordered_map<int, std::unordered_map<int, ContactInformation>> nextcylinderwallspherecontactinformationlist;

    std::unordered_map<int, std::unordered_map<int, ContactInformation>> cylinderwallfibercontactinformationlist;
    std::unordered_map<int, std::unordered_map<int, ContactInformation>> nextcylinderwallfibercontactinformationlist;

    std::map<PropertyTypeID, std::map<PropertyTypeID, double>> effectiveyoungsmodulus;

    std::map<PropertyTypeID, std::map<PropertyTypeID, double>> effectiveshearmodulus;

    std::map<PropertyTypeID, std::map<PropertyTypeID, double>> effectiveslidingfriction;

    std::map<PropertyTypeID, std::map<PropertyTypeID, double>> effectiverollingfriction;

    std::map<PropertyTypeID, std::map<PropertyTypeID, double>> effectiveradius;

    std::map<PropertyTypeID, std::map<PropertyTypeID, double>> effectivemass;

    std::map<PropertyTypeID, std::map<PropertyTypeID, double>> modelparameterbeta;

    double energyDissipation;
};