#pragma once
#include "BoundaryCondition.h"
#include "PlaneWall.h"
#include "ParticlePropertyManager.h"
#include <unordered_map>
#include <Eigen/Dense>

class RectangularContainer : public BoundaryCondition {
public:
    RectangularContainer(){};
    // Constructor
    RectangularContainer(int id, const PropertyTypeID& type, int state, 
                         const Eigen::Vector3d& dimensions, 
                         double rotation,
                         const Eigen::Vector3d &Center,
                         const Eigen::Vector3d& velocity = Eigen::Vector3d::Zero(),
                         const Eigen::Vector3d& force = Eigen::Vector3d::Zero());

    void rotateContainer(double angleDegrees, const Eigen::Vector3d& axis);

    void translateToCenter(const Eigen::Vector3d& simulationDomain);


    void setMass(std::shared_ptr<ParticlePropertyManager> manager);

    void setVelocity(const Eigen::Vector3d& Velocity);

    void setCenter(const Eigen::Vector3d& Center);

    void addForce(const Eigen::Vector3d& Force);

    void resetForce();

    const std::vector<std::unique_ptr<PlaneWall>>& getPlaneWall() const {return walls;} 

    // Override the save_tostring method from the BoundaryCondition class
    virtual std::string save_tostring() const override;

    // Methods to access individual walls
    const std::unique_ptr<PlaneWall>& getTopWall() const { return walls[5];}
    const std::unique_ptr<PlaneWall>& getBottomWall() const {return walls[2];} 
    const std::unique_ptr<PlaneWall>& getLeftWall() const {return walls[0];}
    const std::unique_ptr<PlaneWall>& getRightWall() const {return walls[1];}
    const std::unique_ptr<PlaneWall>& getFrontWall() const {return walls[3];}
    const std::unique_ptr<PlaneWall>& getBackWall() const {return walls[4];}
    double getMass() const{return mass;}
    const Eigen::Vector3d& getVelocity() const {return velocity;}
    const Eigen::Vector3d& getForce() const {return force;}
    const Eigen::Vector3d& getCenter() const {return center;}
    double getRotation() const {return rotation;};
    const Eigen::Vector3d& getDimensions() const {return dimensions;}

private:
    std::vector<std::unique_ptr<PlaneWall>> walls;
    //Eigen::Vector3d lowerCorner;
    Eigen::Vector3d dimensions;
    Eigen::Vector3d center;

    double rotation;
    double mass;
    Eigen::Vector3d force;
    Eigen::Vector3d velocity;

    // Helper methods to create walls
    void createWalls();


};


