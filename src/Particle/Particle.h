#pragma once
#include <Eigen/Dense>
#include "ParticlePropertyManager.h"
#include "PlaneWall.h"
#include <memory>
class Particle
{
public:
    Particle(int id, PropertyTypeID type, int state,
             const std::shared_ptr<ParticlePropertyManager>& manager,
             const Eigen::Vector3d &position,
             const Eigen::Vector3d &velocity = Eigen::Vector3d::Zero(),
             const Eigen::Vector3d &omega = Eigen::Vector3d::Zero(),
             const Eigen::Vector3d &force = Eigen::Vector3d::Zero(),
             const Eigen::Vector3d &torque = Eigen::Vector3d::Zero());
    Particle();
    virtual ~Particle();

    virtual void updateVelocity(double deltaTime, Eigen::Vector3d& gravity, double mass) = 0;

    virtual void updateOmega(double deltaTime,double moment_of_inertia) = 0;

    virtual std::string save_tostring() const = 0;  

    void setPosition(const Eigen::Vector3d &position);

    void setVelocity(const Eigen::Vector3d &velocity);

    void setOmega(const Eigen::Vector3d &omega);

    void setId(int id);


    void addForce(const Eigen::Vector3d &additionalForce);

    void resetForce();

    void addTorque(const Eigen::Vector3d &additionalTorque);

    void resetTorque();

    void updatePosition(double deltaTime);

    // Accessor methods (if needed)
    int getId() const { return id; }    
    const Eigen::Vector3d& getPosition() const { return position; }
    const Eigen::Vector3d& getVelocity() const { return velocity; }
    const Eigen::Vector3d& getOmega() const { return omega; }
    const Eigen::Vector3d& getForce() const { return force; }
    const Eigen::Vector3d& getTorque() const { return torque; }


    const PropertyTypeID& getType() const {return type;}
    int getState() const {return state;}
    std::shared_ptr<ParticlePropertyManager>  getParticlePropertyManager() const { return manager;}


    // ... other accessors ...

protected:
    int id;
    PropertyTypeID type;
    int state;

    Eigen::Vector3d position;
    Eigen::Vector3d velocity;
    Eigen::Vector3d omega;
    Eigen::Vector3d force;
    Eigen::Vector3d torque;

    std::shared_ptr<ParticlePropertyManager> manager;
};
