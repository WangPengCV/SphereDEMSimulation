#pragma once
#include "Particle.h"

class SphereParticle : public Particle
{
public:
    SphereParticle(int id, PropertyTypeID type, int state,
                   const std::shared_ptr<ParticlePropertyManager>& manager,
                   const Eigen::Vector3d &position,
                   const Eigen::Vector3d &velocity = Eigen::Vector3d::Zero(),
                   const Eigen::Vector3d &omega = Eigen::Vector3d::Zero(),
                   const Eigen::Vector3d &force = Eigen::Vector3d::Zero(),
                   const Eigen::Vector3d &torque = Eigen::Vector3d::Zero());
    SphereParticle();
    virtual ~SphereParticle();

    double computeOverlap(const std::shared_ptr<PlaneWall>& planewall);



    virtual std::string save_tostring() const override; 

    std::string save_fibertostring() const; 

    virtual void updateVelocity(double deltaTime, Eigen::Vector3d& gravity, double mass) override;
    virtual void updateOmega(double deltaTime,double moment_of_inertia) override;
};
