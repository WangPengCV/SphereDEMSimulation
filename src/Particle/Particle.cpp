#include "Particle.h"

Particle::Particle(int id, PropertyTypeID type, int state,
                   const std::shared_ptr<ParticlePropertyManager> &manager,
                   const Eigen::Vector3d &position,
                   const Eigen::Vector3d &velocity, const Eigen::Vector3d &omega,
                   const Eigen::Vector3d &force, const Eigen::Vector3d &torque)
    : id(id), type(type), state(state), position(position), manager(manager),
      velocity(velocity), omega(omega), force(force), torque(torque)
{
}
Particle::Particle()
{
}
Particle::~Particle()
{
    // Destructor implementation (if needed)
}

void Particle::setPosition(const Eigen::Vector3d &newPosition)
{
    position = newPosition;
}

void Particle::setVelocity(const Eigen::Vector3d &velocity)
{
    this->velocity = velocity;
}

void Particle::setOmega(const Eigen::Vector3d &omega)
{
    this->omega = omega;
}

void Particle::setId(int newId)
{
    id = newId;
}

void Particle::addForce(const Eigen::Vector3d &additionalForce)
{
    force += additionalForce;
}

void Particle::resetForce()
{
    force.setZero();
}

void Particle::addTorque(const Eigen::Vector3d &additionalTorque)
{
    torque += additionalTorque;
}

void Particle::resetTorque()
{
    torque.setZero();
}

void Particle::updatePosition(double deltaTime)
{
    position += velocity * deltaTime;
}
