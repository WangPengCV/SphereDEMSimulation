#include "SphereParticle.h"

SphereParticle::SphereParticle(int id, PropertyTypeID type, int state,
                               const std::shared_ptr<ParticlePropertyManager>& manager,
                               const Eigen::Vector3d &position,
                               const Eigen::Vector3d &velocity, const Eigen::Vector3d &omega,
                               const Eigen::Vector3d &force, const Eigen::Vector3d &torque)
    : Particle(id, type, state, manager, position, velocity, omega, force, torque)
{
}
SphereParticle::SphereParticle()
{
    
}
SphereParticle::~SphereParticle()
{
    // Custom cleanup for SphereParticle, if needed
}

// void SphereParticle::updateVelocity(double deltaTime, Eigen::Vector3d& gravity)
// {
//     // Implement sphere-specific velocity update logic
//     double mass = manager->getSphereProperties(type)->getMass();
//     Eigen::Vector3d acceleration = force / mass + gravity;
//     velocity += acceleration * deltaTime;
// }

void SphereParticle::updateVelocity(double deltaTime, Eigen::Vector3d& gravity,double mass)
{
    Eigen::Vector3d acceleration = force / mass + gravity;
    velocity += acceleration * deltaTime;
}


void SphereParticle::updateOmega(double deltaTime, double moment_of_inertia)
{
    // Implement sphere-specific angular velocity update logic
    //double moment_of_inertia = manager->getSphereProperties(type)->getMomentOfInertia();
    Eigen::Vector3d angular_acceleration = torque / moment_of_inertia;
    omega += angular_acceleration * deltaTime;
}



double SphereParticle::computeOverlap(const std::shared_ptr<PlaneWall> &planewall)
{
    // Retrieve the sphere's center and radius
    Eigen::Vector3d sphereCenter = this->getPosition();
    double sphereRadius = manager->getSphereProperties(this->getType())->getRadius();

    // Retrieve the plane wall's point and normal
    Eigen::Vector3d planePoint = planewall->getCorner1();
    
    Eigen::Vector3d planeNormal = planewall->getNormal();

    // Compute the vector from a point on the plane to the sphere's center
    Eigen::Vector3d vecToSphereCenter = sphereCenter - planePoint;

    // Compute the distance from the sphere's center to the plane
    double distanceToPlane = vecToSphereCenter.dot(planeNormal);

    // Calculate the overlap (penetration depth)
    double overlap = sphereRadius - std::abs(distanceToPlane);

    return (overlap > 0.0) ? overlap : 0.0;
}



std::string SphereParticle::save_tostring() const {
    std::ostringstream ss;
    ss.precision(std::numeric_limits<double>::digits10 + 1);

    ss << "PARTICLE, "  << "SPHERE, "<< id << ", " << type.getCategory() << ", " << type.getSubType() << ", " << state << ", "
       << position.x() << ", " << position.y() << ", " << position.z() << ", "
       << velocity.x() << ", " << velocity.y() << ", " << velocity.z() << ", "
       << omega.x() << ", " << omega.y() << ", " << omega.z() ;
    return ss.str();
}

std::string SphereParticle::save_fibertostring() const 
{
    std::ostringstream ss;
    ss.precision(std::numeric_limits<double>::digits10 + 1);

    ss << "Fiber, "  << "SPHERE, "<< id << ", " << type.getCategory() << ", " << type.getSubType() << ", " << state << ", "
       << position.x() << ", " << position.y() << ", " << position.z() << ", "
       << velocity.x() << ", " << velocity.y() << ", " << velocity.z() << ", "
       << omega.x() << ", " << omega.y() << ", " << omega.z() ;
    return ss.str();
}