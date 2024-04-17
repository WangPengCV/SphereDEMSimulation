#include "SphereProperties.h"

SphereProperties::SphereProperties(double density, double mass, double radius, double rollingFriction, double slidingFriction,
                                   double youngModulus, double restitution, double poissonRatio, double moment_of_inertia)
    : ParticleProperties(density,  rollingFriction, slidingFriction,
                         youngModulus, restitution, poissonRatio),
                         radius(radius),mass(mass), momentofinertia(moment_of_inertia)
{
}

void SphereProperties::setMass(double Mass)
{
    mass = Mass;
    validateProperty(mass, "Mass");
}

void SphereProperties::setRadius(double Radius)
{
    radius = Radius;
    validateProperty(radius, "Radius");
}

void SphereProperties::setMomentofinertia(double Momentofinertia)
{
    momentofinertia = Momentofinertia;
    validateProperty(radius, "Radius");
}

std::string SphereProperties::save_tostring() const
{
    std::ostringstream ss;
    ss.precision(std::numeric_limits<double>::digits10 + 1);
    ss << density << ", " << radius << ", " << rolling_friction_coefficient << ", "
       << slide_friction_coefficient << ", " << Young_modulus << ", "
       << restitution << ", " << poisson_ratio;
    return ss.str();
}
