#include "CylinderwallProperties.h"

CylinderwallProperties::CylinderwallProperties(double density, double thickness, double rollingFriction, double slidingFriction,
                                               double youngModulus, double restitution, double poissonRatio, bool hollow)
    : ParticleProperties(density, rollingFriction, slidingFriction,
                         youngModulus, restitution, poissonRatio),
      thickness(thickness), hollow(hollow)
{
}

void CylinderwallProperties::setThickness(double thickness)
{
    this->thickness = thickness;
}

std::string CylinderwallProperties::save_tostring() const
{
    std::ostringstream ss;
    ss.precision(std::numeric_limits<double>::digits10 + 1);
    ss << density << ", " << thickness << ", "
       << rolling_friction_coefficient << ", " << slide_friction_coefficient << ", " << Young_modulus << ", "
       << restitution << ", " << poisson_ratio << ", "  << hollow;
    return ss.str();
}
