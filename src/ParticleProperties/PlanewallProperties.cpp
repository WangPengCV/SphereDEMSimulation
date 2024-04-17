#include "PlanewallProperties.h"

PlanewallProperties::PlanewallProperties(double density, double thickness, double rollingFriction, double slidingFriction,
                                         double youngModulus, double restitution, double poissonRatio)
    : ParticleProperties(density, rollingFriction, slidingFriction,
                         youngModulus, restitution, poissonRatio),
      thickness(thickness)
{
}

void PlanewallProperties::setThickness(double thickness)
{
    this->thickness = thickness;
}



std::string PlanewallProperties::save_tostring() const
{
    std::ostringstream ss;
    ss.precision(std::numeric_limits<double>::digits10 + 1);
    ss << density << ", " << thickness << ", "
       << rolling_friction_coefficient << ", " << slide_friction_coefficient << ", " << Young_modulus << ", "
       << restitution << ", " << poisson_ratio;
    return ss.str();
}
