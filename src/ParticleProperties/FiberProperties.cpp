#include "FiberProperties.h"

FiberProperties::FiberProperties(double density, double rollingFriction, double slidingFriction,
                                 double youngModulus, double restitution, double poissonRatio, double radius,
                                 double elementlength, double normalmodulus, double shearmodulus, double twistmodulus, double bendingmodulus,
                                 int nodenumber, double aspectratio, double bonddampingcoefficient) : ParticleProperties(density, rollingFriction, slidingFriction, youngModulus, restitution, poissonRatio),
                                                                                                      radius(radius), elementlength(elementlength),
                                                                                                      normalmodulus(normalmodulus), shearmodulus(shearmodulus), twistmodulus(twistmodulus),
                                                                                                      bendingmodulus(bendingmodulus), nodenumber(nodenumber), aspectratio(aspectratio), bonddampingcoefficient(bonddampingcoefficient)
{
    nodemass = (4.0 / 3.0 * PI * radius * radius * radius + elementlength * (nodenumber - 1) * PI * radius * radius) * density / nodenumber;
    nodemomentofinertia = 2 * nodemass * radius * radius / 5.;


    double sectionarea = PI * radius * radius;
    double node_area_moment_of_inertia = PI * radius * radius * radius * radius / 4.;
    double node_polar_area_moment_of_inertia = node_area_moment_of_inertia * 2;

    normalstiffnesses = normalmodulus * sectionarea / elementlength;
    shearstiffnesses = shearmodulus * sectionarea / elementlength;
    twiststiffnesses = twistmodulus * node_polar_area_moment_of_inertia / elementlength;
    bendingstiffnesses = bendingmodulus * node_area_moment_of_inertia / elementlength;
}

std::string FiberProperties::save_tostring() const
{
    std::ostringstream ss;
    ss.precision(std::numeric_limits<double>::digits10 + 1);
    ss << density << ", " << radius << ", " << rolling_friction_coefficient << ", "
       << slide_friction_coefficient << ", " << Young_modulus << ", "
       << restitution << ", " << poisson_ratio << ", " << normalmodulus << ", " << shearmodulus << ", "
       << twistmodulus << ", " << bendingmodulus << ", " << nodenumber << ", " << aspectratio << ", " << bonddampingcoefficient;
    return ss.str();
}
