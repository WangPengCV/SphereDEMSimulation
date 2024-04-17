#include "CylinderContainer.h"

CylinderContainer::CylinderContainer(int Id, const PropertyTypeID &Type, int State,
                                     double Radius,
                                     const Eigen::Vector3d &Startpoint,
                                     const Eigen::Vector3d &Endpoint,
                                     const Eigen::Vector3d &Velocity,
                                     const Eigen::Vector3d &Force)
    : BoundaryCondition(Id, Type, State), radius(Radius), startpoint(Startpoint), endpoint(Endpoint), force(Force), velocity(Velocity), mass(0)

{
}
void CylinderContainer::setMass(std::shared_ptr<ParticlePropertyManager> manager)
{
    const auto &pro = manager->getCylinderwallProperties(this->getType());
    if (pro->getHollow())
    {
        double thickness = manager->getPlanewallProperties(this->getType())->getThickness();
        double h = (startpoint - endpoint).norm();
        mass = pro->getDensity() * (PI * (radius + thickness) * (radius + thickness) * h - PI * radius * radius * h);
    }
    else
    {
        double h = (startpoint - endpoint).norm();
        mass = pro->getDensity() * PI * radius * radius * h;
    }
}

void CylinderContainer::setVelocity(const Eigen::Vector3d &velocity)
{
    this->velocity = velocity;
}

void CylinderContainer::setStartpoint(const Eigen::Vector3d &startpoint)
{
    this->startpoint = startpoint;
}
void CylinderContainer::setEndpoint(const Eigen::Vector3d &endpoint)
{
    this->endpoint = endpoint;
}
void CylinderContainer::addForce(const Eigen::Vector3d &additionalForce)
{
    force += additionalForce;
}
void CylinderContainer::resetForce()
{
    force.setZero();
}
void CylinderContainer::generateMesh(double meshResolution)
{
    meshResolution = 0.5 * meshResolution;
    double circumference = 2.0 * PI * radius;

    // Calculate the number of angular divisions based on the desired resolution
    int angularDivisions = static_cast<int>(ceil(circumference / meshResolution)) + 3;

    meshVertices.clear(); // Clear existing vertices

    // Calculate the axis direction of the cylinder
    Eigen::Vector3d axis = (endpoint - startpoint).normalized();

    // Calculate the number of divisions along the axis
    double length = (endpoint - startpoint).norm();
    int axialDivisions = static_cast<int>(length / meshResolution);

    // Calculate a basis vector perpendicular to the axis
    Eigen::Vector3d perpendicular;
    if (axis.x() == 0.0 && axis.y() == 0.0)
    {
        perpendicular = Eigen::Vector3d(1.0, 0.0, 0.0); // If axis is along z-axis, choose x-axis as perpendicular
    }
    else
    {
        perpendicular = Eigen::Vector3d(0.0, 0.0, 1.0); // Otherwise choose z-axis as perpendicular
    }
    Eigen::Vector3d u = axis.cross(perpendicular).normalized(); // Basis vector u perpendicular to axis
    Eigen::Vector3d v = axis.cross(u).normalized();             // Basis vector v perpendicular to axis and u

    // Generate grid points for each division along the axis
    for (int i = 0; i <= axialDivisions; ++i)
    {
        double axialPosition = static_cast<double>(i) / static_cast<double>(axialDivisions);

        // Generate grid points for each point on the circumference of the cylinder
        for (int j = 0; j < angularDivisions; ++j)
        {
            double angle = 2.0 * PI * static_cast<double>(j) / static_cast<double>(angularDivisions);
            Eigen::Vector3d point = startpoint + axialPosition * (endpoint - startpoint);
            point += radius * (cos(angle) * u + sin(angle) * v);
            meshVertices.push_back(point);
        }
    }
}
std::string CylinderContainer::save_tostring() const
{
    std::ostringstream ss;
    ss.precision(std::numeric_limits<double>::digits10 + 1);
    ss << "CYLINDERCONTAINER, " << id << ", " << type.getCategory() << ", " << type.getSubType() << ", " << state << ", " << radius << ", "
       << startpoint.x() << ", " << startpoint.y() << ", " << startpoint.z() << ", "
       << endpoint.x() << ", " << endpoint.y() << ", " << endpoint.z() << ", "
       << velocity.x() << ", " << velocity.y() << ", " << velocity.z();
    return ss.str();
}