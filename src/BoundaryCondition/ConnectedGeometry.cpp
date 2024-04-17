#include "ConnectedGeometry.h"

ConnectedGeometry::ConnectedGeometry(int id, Eigen::Vector3d &center, Eigen::Vector3d &velocity, std::bitset<3> dofs)
    : id(id), centerOfMass(center), velocity(velocity), dofs(dofs)
{
}
void ConnectedGeometry::setDegreeOfFreedom(DegreeOfFreedom dof, bool value)
{
    dofs[static_cast<size_t>(dof)] = value;
}
void ConnectedGeometry::setMass(double mass)
{
    this->mass = mass;
}

void ConnectedGeometry::addGeometry(const std::string &geoType, int geometryId)
{
    geometrys.push_back(std::make_tuple(geoType, geometryId));
}

void ConnectedGeometry::setGeometry(const std::vector<std::tuple<std::string, int>> &geometrys)
{
    this->geometrys = geometrys;
}

void ConnectedGeometry::setVelocity(const Eigen::Vector3d &velocity)
{
    this->velocity = velocity;
}

void ConnectedGeometry::setCenterOfMass(const Eigen::Vector3d &centerOfMass)
{
    this->centerOfMass = centerOfMass;
}

void ConnectedGeometry::addForce(const Eigen::Vector3d &force)
{
    this->force = force;
}

int ConnectedGeometry::GeometrySize() const
{
    return static_cast<int>(geometrys.size());
}

void ConnectedGeometry::resetForce()
{
    force = Eigen::Vector3d::Zero();
}
std::string ConnectedGeometry::save_tostring() const
{
    std::ostringstream ss;
    ss.precision(std::numeric_limits<double>::digits10 + 1);
    ss << id << ", " << geometrys.size();
    for (int i = 0; i < geometrys.size(); ++i)
    {
        ss << ", " << std::get<0>(geometrys[i]) << ", " << std::get<1>(geometrys[i]);
    }
    if (getDegreeOfFreedom(DegreeOfFreedom::X) && getDegreeOfFreedom(DegreeOfFreedom::Y) && getDegreeOfFreedom(DegreeOfFreedom::Z))
    {
    }
    else
    {
        if (getDegreeOfFreedom(DegreeOfFreedom::X))
        {
            ss << ", " << "X";
        }
        else if (getDegreeOfFreedom(DegreeOfFreedom::Y))
        {
            ss << ", " << "Y";
        }
        else if (getDegreeOfFreedom(DegreeOfFreedom::Z))
        {
            ss << ", " << "Z";
        }
    }

    return ss.str();
}
