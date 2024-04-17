#include "PlaneWall.h"

void PlaneWall::generateMesh(double meshResolution)
{
    meshResolution = 0.99 * meshResolution;
    meshVertices.clear(); // Clear existing vertices

    // Calculate the basis vectors for the plane
    Eigen::Vector3d u = (corner2 - corner1).normalized();
    Eigen::Vector3d v = (corner3 - corner2).normalized();

    // Determine the number of divisions along each basis vector
    double lengthU = (corner2 - corner1).norm();
    double lengthV = (corner3 - corner2).norm();
    int divisionsU = static_cast<int>(lengthU / meshResolution);
    int divisionsV = static_cast<int>(lengthV / meshResolution);

    // Generate grid points
    for (int i = 0; i <= divisionsU; ++i)
    {
        for (int j = 0; j <= divisionsV; ++j)
        {
            Eigen::Vector3d point = corner1 + u * (meshResolution * i) + v * (meshResolution * j);
            meshVertices.push_back(point);
        }
    }
}
const std::vector<Eigen::Vector3d> &PlaneWall::getMeshVertices() const
{
    return meshVertices;
}
const Eigen::Vector3d &PlaneWall::getNormal() const
{
    return normal;
}
const Eigen::Vector3d &PlaneWall::getCorner1() const
{
    return corner1;
}
const Eigen::Vector3d &PlaneWall::getCorner2() const
{
    return corner2;
}
const Eigen::Vector3d &PlaneWall::getCorner3() const
{
    return corner3;
}
const Eigen::Vector3d &PlaneWall::getVelocity() const
{
    return velocity;
}
void PlaneWall::setNormal(const Eigen::Vector3d &Normal)
{
    normal = Normal;
}
void PlaneWall::setCorner1(const Eigen::Vector3d &Corner1)
{
    corner1 = Corner1;
}
void PlaneWall::setCorner2(const Eigen::Vector3d &Corner2)
{
    corner2 = Corner2;
}
void PlaneWall::setCorner3(const Eigen::Vector3d &Corner3)
{
    corner3 = Corner3;
}
void PlaneWall::setVelocity(const Eigen::Vector3d &Velociy)
{
    velocity = Velociy;
}
void PlaneWall::addForce(const Eigen::Vector3d &additionalForce)
{
    force += additionalForce;
}
void PlaneWall::resetForce()
{
    force.setZero();
}
void PlaneWall::setMass(std::shared_ptr<ParticlePropertyManager> &manager)
{
    double thickness = manager->getPlanewallProperties(this->getType())->getThickness();
    Eigen::Vector3d u = (corner2 - corner1);
    Eigen::Vector3d v = (corner3 - corner2);
    mass = manager->getPlanewallProperties(this->getType())->getDensity() * (u.norm() * v.norm() * thickness);
}
std::string PlaneWall::save_tostring() const
{
    std::ostringstream ss;
    ss.precision(std::numeric_limits<double>::digits10 + 1);
    ss << "PLANEWALL, " << id << ", " << type.getCategory() << ", " << type.getSubType() << ", " << state << ", "
       << normal.x() << ", " << normal.y() << ", " << normal.z() << ", "
       << corner1.x() << ", " << corner1.y() << ", " << corner1.z() << ", "
       << corner2.x() << ", " << corner2.y() << ", " << corner2.z() << ", "
       << corner3.x() << ", " << corner3.y() << ", " << corner3.z() << ", "
       << velocity.x() << ", " << velocity.y() << ", " << velocity.z();
    return ss.str();
}