#include "RectangularContainer.h"

// Implementation of RectangularContainer constructor
RectangularContainer::RectangularContainer(int Id, const PropertyTypeID &Type, int State,
                                           const Eigen::Vector3d &Dimensions,
                                           double Rotation,
                                           const Eigen::Vector3d &Center,
                                           const Eigen::Vector3d &Velocity,
                                           const Eigen::Vector3d &Force)
    : BoundaryCondition(Id, Type, State), dimensions(Dimensions), rotation(Rotation), center(Center),
      force(Force), velocity(Velocity), mass(0)
{
    walls.resize(6);
    createWalls();
    
}

void RectangularContainer::setMass(std::shared_ptr<ParticlePropertyManager> manager)
{
    double thickness = manager->getPlanewallProperties(this->getType())->getThickness();
    mass = manager->getPlanewallProperties(this->getType())->getDensity() * (2 * dimensions.x() * dimensions.z() * thickness + 2 * dimensions.x() * dimensions.y() * thickness + 2 * dimensions.y() * dimensions.z() * thickness);
}

void RectangularContainer::setCenter(const Eigen::Vector3d& newCenter)
{
    center = newCenter;
}

void RectangularContainer::rotateContainer(double angleDegrees, const Eigen::Vector3d &axis)
{
    // Convert angle to radians
    double angleRadians = angleDegrees * PI / 180.0;

    // Create a rotation matrix
    Eigen::AngleAxisd rotationMatrix(angleRadians, axis.normalized());
    // Rotate each wall
    for (auto &wall : walls)
    {
        // Rotate each corner point of the wall
        Eigen::Vector3d corner1Rotated = rotationMatrix * wall->getCorner1();
        Eigen::Vector3d corner2Rotated = rotationMatrix * wall->getCorner2();
        Eigen::Vector3d corner3Rotated = rotationMatrix * wall->getCorner3();
        Eigen::Vector3d normalRotated = rotationMatrix * wall->getNormal();

        // Update the wall with new corners (and possibly normal if needed)
        wall = std::make_unique<PlaneWall>(wall->getId(), wall->getType(), wall->getState(),
                                           normalRotated, // Might need to rotate normal as well
                                           corner1Rotated, corner2Rotated, corner3Rotated,
                                           wall->getVelocity());
    }
    center = rotationMatrix*center;
}

void RectangularContainer::translateToCenter(const Eigen::Vector3d& simulationDomain)
{
    Eigen::Vector3d tranlateVector = simulationDomain / 2 - center;
    for (auto &wall : walls)
    {
        // tanslate each corner point of the wall
        Eigen::Vector3d corner1Rotated = wall->getCorner1() + tranlateVector;
        Eigen::Vector3d corner2Rotated = wall->getCorner2() + tranlateVector;
        Eigen::Vector3d corner3Rotated = wall->getCorner3() + tranlateVector;
       

        // Update the wall with new corners 
        wall = std::make_unique<PlaneWall>(wall->getId(), wall->getType(), wall->getState(),
                                           wall->getNormal(), 
                                           corner1Rotated, corner2Rotated, corner3Rotated,
                                           wall->getVelocity());
    }
    center = simulationDomain / 2;


}

void RectangularContainer::createWalls()
{
    // Assuming y is vertical, x and z are horizontal dimensions
    Eigen::Vector3d topNormal(0, -1, 0);   // Normal pointing upwards
    Eigen::Vector3d bottomNormal(0, 1, 0); // Normal pointing downwards
    Eigen::Vector3d leftNormal(1, 0, 0);   // Normal pointing left
    Eigen::Vector3d rightNormal(-1, 0, 0); // Normal pointing right
    Eigen::Vector3d frontNormal(0, 0, -1); // Normal pointing front
    Eigen::Vector3d backNormal(0, 0, 1);   // Normal pointing back

    // Create each wall with the appropriate corners and normals
    // You'll need to calculate the corners based on the lowerCorner and dimensions

    Eigen::Vector3d lowerCorner = Eigen::Vector3d(center.x() - dimensions.x() / 2.,center.y() - dimensions.y() / 2.,center.z() - dimensions.z() / 2. );
    Eigen::Vector3d corner1 = lowerCorner;
    Eigen::Vector3d corner2 = lowerCorner + Eigen::Vector3d(0, 0, dimensions.z());
    Eigen::Vector3d corner3 = corner2 + Eigen::Vector3d(0, dimensions.y(), 0);
    Eigen::Vector3d corner4 = lowerCorner + Eigen::Vector3d(dimensions.x(), 0, 0);
    Eigen::Vector3d corner5 = corner4 + Eigen::Vector3d(0, 0, dimensions.z());
    Eigen::Vector3d corner6 = corner5 + Eigen::Vector3d(0, dimensions.y(), 0);
    Eigen::Vector3d corner7 = corner4 + Eigen::Vector3d(0, dimensions.y(), 0);
    Eigen::Vector3d corner8 = corner1 + Eigen::Vector3d(0, dimensions.y(), 0);

    walls[0] = std::make_unique<PlaneWall>(0, type, state, leftNormal, corner1, corner2, corner3, Eigen::Vector3d::Zero());
    walls[1] = std::make_unique<PlaneWall>(1, type, state, rightNormal, corner4, corner5, corner6, Eigen::Vector3d::Zero());
    walls[2] = std::make_unique<PlaneWall>(2, type, state, bottomNormal, corner4, corner5, corner2, Eigen::Vector3d::Zero());
    walls[3] = std::make_unique<PlaneWall>(3, type, state, frontNormal, corner5, corner2, corner3, Eigen::Vector3d::Zero());
    walls[4] = std::make_unique<PlaneWall>(4, type, state, backNormal, corner4, corner1, corner8, Eigen::Vector3d::Zero());
    walls[5] = std::make_unique<PlaneWall>(5, type, state, topNormal, corner7, corner6, corner3, Eigen::Vector3d::Zero());
}

void RectangularContainer::addForce(const Eigen::Vector3d &additionForce)
{
    force += additionForce;
}

void RectangularContainer::setVelocity(const Eigen::Vector3d &Velocity)
{
    velocity = Velocity;
}

void RectangularContainer::resetForce()
{
    force.setZero();
}

std::string RectangularContainer::save_tostring() const
{
    std::ostringstream ss;
    ss.precision(std::numeric_limits<double>::digits10 + 1);
    ss << "RECTANGULARCONTAINER, " << id << ", " << type.getCategory() << ", " << type.getSubType() << ", " << state << ", "
       << dimensions.x() << ", " << dimensions.y() << ", " << dimensions.z() << ", "
       << rotation << ", "
       << center.x() << ", " << center.y() << ", " << center.z() << ", "
       << velocity.x() << ", " << velocity.y() << ", " << velocity.z();
    return ss.str();
}