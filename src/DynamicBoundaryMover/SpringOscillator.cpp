#include "SpringOscillator.h"

Eigen::Vector3d SpringOscillator::getForce(double position,double velocity)
{
    Eigen::Vector3d resultantForce = Eigen::Vector3d::Zero();
    double damping_force = -dampingCoefficient * velocity;
    double resistence_force = -springConstant * (position - equilibriumPosition);
    resultantForce.y() = damping_force + resistence_force;

    return resultantForce;
}
void SpringOscillator::setequilibriumPosition(double equilibriumPosition)
{
    this->equilibriumPosition = equilibriumPosition;
}
std::string SpringOscillator::save_tostring() const
{
    std::ostringstream ss;
    ss.precision(std::numeric_limits<double>::digits10 + 1);
    ss << springConstant << ", " << dampingCoefficient << ", " << equilibriumPosition;
    return ss.str();
}
