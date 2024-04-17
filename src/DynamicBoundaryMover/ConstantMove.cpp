#include "ConstantMove.h"

void ConstantMove::setVelocity(Eigen::Vector3d& velocity)
{
    this->velocity = velocity;
}
std::string ConstantMove::save_tostring() const
{
    std::ostringstream ss;
    ss.precision(std::numeric_limits<double>::digits10 + 1);
    ss << velocity.x() << ", " << velocity.y() << ", " << velocity.z() << ", " << period;
    return ss.str();
}