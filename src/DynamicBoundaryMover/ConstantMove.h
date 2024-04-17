#pragma once
#include "DynamicBoundaryMover.h"
#include <Eigen/Dense>

class ConstantMove : public DynamicBoundaryMover
{
public:
    ConstantMove(Eigen::Vector3d &velocity, double period)
        : velocity(velocity), period(period) {}

    Eigen::Vector3d getVelocity() const { return velocity; };
    double getPeriod() const { return period; }
    void setVelocity(Eigen::Vector3d &velocity);
   
    std::string save_tostring() const;

private:
    Eigen::Vector3d velocity;
    double period;
};
