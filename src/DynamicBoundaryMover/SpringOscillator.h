#pragma once
#include "DynamicBoundaryMover.h"
#include <Eigen/Dense>

class SpringOscillator : public DynamicBoundaryMover
 {
public:
    SpringOscillator(double k, double c, double equilibriumPos )
    : springConstant(k), dampingCoefficient(c), equilibriumPosition(equilibriumPos) {}

    Eigen::Vector3d getForce(double position, double velocity);

   double getequilibriumPosition() const {return equilibriumPosition;}
   void setequilibriumPosition(double equilibriumPosition);
   std::string save_tostring() const;

private:
    double springConstant;
    double dampingCoefficient;
    double equilibriumPosition;
};
