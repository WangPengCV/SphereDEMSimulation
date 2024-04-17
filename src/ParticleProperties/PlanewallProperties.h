#pragma once
#include "ParticleProperties.h"
#include <cmath>

class PlanewallProperties : public ParticleProperties {
public:
    // Use the constructor of the base class
    PlanewallProperties(double density, double thickness,  double rollingFriction, double slidingFriction,
                     double youngModulus, double restitution, double poissonRatio);
    virtual std::string save_tostring() const override;
    void setThickness(double thickness );
    double getThickness() const { return thickness; }
   
private:
    double thickness;
};
