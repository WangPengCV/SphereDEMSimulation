
#pragma once
#include "ParticleProperties.h"
#include <cmath>

class CylinderwallProperties : public ParticleProperties {
public:
    // Use the constructor of the base class
    CylinderwallProperties(double density, double thickness,  double rollingFriction, double slidingFriction,
                     double youngModulus, double restitution, double poissonRatio, bool hollow = true);
                     
    virtual std::string save_tostring() const override;

    void setThickness(double thickness );

    double getThickness() const { return thickness; }

    bool getHollow() const { return hollow; }
private:
    double thickness;
    bool hollow;
};