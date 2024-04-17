#pragma once
#include "ParticleProperties.h"

class FiberProperties : public ParticleProperties
{
public:
    FiberProperties(double density, double rollingFriction, double slidingFriction,
                    double youngModulus, double restitution, double poissonRatio, double radius,
                    double elementlength, double normalmodulus, double shearmodulus, double twistmodulus, double bendingmodulus,
                    int nodenumber, double aspectratio, double bonddampingcoefficient);

    double getNodemass() const { return nodemass; }

    double getNodemomentofinertia() const { return nodemomentofinertia; }
    double getElementlength() const { return elementlength; }
    double getAspectratio() const { return aspectratio; }
    double getNormalstiffnesses() const
    {
        return normalstiffnesses;
    }
    double getNormalmodulus() const { return normalmodulus; }

    double getShearstiffnesses() const { return shearstiffnesses; }
    double getShearmodulus() const { return shearmodulus; }

    double getTwiststiffnesses() const { return twiststiffnesses; }
    double getTwistmodulus() const { return twistmodulus; }

    double getBendingstiffnesses() const { return bendingstiffnesses; }
    double getBendingmodulus() const { return bendingmodulus; }

    int getNodenumber() const { return nodenumber; }
    double getRadius() const
    {
        return radius;
    }
    double getBonddampingcoefficient() const { return bonddampingcoefficient; }

    std::string save_tostring() const override;

private:
    double aspectratio;

    double nodemomentofinertia;

    double radius;

    double nodemass;

    double elementlength;

    double normalstiffnesses;

    double normalmodulus;

    double shearstiffnesses;

    double shearmodulus;

    double twiststiffnesses;

    double twistmodulus;

    double bendingstiffnesses;

    double bendingmodulus;

    int nodenumber;

    double bonddampingcoefficient;
};
