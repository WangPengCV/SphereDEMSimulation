#pragma once
#include "ParticlePropertyManager.h"
class Fiber
{
public:
    Fiber(int id, PropertyTypeID type, int startstate, int endstate, std::shared_ptr<ParticlePropertyManager>& manager);
    const std::vector<int>& getSphereCylinderBond() const;
    const PropertyTypeID& getType() const {return type;}

    void setSphereCylinderBond( const std::vector<int>& spherecylinderbond);
    std::string save_tostring() const;
private:
    int id;
    PropertyTypeID type;
    int startstate;
    int endstate;
    std::shared_ptr<ParticlePropertyManager> manager;
    std::vector<int> spherecylinderbond;

};