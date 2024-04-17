#include "Fiber.h"

Fiber::Fiber(int id, PropertyTypeID type, int startstate, int endstate, std::shared_ptr<ParticlePropertyManager>& manager)
: id(id),type(type),startstate(startstate),endstate(endstate),manager(manager)
{

}

void Fiber::setSphereCylinderBond(const std::vector<int>& spherecylinderbond)
{
    this->spherecylinderbond = spherecylinderbond;
}
const std::vector<int>& Fiber::getSphereCylinderBond() const
{
    return spherecylinderbond;
}
std::string Fiber::save_tostring() const
{
    std::ostringstream ss;
    ss.precision(std::numeric_limits<double>::digits10 + 1);

    ss << "Fiber, "  << id << ", " << type.getCategory() << ", " << type.getSubType() << ", " << startstate << ", " << endstate;
    for(int i = 0; i < spherecylinderbond.size(); ++i)
    {
        ss << ", " << spherecylinderbond[i];
    }
    return ss.str();
}

