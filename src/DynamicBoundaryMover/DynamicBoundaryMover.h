#pragma once
#include <string>
class DynamicBoundaryMover {
public:
    virtual ~DynamicBoundaryMover() = default;
    virtual std::string save_tostring() const = 0;  

};
