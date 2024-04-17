#pragma once
#include "BoundaryCondition.h"

class PlaneWall : public BoundaryCondition
{
public:
    PlaneWall(){}
    PlaneWall(int Id, const PropertyTypeID &Type, int State, const Eigen::Vector3d &Normal, const Eigen::Vector3d &Corner1,
              const Eigen::Vector3d &Corner2, const Eigen::Vector3d &Corner3, const Eigen::Vector3d &Velocity,const Eigen::Vector3d &Force = Eigen::Vector3d::Zero())
        : BoundaryCondition(Id, Type,State), normal(Normal), corner1(Corner1), corner2(Corner2), corner3(Corner3), velocity(Velocity), force(Force){}
          
    virtual std::string save_tostring() const override; 
    void setNormal(const Eigen::Vector3d &normal);
    void setCorner1(const Eigen::Vector3d &corner1);
    void setCorner2(const Eigen::Vector3d &corner2);
    void setCorner3(const Eigen::Vector3d &corner3);
    void setVelocity(const Eigen::Vector3d &velociy);
    void setMass(std::shared_ptr<ParticlePropertyManager>& manager);
    void  generateMesh(double meshResolution);
    const std::vector<Eigen::Vector3d> &getMeshVertices() const;
    const Eigen::Vector3d& getNormal() const;
    const Eigen::Vector3d &getCorner1() const;
    const Eigen::Vector3d &getCorner2() const;
    const Eigen::Vector3d &getCorner3() const;
    const Eigen::Vector3d &getVelocity() const;
    double getMass()   const{return mass;};
    void addForce(const Eigen::Vector3d &additionalForce);
    void resetForce();
    const Eigen::Vector3d& getForce() const{
        return force;
    }



private:

    Eigen::Vector3d normal;                    // Normal to the plane
    Eigen::Vector3d corner1, corner2, corner3; // Three corner points of a plane wall (arranged clockwise)
    Eigen::Vector3d velocity;
    Eigen::Vector3d force;
    double mass;
    std::vector<Eigen::Vector3d> meshVertices; // discreted points for contact detection

};
