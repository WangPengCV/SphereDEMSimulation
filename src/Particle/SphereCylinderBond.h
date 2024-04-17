// https://www.geometrictools.com/Documentation/Documentation.html
#pragma once
#include <Eigen/Dense>
#include "ParticlePropertyManager.h"
#include "SphereParticle.h"

class SphereCylinderBond
{
public:
    SphereCylinderBond(int id, PropertyTypeID type,
                       const std::shared_ptr<ParticlePropertyManager> &manager,
                       int fibeid,
                       int node1, int node2, int neighborelement1, int neighborelement2,
                       double energyDissipation = 0,
                       const Eigen::Vector3d &tangentialforce = Eigen::Vector3d::Zero(),
                       const Eigen::Vector3d &twisttorque = Eigen::Vector3d::Zero(),
                       const Eigen::Vector3d &bendtorque = Eigen::Vector3d::Zero());
    SphereCylinderBond();
    int getNode1() const;
    int getNode2() const;
    int getNeighborelement1() const;
    int getNeighborelement2() const;
    int getId() const;
    int getFiberId() const;
    const Eigen::Vector3d& getTangentialforce() const { return tangentialforce; }
    const Eigen::Vector3d& getTwisttorque() const { return twisttorque; }
    const Eigen::Vector3d& getBendtorque() const { return bendtorque; }
    double getEnergydissipation() const {return energyDissipation;} 
    double getPotentialenergy(const std::unique_ptr<SphereParticle> &sphere1, const std::unique_ptr<SphereParticle> &sphere2);
    double getkineticenergy();
    const PropertyTypeID &getType() const;

    void updateBond(std::unique_ptr<SphereParticle> &sphere1, std::unique_ptr<SphereParticle> &sphere2, double timeStep);

    std::string save_tostring() const;

    std::shared_ptr<ParticlePropertyManager>  getParticlePropertyManager() const { return manager;}

    static void computeOverlap(const Eigen::Vector3d &node1position, const Eigen::Vector3d &node2position,
                        const Eigen::Vector3d &sphereposition, double &t,Eigen::Vector3d &projection);
    static void computeOverlap(const Eigen::Vector3d &thisnode1position, const Eigen::Vector3d &thisnode2position, double& s,
                        const Eigen::Vector3d &anothernode1position, const Eigen::Vector3d &anothernode2position, double& t);
    static void computeOverlap(const Eigen::Vector3d &node1position, const Eigen::Vector3d &node2position,
                        const std::unique_ptr<PlaneWall> &planewall, double &distance);

private:
    int id;
    PropertyTypeID type;
    std::shared_ptr<ParticlePropertyManager> manager;
    int fiberid;
    int node1;
    int node2;

    int neighborelement1;
    int neighborelement2;

    Eigen::Vector3d tangentialforce;
    Eigen::Vector3d twisttorque;
    Eigen::Vector3d bendtorque;

    double energyDissipation;

    static double GetClampedRoot(double const &slope, double const &h0, double const &h1);
    // Compute the intersection of the line dR/ds = 0 with the domain
    // [0,1]^2. The direction of the line dR/ds is conjugate to (1,0),
    // so the algorithm for minimization is effectively the conjugate
    // gradient algorithm for a quadratic function.
    static void ComputeIntersection(std::array<double, 2> const &sValue,
                             std::array<int32_t, 2> const &classify, double const &b, double const &f00,
                             double const &f10, std::array<int32_t, 2> &edge,
                             std::array<std::array<double, 2>, 2> &end);
    // Compute the location of the minimum of R on the segment of
    // intersection for the line dR/ds = 0 and the domain [0,1]^2.
    static void ComputeMinimumParameters(std::array<int32_t, 2> const &edge,
                                  std::array<std::array<double, 2>, 2> const &end, double const &b, double const &c,
                                  double const &e, double const &g00, double const &g10, double const &g01, double const &g11,
                                  double &s, double &t);
};