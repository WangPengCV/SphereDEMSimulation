#pragma once
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <SphereParticle.h>
#include <SphereCylinderBond.h>
#include <PlaneWall.h>
#include <CylinderContainer.h>
#include <RectangularContainer.h>
class GridBasedContactDetection
{
public:
    GridBasedContactDetection();

    void initial(double domain_x, double domain_y, double domain_z, double gridsize);

    void ParticleBroadPhase(const std::unordered_map<int, std::unique_ptr<SphereParticle>> &particles,
                            const std::unordered_map<int, std::unique_ptr<SphereCylinderBond>> &fiberbonds,
                            const std::unordered_map<int, std::unique_ptr<SphereParticle>> &fibersphereparticles,
                            std::unordered_map<int, std::unordered_set<int>> &SScontactparis,
                            std::unordered_map<int, std::unordered_set<int>> &SFcontactparis,
                            std::unordered_map<int, std::unordered_set<int>> &FFcontactparis);

    void planewallBroadPhase(const std::unordered_map<int, std::unique_ptr<PlaneWall>> &planewalls,
                             std::unordered_map<int, std::unordered_set<int>> &PScontactparis,
                             std::unordered_map<int, std::unordered_set<int>> &PFcontactpairs);
    void RectangularContainerBroadPhase(const std::unique_ptr<RectangularContainer> &rectangularcontainer,
                                        std::unordered_map<int, std::unordered_set<int>> &PScontactparis,
                                        std::unordered_map<int, std::unordered_set<int>> &PFcontactpairs);
    void cylinderwallBroadPhase(const std::unordered_map<int, std::unique_ptr<CylinderContainer>> &CylinderContainer,
                                std::unordered_map<int, std::unordered_set<int>> &CScontactparis,
                                std::unordered_map<int, std::unordered_set<int>> &CFcontactpairs);
    double getDomainX() const
    {
        return domainX;
    }

    double getDomainY() const
    {
        return domainY;
    }

    double getDomainZ() const
    {
        return domainZ;
    }

    double getGridSizeX() const
    {
        return gridSizeX;
    }

    double getGridSizeY() const
    {
        return gridSizeY;
    }

    double getGridSizeZ() const
    {
        return gridSizeZ;
    }

    int getNumberOfGridX() const
    {
        return numberOfGridX;
    }

    int getNumberOfGridY() const
    {
        return numberOfGridY;
    }

    int getNumberOfGridZ() const
    {
        return numberOfGridZ;
    }
    const std::unordered_map<int, std::unordered_set<int>> &getGridSaveSphereParticles() const
    {
        return gridSaveSphereParticles;
    }
    const std::unordered_map<int, std::unordered_set<int>> &getGridSaveFiberBonds() const
    {
        return gridSaveFiberBonds;
    }
    void assignParticleToGrid(const std::unordered_map<int, std::unique_ptr<SphereParticle>> &sphereparticles,
                              const std::unordered_map<int, std::unique_ptr<SphereCylinderBond>> &fiberbonds,
                              const std::unordered_map<int, std::unique_ptr<SphereParticle>> &fibersphereparticles);
    void ghostParticleBroadPhase(const std::unordered_map<int, std::unique_ptr<SphereParticle>> &ghostsphereparticles,
                                 const std::unordered_map<int, std::unique_ptr<SphereCylinderBond>> &ghostfiberbonds,
                                 const std::unordered_map<int, std::unique_ptr<SphereParticle>> &ghostfibershpereparticles,
                                 const std::unordered_map<int, std::unique_ptr<SphereCylinderBond>> &fiberbonds,
                                 std::unordered_map<int, std::unordered_set<int>> &SScontactparis,
                                 std::unordered_map<int, std::unordered_set<int>> &SFcontactparis,
                                 std::unordered_map<int, std::unordered_set<int>> &FFcontactparis,
                                 std::unordered_map<int, std::unordered_set<int>> &neighborbonds);

private:
    int getGridIndex(double x, double y, double z);

    double domainX;
    double domainY;
    double domainZ;

    double gridSizeX;
    double gridSizeY;
    double gridSizeZ;
    double gridsize;
    int numberOfGridX;
    int numberOfGridY;
    int numberOfGridZ;

    std::unordered_map<int, std::unordered_set<int>> gridSaveSphereParticles;

    std::unordered_map<int, std::unordered_set<int>> gridSaveFiberBonds;

    std::unordered_map<int, std::unordered_set<int>> planeSaveGrid;

    std::unordered_map<int, std::unordered_set<int>> rectangularcontainerSaveGrid;

    std::unordered_map<int, std::unordered_set<int>> cylindercontainerSaveGrid;
};
