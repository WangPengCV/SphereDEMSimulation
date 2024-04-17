#pragma once
#include "ParticlePropertyManager.h"
#include "PlaneWall.h"
#include "SphereParticle.h"
#include "SphereCylinderBond.h"
#include "Fiber.h"
#include "GridBasedContactDetection.h"
#include "ContactForce.h"

#include "SpringOscillator.h"
#include "ConstantMove.h"

#include "CylinderContainer.h"
#include "ConnectedGeometry.h"

#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <random>
#include <iomanip>

class DEMProperties
{

public:
    DEMProperties();

    void loadFromFile(const std::string &filename);

    void saveToFile(const std::string &filename) const;

    std::shared_ptr<ParticlePropertyManager> getParticleManager() const
    {
        return manager;
    }

    std::shared_ptr<ContactForce> getContactForce() const
    {
        return contactforce;
    }

    std::shared_ptr<GridBasedContactDetection> getGridBasedContactDetection() const
    {
        return gridbasedcontactdetection;
    }

    const std::unordered_map<int, std::unique_ptr<PlaneWall>> &getPlaneWall() const
    {
        return planewalls;
    }
    void setplanewalls(std::unordered_map<int, std::unique_ptr<PlaneWall>> &planeWalls);
    void setplanewall(std::unique_ptr<PlaneWall> &planeWall);

    const std::unordered_map<int, std::unique_ptr<CylinderContainer>> &getCylinderWall() const
    {
        return cylinderwalls;
    }
    void setcylinderwalls(std::unordered_map<int, std::unique_ptr<CylinderContainer>> &cylinderWalls);
    void setcylinderwall(std::unique_ptr<CylinderContainer> &cylinderWall);

    const std::unordered_map<int, std::unique_ptr<SphereParticle>> &getsphereParticles() const
    {
        return sphereparticles;
    }
    void setsphereParticles(std::unordered_map<int, std::unique_ptr<SphereParticle>> &sphereParticles);
    void removesphereParticles(int id);
    void addsphereParticles(std::unique_ptr<SphereParticle> &sp);

    const std::unordered_map<int, std::unique_ptr<SphereParticle>> &getfibersphereParticles() const
    {
        return fibershpereparticles;
    }
    void setfibersphereParticles(std::unordered_map<int, std::unique_ptr<SphereParticle>> &fiberSphereParticles);
    void removefibersphereParticles(int id);
    void addfibersphereParticles(std::unique_ptr<SphereParticle> &sp);

    const std::unordered_map<int, std::unique_ptr<SphereCylinderBond>> &getfiberbonds() const
    {
        return fiberbonds;
    }
    void setfiberbonds(std::unordered_map<int, std::unique_ptr<SphereCylinderBond>> &fiberBonds);
    void removefiberbonds(int id);
    void addfiberbonds(std::unique_ptr<SphereCylinderBond> &sp);

    double getTotalTime() const
    {
        return totalTime;
    }
    void setTotalTime(double);

    double getCurrentTime() const
    {
        return currentTime;
    }
    void setCurrentTime(double currentTime);

    double getTimestep() const
    {
        return timestep;
    }
    void setTimestep(double);

    double getTaskTimestep() const
    {
        return taskTimestep;
    }
    void setTaskTimestep(double);

    int getShowInterval() const
    {
        return showInterval;
    }
    void setShowInterval(int);

    int getTaskShowInterval() const
    {
        return taskShowInterval;
    }
    void setTaskShowInterval(double);

    double getCriticalSpeed() const
    {
        return criticalSpeed;
    }
    void setCriticalSpeed(double);

    const std::map<int, Eigen::Vector3d> &getSpecificForces() const
    {
        return specificforces;
    }

    const std::unordered_map<int, std::unique_ptr<ConnectedGeometry>> &getConnectedgeometrys() const
    {
        return connectedgeometrys;
    }

    const Eigen::Vector3d &getGravity() const
    {
        return gravity;
    }
    void setGravity(const Eigen::Vector3d &);

    const Eigen::Vector3d &getGlobalForces() const
    {
        return globalforces;
    }
    void setGlobalForces(const Eigen::Vector3d &);

    const Eigen::Vector3d &getSimulationDimensions() const
    {
        return simulationdimensions;
    }
    void setSimulationDimensions(const Eigen::Vector3d &);

    bool getGernerateSphereFlag() const
    {
        return gernerateSphereFlag;
    }
    void setGernerateSphereFlag(bool);

    bool getgernerateFiberFlag() const
    {
        return gernerateFiberFlag;
    }
    void setGernerateFiberFlag(bool);

    double getAverageVelocity();
    bool isGenerateComplete() const { return (gernerateSphereFlag && gernerateFiberFlag); }
    // function for DEM
    void initialSimulation();
    void generateRemainingParticles();

    void handleCollisions();
    void bondforce();
    void applyExternalForces();
    void motion();

    void initial_task();
    // task function
    void dotask();

private:
    bool parseVector3d(std::istringstream &iss, Eigen::Vector3d &vec);

    void parseLine(const std::string &line);

    void line_process(std::string &line);

    void parseSphereProperties(std::istringstream &iss);

    void parsePlanewallProperties(std::istringstream &iss);

    void parseCylinderwallProperties(std::istringstream &iss);

    void parseFiberProperties(std::istringstream &iss);

    void parsePlaneWall(std::istringstream &iss);

    void parseConnectedGeometry(std::istringstream &iss);

    void parseCylinderWall(std::istringstream &iss);

    void parseRandomParticle(std::istringstream &iss);

    void parseSpecificParticle(std::istringstream &iss);

    void parseFiberParticle(std::istringstream &iss);

    std::shared_ptr<ParticlePropertyManager> manager;

    std::shared_ptr<ContactForce> contactforce;

    std::shared_ptr<GridBasedContactDetection> gridbasedcontactdetection;

    std::unordered_map<int, std::unique_ptr<PlaneWall>> planewalls;

    std::unordered_map<int, std::unique_ptr<ConnectedGeometry>> connectedgeometrys;

    std::unordered_map<int, std::unique_ptr<CylinderContainer>> cylinderwalls;

    std::unordered_map<int, std::unique_ptr<SphereParticle>> sphereparticles;

    std::vector<std::unique_ptr<Fiber>> fibers;

    std::unordered_map<int, std::unique_ptr<SphereCylinderBond>> fiberbonds;

    std::unordered_map<int, std::unique_ptr<SphereParticle>> fibershpereparticles;

    double timestep;

    double taskTimestep;

    double currentTime;

    double totalTime;

    int showInterval;

    double criticalSpeed;

    int taskShowInterval;

    Eigen::Vector3d gravity;

    std::map<int, Eigen::Vector3d> specificforces;

    Eigen::Vector3d globalforces;

    Eigen::Vector3d simulationdimensions;

    bool gernerateSphereFlag;

    bool gernerateFiberFlag;

    std::map<PropertyTypeID, std::string> gernerateInfor;

    std::map<int, std::shared_ptr<DynamicBoundaryMover>> dynamicBoundaryTask;

    double contactEnergydissipation;

    double bondEnergydissipation;
};
