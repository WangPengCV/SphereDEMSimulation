#include <mpi.h>
#include "DEMModel.h"
#include <filesystem>
#include <iostream>
#include <random>
#include <chrono>
class ParallelDEMModel
{

public:
    long long extractNumberFromFilename(const std::string &filename)
    {
        std::regex regex("DEMProperties(\\d+)\\.dat");
        std::smatch match;
        if (std::regex_search(filename, match, regex) && match.size() > 1)
        {

            return std::stoll(match.str(1));
        }
        return 0;
    }
    // Function to check if the simulation can resume from an existing state
    bool canResumeFromState(const std::string &stateFolder)
    {
        std::filesystem::path dirPath(stateFolder);
        return std::filesystem::exists(dirPath) && !std::filesystem::is_empty(dirPath);
    }

    // Function to find the last state file
    std::string findLastStateFile(const std::string &stateFolder)
    {
        std::filesystem::path dirPath(stateFolder);
        long long maxNumber = 0;
        std::string maxFile;
        for (const auto &entry : std::filesystem::directory_iterator(dirPath))
        {
            const auto &path = entry.path();
            if (entry.is_regular_file())
            {
                long long currentNumber = extractNumberFromFilename(path.filename().string());
                if (currentNumber >= maxNumber)
                {
                    maxNumber = currentNumber;
                    maxFile = path.filename().string();
                }
            }
        }

        if (!maxFile.empty())
        {
            return (dirPath / maxFile).string();
        }
        else
        {
            return "";
        }
    }

    struct SphereParticleData
    {
        int id;
        int category;
        int subType;
        int state;
        double position[3];
        double velocity[3];
        double omega[3];
        double force[3];
        double torque[3];
        SphereParticleData() = default;
    };
    struct FiberBondData
    {
        int id;
        int category;
        int subType;
        int fiberid;
        int node1;
        int node2;
        int neighborelement1;
        int neighborelement2;
        double energyDissipation;
        double tangentialforce[3];
        double twisttorque[3];
        double bendtorque[3];
    };
    struct SimulationParameters
    {
        double timestep;

        double taskTimestep;

        double currentTime;

        double totalTime;

        int showInterval;

        double criticalSpeed;

        int taskShowInterval;

        int iternum;

        bool gernerateSphereFlag;

        bool gernerateFiberFlag;

        double gravity[3];

        double globalforces[3];

        double simulationdimensions[3];
    };
    struct SphereParticlePropertiesData
    {
        int category;
        int subType;
        double density;
        double rolling_friction_coefficient;
        double slide_friction_coefficient;
        double Young_modulus;
        double restitution;
        double poisson_ratio;
        double momentofinertia;
        double radius;
        double mass;
        SphereParticlePropertiesData() = default;
    };
    struct FiberPropertiesData
    {
        int category;
        int subType;
        double density;
        double rolling_friction_coefficient;
        double slide_friction_coefficient;
        double Young_modulus;
        double restitution;
        double poisson_ratio;
        double aspectratio;
        double nodemomentofinertia;
        double radius;
        double nodemass;
        double elementlength;
        double normalmodulus;
        double shearmodulus;
        double twistmodulus;
        double bendingmodulus;
        int nodenumber;
        double bonddampingcoefficient;
        FiberPropertiesData() = default;
    };
    struct PlanewallPropertiesData
    {
        int category;
        int subType;
        double density;
        double rolling_friction_coefficient;
        double slide_friction_coefficient;
        double Young_modulus;
        double restitution;
        double poisson_ratio;
        double thickness;
        PlanewallPropertiesData() = default;
    };

    struct CylinderwallPropertiesData
    {
        int category;
        int subType;
        double density;
        double rolling_friction_coefficient;
        double slide_friction_coefficient;
        double Young_modulus;
        double restitution;
        double poisson_ratio;
        double thickness;
        bool hollow;
        CylinderwallPropertiesData() = default;
    };
    struct ContactInformationData
    {
        int particleId1;
        int particleId2;
        double tangentialDisplacement[3];
        double previousTangentialVelocity[3];
        double normalForce[3];
        double tangentialForce[3];
        ContactInformationData() = default;
    };
    struct PlanewallData
    {
        int id;
        int category;
        int subType;
        int state;
        double normal[3];
        double corner1[3];
        double corner2[3];
        double corner3[3];
        double velocity[3];
        double force[3];
        PlanewallData() = default;
    };

    struct CylinderwallData
    {
        int id;
        int category;
        int subType;
        int state;
        double radius;
        double startpoint[3];
        double endpoint[3];
        double velocity[3];
        double force[3];
        CylinderwallData() = default;
    };
    ParallelDEMModel(const std::string &filename);
    void runSimulation();
    void finalizeMPIType();

private:
    void createMPITypeForSimulationParameters();

    void createMPITypeForSphereParticleData();
    void createMPITypeForFiberBondData();

    void createMPITypeForSphereparticleProperties();
    void createMPITypeForFiberProperties();
    void createMPITypeForPlanewallProperties();
    void createMPITypeForCylinderwallProperties();

    void createMPITypeForContactInformation();

    void createMPITypeForPlanewall();
    void createMPITypeForCylinderwall();

    void distributeParticlesToProcesses();
    void collectParticlesFromProcesses();

    void distributePlanewallToProcesses(int id);
    void collectPlanewallfromProcesses(int id);


    void distributeCylinderwallToProcesses(int id);
    // void distributeSphereToProcesses(int id);
    // void distributeFiberbondToProcesses(int id);


    void neighborCommunication();
    void handleCollisions();
    void exchangeSharednodes();
    void motion();


    std::shared_ptr<DEMProperties> DEMproperties;
    std::shared_ptr<Visualization> vis;
    int world_rank, world_size;
    MPI_Datatype mpi_sphereparticleDataType, mpi_simulationParamsType, mpi_fiberbondDataType,
        mpi_sphereparticlepropertiesDataType, mpi_fiberpropertiesDataType, mpi_planewallpropertiesDataType,
        mpi_contactinformationDataType, mpi_planewallDataType, mpi_cylinderwallpropertiesDataType, mpi_cylinderwallDataType;

    SimulationParameters simParams;
    std::vector<int> offset;
    int numberOfGridX;
    int numberOfGridY;
    int numberOfGridZ;

    double gridSizeX;
    double gridSizeY;
    double gridSizeZ;

    std::unordered_set<int> neighborRankId;
    std::unordered_map<int, std::unordered_set<int>> ghostLayer;
    std::vector<std::vector<int>> sphereGrids;
    std::vector<std::vector<int>> fiberGrids;

    std::unordered_map<int, std::unique_ptr<SphereParticle>> ghostsphereparticles;

    std::unordered_map<int, std::unique_ptr<SphereCylinderBond>> ghostfiberbonds;

    std::unordered_map<int, std::unique_ptr<SphereParticle>> ghostfibershpereparticles;

    std::unordered_map<int, std::unordered_set<int>> neighborbonds;

};