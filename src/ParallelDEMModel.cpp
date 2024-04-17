#include "ParallelDEMModel.h"
#include <mpi.h>
#include <filesystem>
#include <iostream>

ParallelDEMModel::ParallelDEMModel(const std::string &filename)
{

    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    createMPITypeForSphereParticleData();
    createMPITypeForFiberBondData();

    createMPITypeForSimulationParameters();

    createMPITypeForSphereparticleProperties();
    createMPITypeForFiberProperties();
    createMPITypeForPlanewallProperties();
    createMPITypeForCylinderwallProperties();

    createMPITypeForPlanewall();
    createMPITypeForCylinderwall();

    createMPITypeForContactInformation();

    offset.resize(world_size * 2, 0);

    std::vector<SphereParticlePropertiesData> sphereparticleproperties;
    std::vector<FiberPropertiesData> fiberproperties;
    std::vector<PlanewallPropertiesData> planewallproperties;
    std::vector<CylinderwallPropertiesData> cylinderwallproperties;

    std::vector<PlanewallData> planewallData;
    std::vector<CylinderwallData> cylinderwallData;

    int sphereparticlepronumber = 0;
    int fiberpronumber = 0;
    int planewallpronumber = 0;
    int cylinderwallpronumber = 0;

    int planewallnumber = 0;
    int cylinderwallnumber = 0;

    // load simulation files by main world rank
    if (world_rank == 0)
    {
        DEMproperties = std::make_shared<DEMProperties>();
        std::string stateFolder = "DEMProperties";
        if (canResumeFromState(stateFolder))
        {
            std::string lastStateFile = findLastStateFile(stateFolder);
            std::cout << "Resuming from state: " << lastStateFile << std::endl;
            DEMproperties->loadFromFile(lastStateFile); // Make sure this function loads all needed state
            simParams.iternum = extractNumberFromFilename(lastStateFile);
        }
        else
        {

            DEMproperties->loadFromFile(filename);

            simParams.iternum = -1;
        }
        vis = std::make_shared<Visualization>(*DEMproperties);

        simParams.timestep = DEMproperties->getTimestep();
        simParams.taskTimestep = DEMproperties->getTaskTimestep();
        simParams.currentTime = DEMproperties->getCurrentTime();
        simParams.totalTime = DEMproperties->getTotalTime();
        simParams.criticalSpeed = DEMproperties->getCriticalSpeed();
        simParams.showInterval = DEMproperties->getShowInterval();
        simParams.taskShowInterval = DEMproperties->getTaskShowInterval();
        simParams.gernerateSphereFlag = DEMproperties->getGernerateSphereFlag();
        simParams.gernerateFiberFlag = DEMproperties->getgernerateFiberFlag();

        for (int i = 0; i < 3; ++i)
        {
            simParams.gravity[i] = DEMproperties->getGravity()[i];
            simParams.globalforces[i] = DEMproperties->getGlobalForces()[i];
            simParams.simulationdimensions[i] = DEMproperties->getSimulationDimensions()[i];
        }
        for (const auto &particleproperties : DEMproperties->getParticleManager()->getParticleProperties())
        {
            if (auto &sphereProps = std::dynamic_pointer_cast<SphereProperties>(particleproperties.second))
            {

                SphereParticlePropertiesData spproperties;
                spproperties.category = particleproperties.first.getCategory();
                spproperties.subType = particleproperties.first.getSubType();
                spproperties.density = sphereProps->getDensity();
                spproperties.rolling_friction_coefficient = sphereProps->getRollingFriction();
                spproperties.slide_friction_coefficient = sphereProps->getSlidingFriction();
                spproperties.Young_modulus = sphereProps->getYoungModulus();
                spproperties.restitution = sphereProps->getRestitution();
                spproperties.poisson_ratio = sphereProps->getPoissonRatio();
                spproperties.momentofinertia = sphereProps->getMomentOfInertia();
                spproperties.radius = sphereProps->getRadius();
                spproperties.mass = sphereProps->getMass();
                sphereparticleproperties.push_back(spproperties);
            }
            else if (auto &fiberProps = std::dynamic_pointer_cast<FiberProperties>(particleproperties.second))
            {
                FiberPropertiesData FProperties;
                FProperties.category = particleproperties.first.getCategory();
                FProperties.subType = particleproperties.first.getSubType();
                FProperties.density = fiberProps->getDensity();
                FProperties.rolling_friction_coefficient = fiberProps->getRollingFriction();
                FProperties.slide_friction_coefficient = fiberProps->getSlidingFriction();
                FProperties.Young_modulus = fiberProps->getYoungModulus();
                FProperties.restitution = fiberProps->getRestitution();
                FProperties.poisson_ratio = fiberProps->getPoissonRatio();
                FProperties.aspectratio = fiberProps->getAspectratio();
                FProperties.nodemomentofinertia = fiberProps->getNodemomentofinertia();
                FProperties.radius = fiberProps->getRadius();
                FProperties.nodemass = fiberProps->getNodemass();
                FProperties.elementlength = fiberProps->getElementlength();
                FProperties.normalmodulus = fiberProps->getNormalmodulus();
                FProperties.shearmodulus = fiberProps->getShearmodulus();
                FProperties.twistmodulus = fiberProps->getTwistmodulus();
                FProperties.bendingmodulus = fiberProps->getBendingmodulus();
                FProperties.nodenumber = fiberProps->getNodenumber();
                FProperties.bonddampingcoefficient = fiberProps->getBonddampingcoefficient();
                fiberproperties.push_back(FProperties);
            }
            else if (auto &planewallProps = std::dynamic_pointer_cast<PlanewallProperties>(particleproperties.second))
            {
                PlanewallPropertiesData properties;
                properties.category = particleproperties.first.getCategory();
                properties.subType = particleproperties.first.getSubType();
                properties.density = planewallProps->getDensity();
                properties.rolling_friction_coefficient = planewallProps->getRollingFriction();
                properties.slide_friction_coefficient = planewallProps->getSlidingFriction();
                properties.Young_modulus = planewallProps->getYoungModulus();
                properties.restitution = planewallProps->getRestitution();
                properties.poisson_ratio = planewallProps->getPoissonRatio();
                properties.thickness = planewallProps->getThickness();
                planewallproperties.push_back(properties);
            }
            else if (auto &cylinderwallProps = std::dynamic_pointer_cast<CylinderwallProperties>(particleproperties.second))
            {
                CylinderwallPropertiesData properties;
                properties.category = particleproperties.first.getCategory();
                properties.subType = particleproperties.first.getSubType();
                properties.density = cylinderwallProps->getDensity();
                properties.rolling_friction_coefficient = cylinderwallProps->getRollingFriction();
                properties.slide_friction_coefficient = cylinderwallProps->getSlidingFriction();
                properties.Young_modulus = cylinderwallProps->getYoungModulus();
                properties.restitution = cylinderwallProps->getRestitution();
                properties.poisson_ratio = cylinderwallProps->getPoissonRatio();
                properties.thickness = cylinderwallProps->getThickness();
                properties.hollow = cylinderwallProps->getHollow();
                cylinderwallproperties.push_back(properties);
            }
        }

        for (const auto &planewall : DEMproperties->getPlaneWall())
        {
            PlanewallData pw;
            pw.id = planewall.second->getId();
            pw.category = planewall.second->getType().getCategory();
            pw.subType = planewall.second->getType().getSubType();
            pw.state = planewall.second->getState();
            pw.normal[0] = planewall.second->getNormal()[0];
            pw.normal[1] = planewall.second->getNormal()[1];
            pw.normal[2] = planewall.second->getNormal()[2];

            pw.corner1[0] = planewall.second->getCorner1()[0];
            pw.corner1[1] = planewall.second->getCorner1()[1];
            pw.corner1[2] = planewall.second->getCorner1()[2];

            pw.corner2[0] = planewall.second->getCorner2()[0];
            pw.corner2[1] = planewall.second->getCorner2()[1];
            pw.corner2[2] = planewall.second->getCorner2()[2];

            pw.corner3[0] = planewall.second->getCorner3()[0];
            pw.corner3[1] = planewall.second->getCorner3()[1];
            pw.corner3[2] = planewall.second->getCorner3()[2];

            pw.velocity[0] = planewall.second->getVelocity()[0];
            pw.velocity[1] = planewall.second->getVelocity()[1];
            pw.velocity[2] = planewall.second->getVelocity()[2];

            pw.force[0] = planewall.second->getForce()[0];
            pw.force[1] = planewall.second->getForce()[1];
            pw.force[2] = planewall.second->getForce()[2];
            planewallData.push_back(pw);
        }
        for (const auto &cylinderwall : DEMproperties->getCylinderWall())
        {
            CylinderwallData cw;
            cw.id = cylinderwall.second->getId();
            cw.category = cylinderwall.second->getType().getCategory();
            cw.subType = cylinderwall.second->getType().getSubType();
            cw.state = cylinderwall.second->getState();
            cw.radius = cylinderwall.second->getRadius();

            cw.startpoint[0] = cylinderwall.second->getStartpoint()[0];
            cw.startpoint[1] = cylinderwall.second->getStartpoint()[1];
            cw.startpoint[2] = cylinderwall.second->getStartpoint()[2];

            cw.endpoint[0] = cylinderwall.second->getEndpoint()[0];
            cw.endpoint[1] = cylinderwall.second->getEndpoint()[1];
            cw.endpoint[2] = cylinderwall.second->getEndpoint()[2];

            cw.velocity[0] = cylinderwall.second->getVelocity()[0];
            cw.velocity[1] = cylinderwall.second->getVelocity()[1];
            cw.velocity[2] = cylinderwall.second->getVelocity()[2];

            cw.force[0] = cylinderwall.second->getForce()[0];
            cw.force[1] = cylinderwall.second->getForce()[1];
            cw.force[2] = cylinderwall.second->getForce()[2];
            cylinderwallData.push_back(cw);
        }
        sphereparticlepronumber = static_cast<int>(sphereparticleproperties.size());
        fiberpronumber = static_cast<int>(fiberproperties.size());
        planewallpronumber = static_cast<int>(planewallproperties.size());
        cylinderwallpronumber = static_cast<int>(cylinderwallproperties.size());
        planewallnumber = static_cast<int>(planewallData.size());
        cylinderwallnumber = static_cast<int>(cylinderwallData.size());
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // transform simulation parameters to anther process
    MPI_Bcast(&simParams, 1, mpi_simulationParamsType, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sphereparticlepronumber, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&fiberpronumber, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&planewallpronumber, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&cylinderwallpronumber, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&planewallnumber, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&cylinderwallnumber, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (world_rank != 0)
    {
        // Initialize DEMProperties for non-zero ranks
        DEMproperties = std::make_shared<DEMProperties>();
        DEMproperties->setTimestep(simParams.timestep);
        DEMproperties->setTaskTimestep(simParams.taskTimestep);
        DEMproperties->setCurrentTime(simParams.currentTime);
        DEMproperties->setTotalTime(simParams.totalTime);
        DEMproperties->setCriticalSpeed(simParams.criticalSpeed);
        DEMproperties->setShowInterval(simParams.showInterval);
        DEMproperties->setTaskShowInterval(simParams.taskShowInterval);
        DEMproperties->setGravity(Eigen::Vector3d(simParams.gravity[0], simParams.gravity[1], simParams.gravity[2]));
        DEMproperties->setGlobalForces(Eigen::Vector3d(simParams.globalforces[0], simParams.globalforces[1], simParams.globalforces[2]));
        DEMproperties->setSimulationDimensions(Eigen::Vector3d(simParams.simulationdimensions[0], simParams.simulationdimensions[1], simParams.simulationdimensions[2]));
        DEMproperties->setGernerateSphereFlag(simParams.gernerateSphereFlag);
        DEMproperties->setGernerateFiberFlag(simParams.gernerateFiberFlag);

        sphereparticleproperties.resize(sphereparticlepronumber);
        fiberproperties.resize(fiberpronumber);
        planewallproperties.resize(planewallpronumber);
        cylinderwallproperties.resize(cylinderwallpronumber);
        planewallData.resize(planewallnumber);
        cylinderwallData.resize(cylinderwallnumber);
    }
    if (sphereparticlepronumber > 0)
    {
        MPI_Bcast(sphereparticleproperties.data(), sphereparticlepronumber, mpi_sphereparticlepropertiesDataType, 0, MPI_COMM_WORLD);
    }
    if (fiberpronumber > 0)
    {
        MPI_Bcast(fiberproperties.data(), fiberpronumber, mpi_fiberpropertiesDataType, 0, MPI_COMM_WORLD);
    }
    if (planewallpronumber > 0)
    {
        MPI_Bcast(planewallproperties.data(), planewallpronumber, mpi_planewallpropertiesDataType, 0, MPI_COMM_WORLD);
    }
    if (cylinderwallpronumber > 0)
    {
        MPI_Bcast(cylinderwallproperties.data(), cylinderwallpronumber, mpi_cylinderwallpropertiesDataType, 0, MPI_COMM_WORLD);
    }
    if (planewallnumber > 0)
    {
        MPI_Bcast(planewallData.data(), planewallnumber, mpi_planewallDataType, 0, MPI_COMM_WORLD);
    }
    if (cylinderwallnumber > 0)
    {
        MPI_Bcast(cylinderwallData.data(), cylinderwallnumber, mpi_cylinderwallDataType, 0, MPI_COMM_WORLD);
    }

    if (world_rank != 0)
    {
        for (int i = 0; i < fiberproperties.size(); ++i)
        {
            auto properties = std::make_shared<FiberProperties>(fiberproperties[i].density, fiberproperties[i].rolling_friction_coefficient, fiberproperties[i].slide_friction_coefficient,
                                                                fiberproperties[i].Young_modulus, fiberproperties[i].restitution, fiberproperties[i].poisson_ratio, fiberproperties[i].radius, fiberproperties[i].elementlength, fiberproperties[i].normalmodulus,
                                                                fiberproperties[i].shearmodulus, fiberproperties[i].twistmodulus, fiberproperties[i].bendingmodulus, fiberproperties[i].nodenumber, fiberproperties[i].aspectratio, fiberproperties[i].bonddampingcoefficient);
            DEMproperties->getParticleManager()->addFiberProperties(PropertyTypeID(fiberproperties[i].category, fiberproperties[i].subType), properties);
        }
        for (int i = 0; i < sphereparticleproperties.size(); ++i)
        {
            auto properties = std::make_shared<SphereProperties>(sphereparticleproperties[i].density, sphereparticleproperties[i].mass, sphereparticleproperties[i].radius, sphereparticleproperties[i].rolling_friction_coefficient,
                                                                 sphereparticleproperties[i].slide_friction_coefficient,
                                                                 sphereparticleproperties[i].Young_modulus, sphereparticleproperties[i].restitution,
                                                                 sphereparticleproperties[i].poisson_ratio, sphereparticleproperties[i].momentofinertia);
            DEMproperties->getParticleManager()->addSphereProperties(PropertyTypeID(sphereparticleproperties[i].category, sphereparticleproperties[i].subType), properties);
        }
        for (int i = 0; i < planewallproperties.size(); ++i)
        {
            auto properties = std::make_shared<PlanewallProperties>(planewallproperties[i].density, planewallproperties[i].thickness, planewallproperties[i].rolling_friction_coefficient, planewallproperties[i].slide_friction_coefficient,
                                                                    planewallproperties[i].Young_modulus, planewallproperties[i].restitution, planewallproperties[i].poisson_ratio);
            DEMproperties->getParticleManager()->addPlanewallProperties(PropertyTypeID(planewallproperties[i].category, planewallproperties[i].subType), properties);
        }
        for (int i = 0; i < cylinderwallproperties.size(); ++i)
        {
            auto properties = std::make_shared<CylinderwallProperties>(cylinderwallproperties[i].density, cylinderwallproperties[i].thickness, cylinderwallproperties[i].rolling_friction_coefficient, cylinderwallproperties[i].slide_friction_coefficient,
                                                                       cylinderwallproperties[i].Young_modulus, cylinderwallproperties[i].restitution, cylinderwallproperties[i].poisson_ratio, cylinderwallproperties[i].hollow);
            DEMproperties->getParticleManager()->addCylinerwallProperties(PropertyTypeID(cylinderwallproperties[i].category, cylinderwallproperties[i].subType), properties);
        }
        std::unordered_map<int, std::unique_ptr<PlaneWall>> planewalls;
        for (const auto &data : planewallData)
        {
            Eigen::Vector3d normal(data.normal[0], data.normal[1], data.normal[2]);
            Eigen::Vector3d corner1(data.corner1[0], data.corner1[1], data.corner1[2]);
            Eigen::Vector3d corner2(data.corner2[0], data.corner2[1], data.corner2[2]);
            Eigen::Vector3d corner3(data.corner3[0], data.corner3[1], data.corner3[2]);
            Eigen::Vector3d velocity(data.velocity[0], data.velocity[1], data.velocity[2]);
            Eigen::Vector3d force(data.force[0], data.force[1], data.force[2]);

            PropertyTypeID type = PropertyTypeID(data.category, data.subType);
            auto planewall = std::make_unique<PlaneWall>(data.id, type, data.state, normal, corner1, corner2, corner3, velocity, force);
            planewalls.insert({planewall->getId(), std::move(planewall)});
        }
        DEMproperties->setplanewalls(planewalls);

        std::unordered_map<int, std::unique_ptr<CylinderContainer>> cylinderwalls;
        for (const auto &data : cylinderwallData)
        {
            Eigen::Vector3d startpoint(data.startpoint[0], data.startpoint[1], data.startpoint[2]);
            Eigen::Vector3d endpoint(data.endpoint[0], data.endpoint[1], data.endpoint[2]);

            Eigen::Vector3d velocity(data.velocity[0], data.velocity[1], data.velocity[2]);
            Eigen::Vector3d force(data.force[0], data.force[1], data.force[2]);

            PropertyTypeID type = PropertyTypeID(data.category, data.subType);
            auto cylinderwall = std::make_unique<CylinderContainer>(data.id, type, data.state, data.radius, startpoint, endpoint, velocity, force);
            cylinderwalls.insert({cylinderwall->getId(), std::move(cylinderwall)});
        }
        DEMproperties->setcylinderwalls(cylinderwalls);
    }
    DEMproperties->initialSimulation();
    numberOfGridX = DEMproperties->getGridBasedContactDetection()->getNumberOfGridX();
    numberOfGridY = DEMproperties->getGridBasedContactDetection()->getNumberOfGridY();
    numberOfGridZ = DEMproperties->getGridBasedContactDetection()->getNumberOfGridZ();

    gridSizeX = DEMproperties->getGridBasedContactDetection()->getGridSizeX();
    gridSizeY = DEMproperties->getGridBasedContactDetection()->getGridSizeY();
    gridSizeZ = DEMproperties->getGridBasedContactDetection()->getGridSizeZ();

    distributeParticlesToProcesses();
}

void ParallelDEMModel::distributeParticlesToProcesses()
{
    // getRankGrid
    if (world_rank == 0)
    {
        int gridNumber = numberOfGridX * numberOfGridZ * numberOfGridY;
        sphereGrids.clear();
        fiberGrids.clear();
        sphereGrids.resize(gridNumber);
        fiberGrids.resize(gridNumber);
        auto &fibershpereparticles = DEMproperties->getfibersphereParticles();
        auto &sphereparticles = DEMproperties->getsphereParticles();
        auto &fiberbonds = DEMproperties->getfiberbonds();
        for (auto &fiberbond : fiberbonds)
        {
            auto &node1 = fibershpereparticles.at(fiberbond.second->getNode1());
            auto &node2 = fibershpereparticles.at(fiberbond.second->getNode2());
            double centerX = 0.5 * (node1->getPosition()[0] + node2->getPosition()[0]);
            double centerY = 0.5 * (node1->getPosition()[1] + node2->getPosition()[1]);
            double centerZ = 0.5 * (node1->getPosition()[2] + node2->getPosition()[2]);

            int cellX = static_cast<int>(centerX / gridSizeX);
            int cellY = static_cast<int>(centerY / gridSizeY);
            int cellZ = static_cast<int>(centerZ / gridSizeZ);

            int cellId = cellX + cellZ * numberOfGridX + cellY * numberOfGridX * numberOfGridZ;

            fiberGrids[cellId].push_back(fiberbond.second->getId());
        }
        for (auto &sphere : sphereparticles)
        {
            int cellX = static_cast<int>(sphere.second->getPosition()[0] / gridSizeX);
            int cellY = static_cast<int>(sphere.second->getPosition()[1] / gridSizeY);
            int cellZ = static_cast<int>(sphere.second->getPosition()[2] / gridSizeZ);

            sphereGrids[cellX + cellZ * numberOfGridX + cellY * numberOfGridX * numberOfGridZ].push_back(sphere.second->getId());
        }

        int averageParticleNumber = floor((fiberbonds.size() + sphereparticles.size()) / world_size);

        int count_number = 0;
        int rank_index = 0;

        for (int k = 0; k < numberOfGridY; ++k)
        {
            for (int j = 0; j < numberOfGridZ; ++j)
            {
                for (int i = 0; i < numberOfGridX; ++i)
                {
                    int gridIndex = i + j * numberOfGridX + k * numberOfGridX * numberOfGridZ;
                    count_number += static_cast<int>(sphereGrids[gridIndex].size());
                    count_number += static_cast<int>(fiberGrids[gridIndex].size());

                    if (0.99 * averageParticleNumber <= count_number || gridIndex == gridNumber - 1)
                    {
                        if (rank_index == world_size - 1)
                        {
                            i = numberOfGridX;
                            j = numberOfGridZ;
                            k = numberOfGridY;
                            offset[rank_index * 2 + 1] = gridNumber - 1;
                        }
                        else
                        {
                            offset[rank_index * 2 + 1] = gridIndex;

                            count_number = 0;
                            rank_index++;
                            offset[rank_index * 2] = gridIndex + 1;
                        }
                    }
                }
            }
        }
    }
    MPI_Bcast(offset.data(), world_size * 2, MPI_INT, 0, MPI_COMM_WORLD);
    // getGhostLayer
    neighborRankId.clear();
    ghostLayer.clear();

    for (int gridID = offset[world_rank * 2]; gridID <= offset[world_rank * 2 + 1]; ++gridID)
    {
        int gridy = gridID / (numberOfGridX * numberOfGridZ);
        int temp = gridID % (numberOfGridX * numberOfGridZ);
        int gridz = temp / numberOfGridX;
        int gridx = temp % numberOfGridX;
        for (int dy = -1; dy <= 1; dy++)
        {
            for (int dz = -1; dz <= 1; dz++)
            {
                for (int dx = -1; dx <= 1; dx++)
                {
                    if (dx == 0 && dy == 0 && dz == 0)
                        continue; // Skip the current grid itself

                    int xNeighbor = gridx + dx;
                    int zNeighbor = gridz + dz;
                    int yNeighbor = gridy + dy;

                    // Check for boundary conditions
                    if (xNeighbor >= 0 && xNeighbor < numberOfGridX &&
                        zNeighbor >= 0 && zNeighbor < numberOfGridZ &&
                        yNeighbor >= 0 && yNeighbor < numberOfGridY)
                    {
                        int neighborId = (yNeighbor * numberOfGridZ + zNeighbor) * numberOfGridX + xNeighbor;
                        for (int otherRank = 0; otherRank < world_size; ++otherRank)
                        {
                            if (otherRank == world_rank)
                                continue; // Skip the current rank itself
                            if (neighborId >= offset[otherRank * 2] && neighborId <= offset[otherRank * 2 + 1])
                            {
                                neighborRankId.insert(otherRank);
                                ghostLayer[otherRank].insert(gridID);
                            }
                        }
                    }
                }
            }
        }
    }
    // main rank
    std::vector<SphereParticleData> sendsphereData;
    std::vector<SphereParticleData> sendfibersphereData;
    std::vector<FiberBondData> sendfiberbondData;

    std::vector<ContactInformationData> spherespherecontactinformationData;
    std::vector<ContactInformationData> spherefibercontactinformationData;
    std::vector<ContactInformationData> sphereplanewallcontactinformationData;
    std::vector<ContactInformationData> spherecylinderwallcontactinformationData;

    std::vector<ContactInformationData> fiberfibercontactinformationData;
    std::vector<ContactInformationData> fiberplanewallcontactinformationData;
    std::vector<ContactInformationData> fibercylinderwallcontactinformationData;

    std::vector<int> sendsphereCounts(world_size, 0);
    std::vector<int> sendfibersphereCounts(world_size, 0);
    std::vector<int> sendfiberbondCounts(world_size, 0);

    std::vector<int> sendsphereDispls(world_size, 0);
    std::vector<int> sendfibersphereDispls(world_size, 0);
    std::vector<int> sendfiberbondDispls(world_size, 0);

    std::vector<int> spherespherecontactCounts(world_size, 0);
    std::vector<int> spherefibercontactCounts(world_size, 0);
    std::vector<int> sphereplanewallcontactCounts(world_size, 0);
    std::vector<int> spherecylinderwallcontactCounts(world_size, 0);

    std::vector<int> spherespherecontactDispls(world_size, 0);
    std::vector<int> spherefibercontactDispls(world_size, 0);
    std::vector<int> sphereplanewallcontactDispls(world_size, 0);
    std::vector<int> spherecylinderwallcontactDispls(world_size, 0);

    std::vector<int> fiberfibercontactCounts(world_size, 0);
    std::vector<int> fiberplanewallcontactCounts(world_size, 0);
    std::vector<int> fibercylinderwallcontactCounts(world_size, 0);

    std::vector<int> fiberfibercontactDispls(world_size, 0);
    std::vector<int> fiberplanewallcontactDispls(world_size, 0);
    std::vector<int> fibercylinderwallcontactDispls(world_size, 0);

    if (world_rank == 0)
    {
        auto &fibersphereparticles = DEMproperties->getfibersphereParticles();
        auto &sphereparticles = DEMproperties->getsphereParticles();
        auto &fiberbonds = DEMproperties->getfiberbonds();

        auto &spheresphereContactInformationList = DEMproperties->getContactForce()->getSphereSphereContactInformationList();
        auto &spherefiberContactInformationList = DEMproperties->getContactForce()->getSphereFiberContactInformationList();
        auto &sphereplanewallContactInformationList = DEMproperties->getContactForce()->getplanewallSphereContactInformationList();
        auto &spherecylinderwallContactInformationList = DEMproperties->getContactForce()->getcylinderwallSphereContactInformationList();

        auto &fiberfiberContactInformationList = DEMproperties->getContactForce()->getFiberFiberContactInformationList();
        auto &fiberplanewallContactInformationList = DEMproperties->getContactForce()->getplanewallFiberContactInformationList();
        auto &fibercylinderwallContactInformationList = DEMproperties->getContactForce()->getcylinderwallFiberContactInformationList();

        for (int rankId = 0; rankId < world_size; ++rankId)
        {
            int sendsphere = 0;
            int sendfibersphere = 0;
            int sendfiberbond = 0;

            int sendspheresphereContact = 0;
            int sendspherefiberContact = 0;
            int sendsphereplanewallContact = 0;
            int sendspherecylinderwallContact = 0;

            int sendfiberfiberContact = 0;
            int sendfiberplanewallContact = 0;
            int sendfibercylinderwallContact = 0;

            std::unordered_set<int> nodes;

            for (int gridId = offset[rankId * 2]; gridId <= offset[rankId * 2 + 1]; ++gridId)
            {

                for (auto i : sphereGrids[gridId])
                {
                    SphereParticleData sp;
                    sp.id = sphereparticles.at(i)->getId();
                    sp.category = sphereparticles.at(i)->getType().getCategory();
                    sp.subType = sphereparticles.at(i)->getType().getSubType();
                    sp.state = sphereparticles.at(i)->getState();

                    sp.position[0] = sphereparticles.at(i)->getPosition()[0];
                    sp.position[1] = sphereparticles.at(i)->getPosition()[1];
                    sp.position[2] = sphereparticles.at(i)->getPosition()[2];

                    sp.velocity[0] = sphereparticles.at(i)->getVelocity()[0];
                    sp.velocity[1] = sphereparticles.at(i)->getVelocity()[1];
                    sp.velocity[2] = sphereparticles.at(i)->getVelocity()[2];

                    sp.omega[0] = sphereparticles.at(i)->getOmega()[0];
                    sp.omega[1] = sphereparticles.at(i)->getOmega()[1];
                    sp.omega[2] = sphereparticles.at(i)->getOmega()[2];

                    sp.force[0] = sphereparticles.at(i)->getForce()[0];
                    sp.force[1] = sphereparticles.at(i)->getForce()[1];
                    sp.force[2] = sphereparticles.at(i)->getForce()[2];

                    sp.torque[0] = sphereparticles.at(i)->getTorque()[0];
                    sp.torque[1] = sphereparticles.at(i)->getTorque()[1];
                    sp.torque[2] = sphereparticles.at(i)->getTorque()[2];
                    sendsphereData.push_back(sp);
                    sendsphere++;
                    for (auto &contactInfor : spheresphereContactInformationList[i])
                    {
                        ContactInformationData ssc;
                        ssc.particleId1 = contactInfor.second.particleId1;
                        ssc.particleId2 = contactInfor.second.particleId2;

                        ssc.previousTangentialVelocity[0] = contactInfor.second.previousTangentialVelocity.x();
                        ssc.previousTangentialVelocity[1] = contactInfor.second.previousTangentialVelocity.y();
                        ssc.previousTangentialVelocity[2] = contactInfor.second.previousTangentialVelocity.z();

                        ssc.tangentialDisplacement[0] = contactInfor.second.tangentialDisplacement.x();
                        ssc.tangentialDisplacement[1] = contactInfor.second.tangentialDisplacement.y();
                        ssc.tangentialDisplacement[2] = contactInfor.second.tangentialDisplacement.z();

                        ssc.normalForce[0] = contactInfor.second.normalForce.x();
                        ssc.normalForce[1] = contactInfor.second.normalForce.y();
                        ssc.normalForce[2] = contactInfor.second.normalForce.z();

                        ssc.tangentialForce[0] = contactInfor.second.tangentialForce.x();
                        ssc.tangentialForce[1] = contactInfor.second.tangentialForce.y();
                        ssc.tangentialForce[2] = contactInfor.second.tangentialForce.z();
                        spherespherecontactinformationData.push_back(ssc);
                        sendspheresphereContact++;
                    }
                    for (auto &spherespherecontacts : spheresphereContactInformationList)
                    {
                        if (spherespherecontacts.second.find(i) != spherespherecontacts.second.end())
                        {
                            ContactInformationData ssc;
                            ssc.particleId1 = spherespherecontacts.second[i].particleId1;
                            ssc.particleId2 = spherespherecontacts.second[i].particleId2;

                            ssc.previousTangentialVelocity[0] = spherespherecontacts.second[i].previousTangentialVelocity.x();
                            ssc.previousTangentialVelocity[1] = spherespherecontacts.second[i].previousTangentialVelocity.y();
                            ssc.previousTangentialVelocity[2] = spherespherecontacts.second[i].previousTangentialVelocity.z();

                            ssc.tangentialDisplacement[0] = spherespherecontacts.second[i].tangentialDisplacement.x();
                            ssc.tangentialDisplacement[1] = spherespherecontacts.second[i].tangentialDisplacement.y();
                            ssc.tangentialDisplacement[2] = spherespherecontacts.second[i].tangentialDisplacement.z();

                            ssc.normalForce[0] = spherespherecontacts.second[i].normalForce.x();
                            ssc.normalForce[1] = spherespherecontacts.second[i].normalForce.y();
                            ssc.normalForce[2] = spherespherecontacts.second[i].normalForce.z();

                            ssc.tangentialForce[0] = spherespherecontacts.second[i].tangentialForce.x();
                            ssc.tangentialForce[1] = spherespherecontacts.second[i].tangentialForce.y();
                            ssc.tangentialForce[2] = spherespherecontacts.second[i].tangentialForce.z();
                            spherespherecontactinformationData.push_back(ssc);
                            sendspheresphereContact++;
                        }
                    }

                    for (auto &contactInfor : spherefiberContactInformationList[i])
                    {
                        ContactInformationData ssc;
                        ssc.particleId1 = contactInfor.second.particleId1;
                        ssc.particleId2 = contactInfor.second.particleId2;

                        ssc.previousTangentialVelocity[0] = contactInfor.second.previousTangentialVelocity.x();
                        ssc.previousTangentialVelocity[1] = contactInfor.second.previousTangentialVelocity.y();
                        ssc.previousTangentialVelocity[2] = contactInfor.second.previousTangentialVelocity.z();

                        ssc.tangentialDisplacement[0] = contactInfor.second.tangentialDisplacement.x();
                        ssc.tangentialDisplacement[1] = contactInfor.second.tangentialDisplacement.y();
                        ssc.tangentialDisplacement[2] = contactInfor.second.tangentialDisplacement.z();

                        ssc.normalForce[0] = contactInfor.second.normalForce.x();
                        ssc.normalForce[1] = contactInfor.second.normalForce.y();
                        ssc.normalForce[2] = contactInfor.second.normalForce.z();

                        ssc.tangentialForce[0] = contactInfor.second.tangentialForce.x();
                        ssc.tangentialForce[1] = contactInfor.second.tangentialForce.y();
                        ssc.tangentialForce[2] = contactInfor.second.tangentialForce.z();
                        spherefibercontactinformationData.push_back(ssc);
                        sendspherefiberContact++;
                    }
                    for (auto &sphereplanewllcontacts : sphereplanewallContactInformationList)
                    {
                        if (sphereplanewllcontacts.second.find(i) != sphereplanewllcontacts.second.end())
                        {
                            ContactInformationData ssc;
                            ssc.particleId1 = sphereplanewllcontacts.second[i].particleId1;
                            ssc.particleId2 = sphereplanewllcontacts.second[i].particleId2;

                            ssc.previousTangentialVelocity[0] = sphereplanewllcontacts.second[i].previousTangentialVelocity.x();
                            ssc.previousTangentialVelocity[1] = sphereplanewllcontacts.second[i].previousTangentialVelocity.y();
                            ssc.previousTangentialVelocity[2] = sphereplanewllcontacts.second[i].previousTangentialVelocity.z();

                            ssc.tangentialDisplacement[0] = sphereplanewllcontacts.second[i].tangentialDisplacement.x();
                            ssc.tangentialDisplacement[1] = sphereplanewllcontacts.second[i].tangentialDisplacement.y();
                            ssc.tangentialDisplacement[2] = sphereplanewllcontacts.second[i].tangentialDisplacement.z();

                            ssc.normalForce[0] = sphereplanewllcontacts.second[i].normalForce.x();
                            ssc.normalForce[1] = sphereplanewllcontacts.second[i].normalForce.y();
                            ssc.normalForce[2] = sphereplanewllcontacts.second[i].normalForce.z();

                            ssc.tangentialForce[0] = sphereplanewllcontacts.second[i].tangentialForce.x();
                            ssc.tangentialForce[1] = sphereplanewllcontacts.second[i].tangentialForce.y();
                            ssc.tangentialForce[2] = sphereplanewllcontacts.second[i].tangentialForce.z();
                            sphereplanewallcontactinformationData.push_back(ssc);
                            sendsphereplanewallContact++;
                        }
                    }
                    for (auto &spherecylinderwllcontacts : spherecylinderwallContactInformationList)
                    {
                        if (spherecylinderwllcontacts.second.find(i) != spherecylinderwllcontacts.second.end())
                        {
                            ContactInformationData ssc;
                            ssc.particleId1 = spherecylinderwllcontacts.second[i].particleId1;
                            ssc.particleId2 = spherecylinderwllcontacts.second[i].particleId2;

                            ssc.previousTangentialVelocity[0] = spherecylinderwllcontacts.second[i].previousTangentialVelocity.x();
                            ssc.previousTangentialVelocity[1] = spherecylinderwllcontacts.second[i].previousTangentialVelocity.y();
                            ssc.previousTangentialVelocity[2] = spherecylinderwllcontacts.second[i].previousTangentialVelocity.z();

                            ssc.tangentialDisplacement[0] = spherecylinderwllcontacts.second[i].tangentialDisplacement.x();
                            ssc.tangentialDisplacement[1] = spherecylinderwllcontacts.second[i].tangentialDisplacement.y();
                            ssc.tangentialDisplacement[2] = spherecylinderwllcontacts.second[i].tangentialDisplacement.z();

                            ssc.normalForce[0] = spherecylinderwllcontacts.second[i].normalForce.x();
                            ssc.normalForce[1] = spherecylinderwllcontacts.second[i].normalForce.y();
                            ssc.normalForce[2] = spherecylinderwllcontacts.second[i].normalForce.z();

                            ssc.tangentialForce[0] = spherecylinderwllcontacts.second[i].tangentialForce.x();
                            ssc.tangentialForce[1] = spherecylinderwllcontacts.second[i].tangentialForce.y();
                            ssc.tangentialForce[2] = spherecylinderwllcontacts.second[i].tangentialForce.z();
                            spherecylinderwallcontactinformationData.push_back(ssc);
                            sendspherecylinderwallContact++;
                        }
                    }
                }

                for (auto i : fiberGrids[gridId])
                {
                    FiberBondData fb;
                    nodes.insert(fiberbonds.at(i)->getNode1());
                    nodes.insert(fiberbonds.at(i)->getNode2());
                    fb.id = fiberbonds.at(i)->getId();
                    fb.category = fiberbonds.at(i)->getType().getCategory();
                    fb.subType = fiberbonds.at(i)->getType().getSubType();
                    fb.fiberid = fiberbonds.at(i)->getFiberId();
                    fb.node1 = fiberbonds.at(i)->getNode1();
                    fb.node2 = fiberbonds.at(i)->getNode2();
                    fb.neighborelement1 = fiberbonds.at(i)->getNeighborelement1();
                    fb.neighborelement2 = fiberbonds.at(i)->getNeighborelement2();
                    fb.energyDissipation = fiberbonds.at(i)->getEnergydissipation();
                    fb.tangentialforce[0] = fiberbonds.at(i)->getTangentialforce()[0];
                    fb.tangentialforce[1] = fiberbonds.at(i)->getTangentialforce()[1];
                    fb.tangentialforce[2] = fiberbonds.at(i)->getTangentialforce()[2];
                    fb.twisttorque[0] = fiberbonds.at(i)->getTwisttorque()[0];
                    fb.twisttorque[1] = fiberbonds.at(i)->getTwisttorque()[1];
                    fb.twisttorque[2] = fiberbonds.at(i)->getTwisttorque()[2];

                    fb.bendtorque[0] = fiberbonds.at(i)->getBendtorque()[0];
                    fb.bendtorque[1] = fiberbonds.at(i)->getBendtorque()[1];
                    fb.bendtorque[2] = fiberbonds.at(i)->getBendtorque()[2];
                    sendfiberbondData.push_back(fb);
                    sendfiberbond++;
                    for (auto &contactInfor : fiberfiberContactInformationList[i])
                    {
                        ContactInformationData ssc;
                        ssc.particleId1 = contactInfor.second.particleId1;
                        ssc.particleId2 = contactInfor.second.particleId2;

                        ssc.previousTangentialVelocity[0] = contactInfor.second.previousTangentialVelocity.x();
                        ssc.previousTangentialVelocity[1] = contactInfor.second.previousTangentialVelocity.y();
                        ssc.previousTangentialVelocity[2] = contactInfor.second.previousTangentialVelocity.z();

                        ssc.tangentialDisplacement[0] = contactInfor.second.tangentialDisplacement.x();
                        ssc.tangentialDisplacement[1] = contactInfor.second.tangentialDisplacement.y();
                        ssc.tangentialDisplacement[2] = contactInfor.second.tangentialDisplacement.z();

                        ssc.normalForce[0] = contactInfor.second.normalForce.x();
                        ssc.normalForce[1] = contactInfor.second.normalForce.y();
                        ssc.normalForce[2] = contactInfor.second.normalForce.z();

                        ssc.tangentialForce[0] = contactInfor.second.tangentialForce.x();
                        ssc.tangentialForce[1] = contactInfor.second.tangentialForce.y();
                        ssc.tangentialForce[2] = contactInfor.second.tangentialForce.z();
                        fiberfibercontactinformationData.push_back(ssc);
                        sendfiberfiberContact++;
                    }
                    for (auto &fiberfibercontacts : fiberfiberContactInformationList)
                    {
                        if (fiberfibercontacts.second.find(i) != fiberfibercontacts.second.end())
                        {
                            ContactInformationData ssc;
                            ssc.particleId1 = fiberfibercontacts.second[i].particleId1;
                            ssc.particleId2 = fiberfibercontacts.second[i].particleId2;

                            ssc.previousTangentialVelocity[0] = fiberfibercontacts.second[i].previousTangentialVelocity.x();
                            ssc.previousTangentialVelocity[1] = fiberfibercontacts.second[i].previousTangentialVelocity.y();
                            ssc.previousTangentialVelocity[2] = fiberfibercontacts.second[i].previousTangentialVelocity.z();

                            ssc.tangentialDisplacement[0] = fiberfibercontacts.second[i].tangentialDisplacement.x();
                            ssc.tangentialDisplacement[1] = fiberfibercontacts.second[i].tangentialDisplacement.y();
                            ssc.tangentialDisplacement[2] = fiberfibercontacts.second[i].tangentialDisplacement.z();

                            ssc.normalForce[0] = fiberfibercontacts.second[i].normalForce.x();
                            ssc.normalForce[1] = fiberfibercontacts.second[i].normalForce.y();
                            ssc.normalForce[2] = fiberfibercontacts.second[i].normalForce.z();

                            ssc.tangentialForce[0] = fiberfibercontacts.second[i].tangentialForce.x();
                            ssc.tangentialForce[1] = fiberfibercontacts.second[i].tangentialForce.y();
                            ssc.tangentialForce[2] = fiberfibercontacts.second[i].tangentialForce.z();
                            fiberfibercontactinformationData.push_back(ssc);
                            sendfiberfiberContact++;
                        }
                    }
                    for (auto &spherefibercontacts : spherefiberContactInformationList)
                    {
                        if (spherefibercontacts.second.find(i) != spherefibercontacts.second.end())
                        {
                            ContactInformationData ssc;
                            ssc.particleId1 = spherefibercontacts.second[i].particleId1;
                            ssc.particleId2 = spherefibercontacts.second[i].particleId2;

                            ssc.previousTangentialVelocity[0] = spherefibercontacts.second[i].previousTangentialVelocity.x();
                            ssc.previousTangentialVelocity[1] = spherefibercontacts.second[i].previousTangentialVelocity.y();
                            ssc.previousTangentialVelocity[2] = spherefibercontacts.second[i].previousTangentialVelocity.z();

                            ssc.tangentialDisplacement[0] = spherefibercontacts.second[i].tangentialDisplacement.x();
                            ssc.tangentialDisplacement[1] = spherefibercontacts.second[i].tangentialDisplacement.y();
                            ssc.tangentialDisplacement[2] = spherefibercontacts.second[i].tangentialDisplacement.z();

                            ssc.normalForce[0] = spherefibercontacts.second[i].normalForce.x();
                            ssc.normalForce[1] = spherefibercontacts.second[i].normalForce.y();
                            ssc.normalForce[2] = spherefibercontacts.second[i].normalForce.z();

                            ssc.tangentialForce[0] = spherefibercontacts.second[i].tangentialForce.x();
                            ssc.tangentialForce[1] = spherefibercontacts.second[i].tangentialForce.y();
                            ssc.tangentialForce[2] = spherefibercontacts.second[i].tangentialForce.z();
                            spherefibercontactinformationData.push_back(ssc);
                            sendspherefiberContact++;
                        }
                    }
                    for (auto &fiberplanewllcontacts : fiberplanewallContactInformationList)
                    {
                        if (fiberplanewllcontacts.second.find(i) != fiberplanewllcontacts.second.end())
                        {
                            ContactInformationData ssc;
                            ssc.particleId1 = fiberplanewllcontacts.second[i].particleId1;
                            ssc.particleId2 = fiberplanewllcontacts.second[i].particleId2;

                            ssc.previousTangentialVelocity[0] = fiberplanewllcontacts.second[i].previousTangentialVelocity.x();
                            ssc.previousTangentialVelocity[1] = fiberplanewllcontacts.second[i].previousTangentialVelocity.y();
                            ssc.previousTangentialVelocity[2] = fiberplanewllcontacts.second[i].previousTangentialVelocity.z();

                            ssc.tangentialDisplacement[0] = fiberplanewllcontacts.second[i].tangentialDisplacement.x();
                            ssc.tangentialDisplacement[1] = fiberplanewllcontacts.second[i].tangentialDisplacement.y();
                            ssc.tangentialDisplacement[2] = fiberplanewllcontacts.second[i].tangentialDisplacement.z();

                            ssc.normalForce[0] = fiberplanewllcontacts.second[i].normalForce.x();
                            ssc.normalForce[1] = fiberplanewllcontacts.second[i].normalForce.y();
                            ssc.normalForce[2] = fiberplanewllcontacts.second[i].normalForce.z();

                            ssc.tangentialForce[0] = fiberplanewllcontacts.second[i].tangentialForce.x();
                            ssc.tangentialForce[1] = fiberplanewllcontacts.second[i].tangentialForce.y();
                            ssc.tangentialForce[2] = fiberplanewllcontacts.second[i].tangentialForce.z();
                            fiberplanewallcontactinformationData.push_back(ssc);
                            sendfiberplanewallContact++;
                        }
                    }
                    for (auto &fibercylinderwllcontacts : fibercylinderwallContactInformationList)
                    {
                        if (fibercylinderwllcontacts.second.find(i) != fibercylinderwllcontacts.second.end())
                        {
                            ContactInformationData ssc;
                            ssc.particleId1 = fibercylinderwllcontacts.second[i].particleId1;
                            ssc.particleId2 = fibercylinderwllcontacts.second[i].particleId2;

                            ssc.previousTangentialVelocity[0] = fibercylinderwllcontacts.second[i].previousTangentialVelocity.x();
                            ssc.previousTangentialVelocity[1] = fibercylinderwllcontacts.second[i].previousTangentialVelocity.y();
                            ssc.previousTangentialVelocity[2] = fibercylinderwllcontacts.second[i].previousTangentialVelocity.z();

                            ssc.tangentialDisplacement[0] = fibercylinderwllcontacts.second[i].tangentialDisplacement.x();
                            ssc.tangentialDisplacement[1] = fibercylinderwllcontacts.second[i].tangentialDisplacement.y();
                            ssc.tangentialDisplacement[2] = fibercylinderwllcontacts.second[i].tangentialDisplacement.z();

                            ssc.normalForce[0] = fibercylinderwllcontacts.second[i].normalForce.x();
                            ssc.normalForce[1] = fibercylinderwllcontacts.second[i].normalForce.y();
                            ssc.normalForce[2] = fibercylinderwllcontacts.second[i].normalForce.z();

                            ssc.tangentialForce[0] = fibercylinderwllcontacts.second[i].tangentialForce.x();
                            ssc.tangentialForce[1] = fibercylinderwllcontacts.second[i].tangentialForce.y();
                            ssc.tangentialForce[2] = fibercylinderwllcontacts.second[i].tangentialForce.z();
                            fibercylinderwallcontactinformationData.push_back(ssc);
                            sendfibercylinderwallContact++;
                        }
                    }
                }
            }
            for (auto i : nodes)
            {
                SphereParticleData sp;
                sp.id = fibersphereparticles.at(i)->getId();
                sp.category = fibersphereparticles.at(i)->getType().getCategory();
                sp.subType = fibersphereparticles.at(i)->getType().getSubType();
                sp.state = fibersphereparticles.at(i)->getState();

                sp.position[0] = fibersphereparticles.at(i)->getPosition()[0];
                sp.position[1] = fibersphereparticles.at(i)->getPosition()[1];
                sp.position[2] = fibersphereparticles.at(i)->getPosition()[2];

                sp.velocity[0] = fibersphereparticles.at(i)->getVelocity()[0];
                sp.velocity[1] = fibersphereparticles.at(i)->getVelocity()[1];
                sp.velocity[2] = fibersphereparticles.at(i)->getVelocity()[2];

                sp.omega[0] = fibersphereparticles.at(i)->getOmega()[0];
                sp.omega[1] = fibersphereparticles.at(i)->getOmega()[1];
                sp.omega[2] = fibersphereparticles.at(i)->getOmega()[2];

                sp.force[0] = fibersphereparticles.at(i)->getForce()[0];
                sp.force[1] = fibersphereparticles.at(i)->getForce()[1];
                sp.force[2] = fibersphereparticles.at(i)->getForce()[2];

                sp.torque[0] = fibersphereparticles.at(i)->getTorque()[0];
                sp.torque[1] = fibersphereparticles.at(i)->getTorque()[1];
                sp.torque[2] = fibersphereparticles.at(i)->getTorque()[2];

                sendfibersphereData.push_back(sp);
                sendfibersphere++;
            }

            sendsphereCounts[rankId] = sendsphere;
            sendfiberbondCounts[rankId] = sendfiberbond;
            sendfibersphereCounts[rankId] = sendfibersphere;

            spherespherecontactCounts[rankId] = sendspheresphereContact;
            spherefibercontactCounts[rankId] = sendspherefiberContact;
            sphereplanewallcontactCounts[rankId] = sendsphereplanewallContact;
            spherecylinderwallcontactCounts[rankId] = sendspherecylinderwallContact;

            fiberfibercontactCounts[rankId] = sendspherefiberContact;
            fiberplanewallcontactCounts[rankId] = sendfiberplanewallContact;
            fibercylinderwallcontactCounts[rankId] = sendfibercylinderwallContact;

            if (rankId + 1 < world_size)
            {
                sendsphereDispls[rankId + 1] = static_cast<int>(sendsphereData.size());
                sendfibersphereDispls[rankId + 1] = static_cast<int>(sendfibersphereData.size());
                sendfiberbondDispls[rankId + 1] = static_cast<int>(sendfiberbondData.size());

                spherespherecontactDispls[rankId + 1] = static_cast<int>(spherespherecontactinformationData.size());
                spherefibercontactDispls[rankId + 1] = static_cast<int>(spherefibercontactinformationData.size());
                sphereplanewallcontactDispls[rankId + 1] = static_cast<int>(sphereplanewallcontactinformationData.size());
                spherecylinderwallcontactDispls[rankId + 1] = static_cast<int>(spherecylinderwallcontactinformationData.size());

                fiberfibercontactDispls[rankId + 1] = static_cast<int>(fiberfibercontactinformationData.size());
                fiberplanewallcontactDispls[rankId + 1] = static_cast<int>(fiberplanewallcontactinformationData.size());
                fibercylinderwallcontactDispls[rankId + 1] = static_cast<int>(fibercylinderwallcontactinformationData.size());
            }
        }
    }
    // all rank
    int recvsphereCount;
    int recvfibersphereCount;
    int recvfiberbondCount;

    int recvspherespherecontactCount;
    int recvspherefibercontactCount;
    int recvsphereplanewallcontactCount;
    int recvspherecylinderwallcontactCount;

    int recvfiberfibercontactCount;
    int recvfiberplanewallcontactCount;
    int recvfibercylinderwallcontactCount;

    MPI_Scatter(sendsphereCounts.data(), 1, MPI_INT, &recvsphereCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(sendfiberbondCounts.data(), 1, MPI_INT, &recvfiberbondCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(sendfibersphereCounts.data(), 1, MPI_INT, &recvfibersphereCount, 1, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Scatter(spherespherecontactCounts.data(), 1, MPI_INT, &recvspherespherecontactCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(spherefibercontactCounts.data(), 1, MPI_INT, &recvspherefibercontactCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(sphereplanewallcontactCounts.data(), 1, MPI_INT, &recvsphereplanewallcontactCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(spherecylinderwallcontactCounts.data(), 1, MPI_INT, &recvspherecylinderwallcontactCount, 1, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Scatter(fiberfibercontactCounts.data(), 1, MPI_INT, &recvfiberfibercontactCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(fiberplanewallcontactCounts.data(), 1, MPI_INT, &recvfiberplanewallcontactCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(fibercylinderwallcontactCounts.data(), 1, MPI_INT, &recvfibercylinderwallcontactCount, 1, MPI_INT, 0, MPI_COMM_WORLD);

    std::vector<SphereParticleData> receivesphereData(recvsphereCount);
    std::vector<SphereParticleData> receivefibersphereData(recvfibersphereCount);
    std::vector<FiberBondData> receivefiberbondData(recvfiberbondCount);

    std::vector<ContactInformationData> receivespherespherecontactData(recvspherespherecontactCount);
    std::vector<ContactInformationData> receivespherefibercontactData(recvspherefibercontactCount);
    std::vector<ContactInformationData> receivesphereplanewallcontactData(recvsphereplanewallcontactCount);
    std::vector<ContactInformationData> receivespherecylinderwallcontactData(recvspherecylinderwallcontactCount);

    std::vector<ContactInformationData> receivefiberfibercontactData(recvfiberfibercontactCount);
    std::vector<ContactInformationData> receivefiberplanewallcontactData(recvfiberplanewallcontactCount);
    std::vector<ContactInformationData> receivefibercylinderwallcontactData(recvfibercylinderwallcontactCount);

    MPI_Scatterv(sendsphereData.data(), sendsphereCounts.data(), sendsphereDispls.data(), mpi_sphereparticleDataType, receivesphereData.data(), recvsphereCount, mpi_sphereparticleDataType, 0, MPI_COMM_WORLD);
    MPI_Scatterv(sendfibersphereData.data(), sendfibersphereCounts.data(), sendfibersphereDispls.data(), mpi_sphereparticleDataType, receivefibersphereData.data(), recvfibersphereCount, mpi_sphereparticleDataType, 0, MPI_COMM_WORLD);
    MPI_Scatterv(sendfiberbondData.data(), sendfiberbondCounts.data(), sendfiberbondDispls.data(), mpi_fiberbondDataType, receivefiberbondData.data(), recvfiberbondCount, mpi_fiberbondDataType, 0, MPI_COMM_WORLD);

    MPI_Scatterv(spherespherecontactinformationData.data(), spherespherecontactCounts.data(), spherespherecontactDispls.data(),
                 mpi_contactinformationDataType, receivespherespherecontactData.data(), recvspherespherecontactCount,
                 mpi_contactinformationDataType, 0, MPI_COMM_WORLD);
    MPI_Scatterv(spherefibercontactinformationData.data(), spherefibercontactCounts.data(), spherefibercontactDispls.data(),
                 mpi_contactinformationDataType, receivespherefibercontactData.data(), recvspherefibercontactCount,
                 mpi_contactinformationDataType, 0, MPI_COMM_WORLD);
    MPI_Scatterv(sphereplanewallcontactinformationData.data(), sphereplanewallcontactCounts.data(), sphereplanewallcontactDispls.data(),
                 mpi_contactinformationDataType, receivesphereplanewallcontactData.data(), recvsphereplanewallcontactCount,
                 mpi_contactinformationDataType, 0, MPI_COMM_WORLD);
    MPI_Scatterv(spherecylinderwallcontactinformationData.data(), spherecylinderwallcontactCounts.data(), spherecylinderwallcontactDispls.data(),
                 mpi_contactinformationDataType, receivespherecylinderwallcontactData.data(), recvspherecylinderwallcontactCount,
                 mpi_contactinformationDataType, 0, MPI_COMM_WORLD);

    MPI_Scatterv(fiberfibercontactinformationData.data(), fiberfibercontactCounts.data(), fiberfibercontactDispls.data(),
                 mpi_contactinformationDataType, receivefiberfibercontactData.data(), recvfiberfibercontactCount,
                 mpi_contactinformationDataType, 0, MPI_COMM_WORLD);
    MPI_Scatterv(fiberplanewallcontactinformationData.data(), fiberplanewallcontactCounts.data(), fiberplanewallcontactDispls.data(),
                 mpi_contactinformationDataType, receivefiberplanewallcontactData.data(), recvfiberplanewallcontactCount,
                 mpi_contactinformationDataType, 0, MPI_COMM_WORLD);
    MPI_Scatterv(fibercylinderwallcontactinformationData.data(), fibercylinderwallcontactCounts.data(), fibercylinderwallcontactDispls.data(),
                 mpi_contactinformationDataType, receivefibercylinderwallcontactData.data(), recvfibercylinderwallcontactCount,
                 mpi_contactinformationDataType, 0, MPI_COMM_WORLD);
    std::unordered_map<int, std::unique_ptr<SphereParticle>> subsphereparticles;
    std::unordered_map<int, std::unique_ptr<SphereParticle>> subfibersphereparticles;
    std::unordered_map<int, std::unique_ptr<SphereCylinderBond>> subfiberbonds;

    for (const auto &data : receivesphereData)
    {
        Eigen::Vector3d position(data.position[0], data.position[1], data.position[2]);
        Eigen::Vector3d velocity(data.velocity[0], data.velocity[1], data.velocity[2]);
        Eigen::Vector3d omega(data.omega[0], data.omega[1], data.omega[2]);
        Eigen::Vector3d force(data.force[0], data.force[1], data.force[2]);
        Eigen::Vector3d torque(data.torque[0], data.torque[1], data.torque[2]);

        PropertyTypeID type = PropertyTypeID(data.category, data.subType);
        auto particle = std::make_unique<SphereParticle>(data.id, type, data.state, DEMproperties->getParticleManager(), position, velocity, omega, force, torque);
        subsphereparticles.insert({particle->getId(), std::move(particle)});
    }
    for (const auto &data : receivefibersphereData)
    {
        Eigen::Vector3d position(data.position[0], data.position[1], data.position[2]);
        Eigen::Vector3d velocity(data.velocity[0], data.velocity[1], data.velocity[2]);
        Eigen::Vector3d omega(data.omega[0], data.omega[1], data.omega[2]);
        Eigen::Vector3d force(data.force[0], data.force[1], data.force[2]);
        Eigen::Vector3d torque(data.torque[0], data.torque[1], data.torque[2]);

        PropertyTypeID type = PropertyTypeID(data.category, data.subType);
        auto particle = std::make_unique<SphereParticle>(
            data.id, type, data.state, DEMproperties->getParticleManager(), position, velocity, omega, force, torque);
        subfibersphereparticles.insert({particle->getId(), std::move(particle)});
    }
    for (const auto &data : receivefiberbondData)
    {
        Eigen::Vector3d tangentialforce(data.tangentialforce[0], data.tangentialforce[1], data.tangentialforce[2]);
        Eigen::Vector3d twisttorque(data.twisttorque[0], data.twisttorque[1], data.twisttorque[2]);
        Eigen::Vector3d bendtorque(data.bendtorque[0], data.bendtorque[1], data.bendtorque[2]);

        auto bond = std::make_unique<SphereCylinderBond>(
            data.id,
            PropertyTypeID(data.category, data.subType),
            DEMproperties->getParticleManager(),
            data.fiberid,
            data.node1,
            data.node2,
            data.neighborelement1,
            data.neighborelement2,
            data.energyDissipation);
        subfiberbonds.insert({bond->getId(), std::move(bond)});
    }
    DEMproperties->setsphereParticles(subsphereparticles);
    DEMproperties->setfibersphereParticles(subfibersphereparticles);
    DEMproperties->setfiberbonds(subfiberbonds);

    for (const auto &data : receivespherespherecontactData)
    {
        ContactInformation info;
        info.particleId1 = data.particleId1;
        info.particleId2 = data.particleId2;
        info.normalForce = Eigen::Vector3d(data.normalForce[0], data.normalForce[1], data.normalForce[2]);
        info.tangentialForce = Eigen::Vector3d(data.tangentialForce[0], data.tangentialForce[1], data.tangentialForce[2]);
        info.tangentialDisplacement = Eigen::Vector3d(data.tangentialDisplacement[0], data.tangentialDisplacement[1], data.tangentialDisplacement[2]);
        info.previousTangentialVelocity = Eigen::Vector3d(data.previousTangentialVelocity[0], data.previousTangentialVelocity[1], data.previousTangentialVelocity[2]);
        DEMproperties->getContactForce()->getSphereSphereContactInformationList()[info.particleId1][info.particleId2] = info;
    }
    for (const auto &data : receivespherefibercontactData)
    {
        ContactInformation info;
        info.particleId1 = data.particleId1;
        info.particleId2 = data.particleId2;
        info.normalForce = Eigen::Vector3d(data.normalForce[0], data.normalForce[1], data.normalForce[2]);
        info.tangentialForce = Eigen::Vector3d(data.tangentialForce[0], data.tangentialForce[1], data.tangentialForce[2]);
        info.tangentialDisplacement = Eigen::Vector3d(data.tangentialDisplacement[0], data.tangentialDisplacement[1], data.tangentialDisplacement[2]);
        info.previousTangentialVelocity = Eigen::Vector3d(data.previousTangentialVelocity[0], data.previousTangentialVelocity[1], data.previousTangentialVelocity[2]);
        DEMproperties->getContactForce()->getSphereFiberContactInformationList()[info.particleId1][info.particleId2] = info;
    }
    for (const auto &data : receivesphereplanewallcontactData)
    {
        ContactInformation info;
        info.particleId1 = data.particleId1;
        info.particleId2 = data.particleId2;
        info.normalForce = Eigen::Vector3d(data.normalForce[0], data.normalForce[1], data.normalForce[2]);
        info.tangentialForce = Eigen::Vector3d(data.tangentialForce[0], data.tangentialForce[1], data.tangentialForce[2]);
        info.tangentialDisplacement = Eigen::Vector3d(data.tangentialDisplacement[0], data.tangentialDisplacement[1], data.tangentialDisplacement[2]);
        info.previousTangentialVelocity = Eigen::Vector3d(data.previousTangentialVelocity[0], data.previousTangentialVelocity[1], data.previousTangentialVelocity[2]);
        DEMproperties->getContactForce()->getplanewallSphereContactInformationList()[info.particleId1][info.particleId2] = info;
    }
    for (const auto &data : receivespherecylinderwallcontactData)
    {
        ContactInformation info;
        info.particleId1 = data.particleId1;
        info.particleId2 = data.particleId2;
        info.normalForce = Eigen::Vector3d(data.normalForce[0], data.normalForce[1], data.normalForce[2]);
        info.tangentialForce = Eigen::Vector3d(data.tangentialForce[0], data.tangentialForce[1], data.tangentialForce[2]);
        info.tangentialDisplacement = Eigen::Vector3d(data.tangentialDisplacement[0], data.tangentialDisplacement[1], data.tangentialDisplacement[2]);
        info.previousTangentialVelocity = Eigen::Vector3d(data.previousTangentialVelocity[0], data.previousTangentialVelocity[1], data.previousTangentialVelocity[2]);
        DEMproperties->getContactForce()->getcylinderwallSphereContactInformationList()[info.particleId1][info.particleId2] = info;
    }

    for (const auto &data : receivefiberfibercontactData)
    {
        ContactInformation info;
        info.particleId1 = data.particleId1;
        info.particleId2 = data.particleId2;
        info.normalForce = Eigen::Vector3d(data.normalForce[0], data.normalForce[1], data.normalForce[2]);
        info.tangentialForce = Eigen::Vector3d(data.tangentialForce[0], data.tangentialForce[1], data.tangentialForce[2]);
        info.tangentialDisplacement = Eigen::Vector3d(data.tangentialDisplacement[0], data.tangentialDisplacement[1], data.tangentialDisplacement[2]);
        info.previousTangentialVelocity = Eigen::Vector3d(data.previousTangentialVelocity[0], data.previousTangentialVelocity[1], data.previousTangentialVelocity[2]);
        DEMproperties->getContactForce()->getFiberFiberContactInformationList()[info.particleId1][info.particleId2] = info;
    }
    for (const auto &data : receivefiberplanewallcontactData)
    {
        ContactInformation info;
        info.particleId1 = data.particleId1;
        info.particleId2 = data.particleId2;
        info.normalForce = Eigen::Vector3d(data.normalForce[0], data.normalForce[1], data.normalForce[2]);
        info.tangentialForce = Eigen::Vector3d(data.tangentialForce[0], data.tangentialForce[1], data.tangentialForce[2]);
        info.tangentialDisplacement = Eigen::Vector3d(data.tangentialDisplacement[0], data.tangentialDisplacement[1], data.tangentialDisplacement[2]);
        info.previousTangentialVelocity = Eigen::Vector3d(data.previousTangentialVelocity[0], data.previousTangentialVelocity[1], data.previousTangentialVelocity[2]);
        DEMproperties->getContactForce()->getplanewallFiberContactInformationList()[info.particleId1][info.particleId2] = info;
    }
    for (const auto &data : receivespherecylinderwallcontactData)
    {
        ContactInformation info;
        info.particleId1 = data.particleId1;
        info.particleId2 = data.particleId2;
        info.normalForce = Eigen::Vector3d(data.normalForce[0], data.normalForce[1], data.normalForce[2]);
        info.tangentialForce = Eigen::Vector3d(data.tangentialForce[0], data.tangentialForce[1], data.tangentialForce[2]);
        info.tangentialDisplacement = Eigen::Vector3d(data.tangentialDisplacement[0], data.tangentialDisplacement[1], data.tangentialDisplacement[2]);
        info.previousTangentialVelocity = Eigen::Vector3d(data.previousTangentialVelocity[0], data.previousTangentialVelocity[1], data.previousTangentialVelocity[2]);
        DEMproperties->getContactForce()->getcylinderwallFiberContactInformationList()[info.particleId1][info.particleId2] = info;
    }
}
void ParallelDEMModel::collectParticlesFromProcesses()
{
    std::vector<int> sendsphereCounts(world_size, 0);
    std::vector<int> sendfibersphereCounts(world_size, 0);
    std::vector<int> sendfiberbondCounts(world_size, 0);

    std::vector<int> sendsphereDispls(world_size, 0);
    std::vector<int> sendfibersphereDispls(world_size, 0);
    std::vector<int> sendfiberbondDispls(world_size, 0);

    std::vector<int> spherespherecontactCounts(world_size, 0);
    std::vector<int> spherefibercontactCounts(world_size, 0);
    std::vector<int> sphereplanewallcontactCounts(world_size, 0);
    std::vector<int> spherecylinderwallcontactCounts(world_size, 0);

    std::vector<int> spherespherecontactDispls(world_size, 0);
    std::vector<int> spherefibercontactDispls(world_size, 0);
    std::vector<int> sphereplanewallcontactDispls(world_size, 0);
    std::vector<int> spherecylinderwallcontactDispls(world_size, 0);

    std::vector<int> fiberfibercontactCounts(world_size, 0);
    std::vector<int> fiberplanewallcontactCounts(world_size, 0);
    std::vector<int> fibercylinderwallcontactCounts(world_size, 0);

    std::vector<int> fiberfibercontactDispls(world_size, 0);
    std::vector<int> fiberplanewallcontactDispls(world_size, 0);
    std::vector<int> fibercylinderwallcontactDispls(world_size, 0);

    std::vector<SphereParticleData> receivedSphereData;
    std::vector<SphereParticleData> receivedFiberSphereData;
    std::vector<FiberBondData> receivedFiberBondData;
    std::vector<ContactInformationData> receivedspherespherecontactinformationData;
    std::vector<ContactInformationData> receivedspherefibercontactinformationData;
    std::vector<ContactInformationData> receivedsphereplanewallcontactinformationData;
    std::vector<ContactInformationData> receivedspherecylinderwallcontactinformationData;
    std::vector<ContactInformationData> receivedfiberfibercontactinformationData;
    std::vector<ContactInformationData> receivedfiberplanewallcontactinformationData;
    std::vector<ContactInformationData> receivedfibercylinderwallcontactinformationData;

    std::unordered_set<int> nodes;
    std::vector<SphereParticleData> sendsphereData;
    std::vector<SphereParticleData> sendfibersphereData;
    std::vector<FiberBondData> sendfiberbondData;
    std::vector<ContactInformationData> sendspherespherecontactinformationData;
    std::vector<ContactInformationData> sendspherefibercontactinformationData;
    std::vector<ContactInformationData> sendsphereplanewallcontactinformationData;
    std::vector<ContactInformationData> sendspherecylinderwallcontactinformationData;
    std::vector<ContactInformationData> sendfiberfibercontactinformationData;
    std::vector<ContactInformationData> sendfiberplanewallcontactinformationData;
    std::vector<ContactInformationData> sendfibercylinderwallcontactinformationData;

    auto &sphereparticles = DEMproperties->getsphereParticles();
    auto &fiberbonds = DEMproperties->getfiberbonds();
    auto &fibersphereparticles = DEMproperties->getfibersphereParticles();

    auto &spheresphereContactInformationList = DEMproperties->getContactForce()->getSphereSphereContactInformationList();
    auto &spherefiberContactInformationList = DEMproperties->getContactForce()->getSphereFiberContactInformationList();
    auto &sphereplanewallContactInformationList = DEMproperties->getContactForce()->getplanewallSphereContactInformationList();
    auto &spherecylinderwallContactInformationList = DEMproperties->getContactForce()->getcylinderwallSphereContactInformationList();

    auto &fiberfiberContactInformationList = DEMproperties->getContactForce()->getFiberFiberContactInformationList();
    auto &fiberplanewallContactInformationList = DEMproperties->getContactForce()->getplanewallFiberContactInformationList();
    auto &fibercylinderwallContactInformationList = DEMproperties->getContactForce()->getcylinderwallFiberContactInformationList();

    for (auto &particle : sphereparticles)
    {
        SphereParticleData sp;
        sp.id = particle.second->getId();
        sp.category = particle.second->getType().getCategory();
        sp.subType = particle.second->getType().getSubType();
        sp.state = particle.second->getState();

        sp.position[0] = particle.second->getPosition()[0];
        sp.position[1] = particle.second->getPosition()[1];
        sp.position[2] = particle.second->getPosition()[2];

        sp.velocity[0] = particle.second->getVelocity()[0];
        sp.velocity[1] = particle.second->getVelocity()[1];
        sp.velocity[2] = particle.second->getVelocity()[2];

        sp.omega[0] = particle.second->getOmega()[0];
        sp.omega[1] = particle.second->getOmega()[1];
        sp.omega[2] = particle.second->getOmega()[2];

        sp.force[0] = particle.second->getForce()[0];
        sp.force[1] = particle.second->getForce()[1];
        sp.force[2] = particle.second->getForce()[2];

        sp.torque[0] = particle.second->getTorque()[0];
        sp.torque[1] = particle.second->getTorque()[1];
        sp.torque[2] = particle.second->getTorque()[2];
        sendsphereData.push_back(sp);
    }
    for (auto &particle : fiberbonds)
    {
        FiberBondData fb;
        nodes.insert(particle.second->getNode1());
        nodes.insert(particle.second->getNode2());
        fb.id = particle.second->getId();
        fb.category = particle.second->getType().getCategory();
        fb.subType = particle.second->getType().getSubType();
        fb.fiberid = particle.second->getFiberId();
        fb.node1 = particle.second->getNode1();
        fb.node2 = particle.second->getNode2();
        fb.neighborelement1 = particle.second->getNeighborelement1();
        fb.neighborelement2 = particle.second->getNeighborelement2();
        fb.energyDissipation = particle.second->getEnergydissipation();
        fb.tangentialforce[0] = particle.second->getTangentialforce()[0];
        fb.tangentialforce[1] = particle.second->getTangentialforce()[1];
        fb.tangentialforce[2] = particle.second->getTangentialforce()[2];
        fb.twisttorque[0] = particle.second->getTwisttorque()[0];
        fb.twisttorque[1] = particle.second->getTwisttorque()[1];
        fb.twisttorque[2] = particle.second->getTwisttorque()[2];

        fb.bendtorque[0] = particle.second->getBendtorque()[0];
        fb.bendtorque[1] = particle.second->getBendtorque()[1];
        fb.bendtorque[2] = particle.second->getBendtorque()[2];
        sendfiberbondData.push_back(fb);
    }
    for (auto i : nodes)
    {
        SphereParticleData sp;
        sp.id = fibersphereparticles.at(i)->getId();
        sp.category = fibersphereparticles.at(i)->getType().getCategory();
        sp.subType = fibersphereparticles.at(i)->getType().getSubType();
        sp.state = fibersphereparticles.at(i)->getState();

        sp.position[0] = fibersphereparticles.at(i)->getPosition()[0];
        sp.position[1] = fibersphereparticles.at(i)->getPosition()[1];
        sp.position[2] = fibersphereparticles.at(i)->getPosition()[2];

        sp.velocity[0] = fibersphereparticles.at(i)->getVelocity()[0];
        sp.velocity[1] = fibersphereparticles.at(i)->getVelocity()[1];
        sp.velocity[2] = fibersphereparticles.at(i)->getVelocity()[2];

        sp.omega[0] = fibersphereparticles.at(i)->getOmega()[0];
        sp.omega[1] = fibersphereparticles.at(i)->getOmega()[1];
        sp.omega[2] = fibersphereparticles.at(i)->getOmega()[2];

        sp.force[0] = fibersphereparticles.at(i)->getForce()[0];
        sp.force[1] = fibersphereparticles.at(i)->getForce()[1];
        sp.force[2] = fibersphereparticles.at(i)->getForce()[2];

        sp.torque[0] = fibersphereparticles.at(i)->getTorque()[0];
        sp.torque[1] = fibersphereparticles.at(i)->getTorque()[1];
        sp.torque[2] = fibersphereparticles.at(i)->getTorque()[2];

        sendfibersphereData.push_back(sp);
    }

    for (auto &spherespherecontacts : spheresphereContactInformationList)
    {
        std::unordered_map<int, ContactInformation> &innerMap = spherespherecontacts.second;
        for (auto &innerPair : innerMap)
        {
            ContactInformationData ssc;
            ssc.particleId1 = innerPair.second.particleId1;
            ssc.particleId2 = innerPair.second.particleId2;

            ssc.previousTangentialVelocity[0] = innerPair.second.previousTangentialVelocity.x();
            ssc.previousTangentialVelocity[1] = innerPair.second.previousTangentialVelocity.y();
            ssc.previousTangentialVelocity[2] = innerPair.second.previousTangentialVelocity.z();

            ssc.tangentialDisplacement[0] = innerPair.second.tangentialDisplacement.x();
            ssc.tangentialDisplacement[1] = innerPair.second.tangentialDisplacement.y();
            ssc.tangentialDisplacement[2] = innerPair.second.tangentialDisplacement.z();

            ssc.normalForce[0] = innerPair.second.normalForce.x();
            ssc.normalForce[1] = innerPair.second.normalForce.y();
            ssc.normalForce[2] = innerPair.second.normalForce.z();

            ssc.tangentialForce[0] = innerPair.second.tangentialForce.x();
            ssc.tangentialForce[1] = innerPair.second.tangentialForce.y();
            ssc.tangentialForce[2] = innerPair.second.tangentialForce.z();
            sendspherespherecontactinformationData.push_back(ssc);
        }
    }
    for (auto &spherefibercontacts : spherefiberContactInformationList)
    {
        std::unordered_map<int, ContactInformation> &innerMap = spherefibercontacts.second;
        for (auto &innerPair : innerMap)
        {
            ContactInformationData ssc;
            ssc.particleId1 = innerPair.second.particleId1;
            ssc.particleId2 = innerPair.second.particleId2;

            ssc.previousTangentialVelocity[0] = innerPair.second.previousTangentialVelocity.x();
            ssc.previousTangentialVelocity[1] = innerPair.second.previousTangentialVelocity.y();
            ssc.previousTangentialVelocity[2] = innerPair.second.previousTangentialVelocity.z();

            ssc.tangentialDisplacement[0] = innerPair.second.tangentialDisplacement.x();
            ssc.tangentialDisplacement[1] = innerPair.second.tangentialDisplacement.y();
            ssc.tangentialDisplacement[2] = innerPair.second.tangentialDisplacement.z();

            ssc.normalForce[0] = innerPair.second.normalForce.x();
            ssc.normalForce[1] = innerPair.second.normalForce.y();
            ssc.normalForce[2] = innerPair.second.normalForce.z();

            ssc.tangentialForce[0] = innerPair.second.tangentialForce.x();
            ssc.tangentialForce[1] = innerPair.second.tangentialForce.y();
            ssc.tangentialForce[2] = innerPair.second.tangentialForce.z();
            sendspherefibercontactinformationData.push_back(ssc);
        }
    }
    for (auto &sphereplanewllcontacts : sphereplanewallContactInformationList)
    {
        std::unordered_map<int, ContactInformation> &innerMap = sphereplanewllcontacts.second;
        for (auto &innerPair : innerMap)
        {
            ContactInformationData ssc;
            ssc.particleId1 = innerPair.second.particleId1;
            ssc.particleId2 = innerPair.second.particleId2;

            ssc.previousTangentialVelocity[0] = innerPair.second.previousTangentialVelocity.x();
            ssc.previousTangentialVelocity[1] = innerPair.second.previousTangentialVelocity.y();
            ssc.previousTangentialVelocity[2] = innerPair.second.previousTangentialVelocity.z();

            ssc.tangentialDisplacement[0] = innerPair.second.tangentialDisplacement.x();
            ssc.tangentialDisplacement[1] = innerPair.second.tangentialDisplacement.y();
            ssc.tangentialDisplacement[2] = innerPair.second.tangentialDisplacement.z();

            ssc.normalForce[0] = innerPair.second.normalForce.x();
            ssc.normalForce[1] = innerPair.second.normalForce.y();
            ssc.normalForce[2] = innerPair.second.normalForce.z();

            ssc.tangentialForce[0] = innerPair.second.tangentialForce.x();
            ssc.tangentialForce[1] = innerPair.second.tangentialForce.y();
            ssc.tangentialForce[2] = innerPair.second.tangentialForce.z();
            sendsphereplanewallcontactinformationData.push_back(ssc);
        }
    }
    for (auto &spherecylinderwllcontacts : spherecylinderwallContactInformationList)
    {
        std::unordered_map<int, ContactInformation> &innerMap = spherecylinderwllcontacts.second;
        for (auto &innerPair : innerMap)
        {
            ContactInformationData ssc;
            ssc.particleId1 = innerPair.second.particleId1;
            ssc.particleId2 = innerPair.second.particleId2;

            ssc.previousTangentialVelocity[0] = innerPair.second.previousTangentialVelocity.x();
            ssc.previousTangentialVelocity[1] = innerPair.second.previousTangentialVelocity.y();
            ssc.previousTangentialVelocity[2] = innerPair.second.previousTangentialVelocity.z();

            ssc.tangentialDisplacement[0] = innerPair.second.tangentialDisplacement.x();
            ssc.tangentialDisplacement[1] = innerPair.second.tangentialDisplacement.y();
            ssc.tangentialDisplacement[2] = innerPair.second.tangentialDisplacement.z();

            ssc.normalForce[0] = innerPair.second.normalForce.x();
            ssc.normalForce[1] = innerPair.second.normalForce.y();
            ssc.normalForce[2] = innerPair.second.normalForce.z();

            ssc.tangentialForce[0] = innerPair.second.tangentialForce.x();
            ssc.tangentialForce[1] = innerPair.second.tangentialForce.y();
            ssc.tangentialForce[2] = innerPair.second.tangentialForce.z();

            sendspherecylinderwallcontactinformationData.push_back(ssc);
        }
    }
    for (auto &fiberfibercontacts : fiberfiberContactInformationList)
    {
        std::unordered_map<int, ContactInformation> &innerMap = fiberfibercontacts.second;
        for (auto &innerPair : innerMap)
        {
            ContactInformationData ssc;
            ssc.particleId1 = innerPair.second.particleId1;
            ssc.particleId2 = innerPair.second.particleId2;

            ssc.previousTangentialVelocity[0] = innerPair.second.previousTangentialVelocity.x();
            ssc.previousTangentialVelocity[1] = innerPair.second.previousTangentialVelocity.y();
            ssc.previousTangentialVelocity[2] = innerPair.second.previousTangentialVelocity.z();

            ssc.tangentialDisplacement[0] = innerPair.second.tangentialDisplacement.x();
            ssc.tangentialDisplacement[1] = innerPair.second.tangentialDisplacement.y();
            ssc.tangentialDisplacement[2] = innerPair.second.tangentialDisplacement.z();

            ssc.normalForce[0] = innerPair.second.normalForce.x();
            ssc.normalForce[1] = innerPair.second.normalForce.y();
            ssc.normalForce[2] = innerPair.second.normalForce.z();

            ssc.tangentialForce[0] = innerPair.second.tangentialForce.x();
            ssc.tangentialForce[1] = innerPair.second.tangentialForce.y();
            ssc.tangentialForce[2] = innerPair.second.tangentialForce.z();
            sendfiberfibercontactinformationData.push_back(ssc);
        }
    }
    for (auto &fiberplanewllcontacts : fiberplanewallContactInformationList)
    {
        std::unordered_map<int, ContactInformation> &innerMap = fiberplanewllcontacts.second;
        for (auto &innerPair : innerMap)
        {
            ContactInformationData ssc;
            ssc.particleId1 = innerPair.second.particleId1;
            ssc.particleId2 = innerPair.second.particleId2;

            ssc.previousTangentialVelocity[0] = innerPair.second.previousTangentialVelocity.x();
            ssc.previousTangentialVelocity[1] = innerPair.second.previousTangentialVelocity.y();
            ssc.previousTangentialVelocity[2] = innerPair.second.previousTangentialVelocity.z();

            ssc.tangentialDisplacement[0] = innerPair.second.tangentialDisplacement.x();
            ssc.tangentialDisplacement[1] = innerPair.second.tangentialDisplacement.y();
            ssc.tangentialDisplacement[2] = innerPair.second.tangentialDisplacement.z();

            ssc.normalForce[0] = innerPair.second.normalForce.x();
            ssc.normalForce[1] = innerPair.second.normalForce.y();
            ssc.normalForce[2] = innerPair.second.normalForce.z();

            ssc.tangentialForce[0] = innerPair.second.tangentialForce.x();
            ssc.tangentialForce[1] = innerPair.second.tangentialForce.y();
            ssc.tangentialForce[2] = innerPair.second.tangentialForce.z();
            sendfiberplanewallcontactinformationData.push_back(ssc);
        }
    }
    for (auto &fibercylinderwllcontacts : fibercylinderwallContactInformationList)
    {
        std::unordered_map<int, ContactInformation> &innerMap = fibercylinderwllcontacts.second;
        for (auto &innerPair : innerMap)
        {
            ContactInformationData ssc;
            ssc.particleId1 = innerPair.second.particleId1;
            ssc.particleId2 = innerPair.second.particleId2;

            ssc.previousTangentialVelocity[0] = innerPair.second.previousTangentialVelocity.x();
            ssc.previousTangentialVelocity[1] = innerPair.second.previousTangentialVelocity.y();
            ssc.previousTangentialVelocity[2] = innerPair.second.previousTangentialVelocity.z();

            ssc.tangentialDisplacement[0] = innerPair.second.tangentialDisplacement.x();
            ssc.tangentialDisplacement[1] = innerPair.second.tangentialDisplacement.y();
            ssc.tangentialDisplacement[2] = innerPair.second.tangentialDisplacement.z();

            ssc.normalForce[0] = innerPair.second.normalForce.x();
            ssc.normalForce[1] = innerPair.second.normalForce.y();
            ssc.normalForce[2] = innerPair.second.normalForce.z();

            ssc.tangentialForce[0] = innerPair.second.tangentialForce.x();
            ssc.tangentialForce[1] = innerPair.second.tangentialForce.y();
            ssc.tangentialForce[2] = innerPair.second.tangentialForce.z();

            sendfibercylinderwallcontactinformationData.push_back(ssc);
        }
    }

    int sendSphereCount = static_cast<int>(sendsphereData.size());
    int sendFiberSphereCount = static_cast<int>(sendfibersphereData.size());
    int sendFiberBondCount = static_cast<int>(sendfiberbondData.size());

    int sendspherespherecontactCount = static_cast<int>(sendspherespherecontactinformationData.size());
    int sendspherefibercontactCount = static_cast<int>(sendspherefibercontactinformationData.size());
    int sendsphereplanewallcontactCount = static_cast<int>(sendsphereplanewallcontactinformationData.size());
    int sendspherecylinderwallcontactCount = static_cast<int>(sendspherecylinderwallcontactinformationData.size());
    int sendfiberfibercontactCount = static_cast<int>(sendfiberfibercontactinformationData.size());
    int sendfiberplanewallcontactCount = static_cast<int>(sendfiberplanewallcontactinformationData.size());
    int sendfibercylinderwallcontactCount = static_cast<int>(sendfibercylinderwallcontactinformationData.size());

    MPI_Gather(&sendSphereCount, 1, MPI_INT, sendsphereCounts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&sendFiberSphereCount, 1, MPI_INT, sendfibersphereCounts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&sendFiberBondCount, 1, MPI_INT, sendfiberbondCounts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Gather(&sendspherespherecontactCount, 1, MPI_INT, spherespherecontactCounts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&sendspherefibercontactCount, 1, MPI_INT, spherefibercontactCounts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&sendsphereplanewallcontactCount, 1, MPI_INT, sphereplanewallcontactCounts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&sendspherecylinderwallcontactCount, 1, MPI_INT, spherecylinderwallcontactCounts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Gather(&sendfiberfibercontactCount, 1, MPI_INT, fiberfibercontactCounts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&sendfiberplanewallcontactCount, 1, MPI_INT, fiberplanewallcontactCounts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&sendfibercylinderwallcontactCount, 1, MPI_INT, fibercylinderwallcontactCounts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (world_rank == 0)
    {
        int totalSphereCount = sendsphereCounts[0];
        int totalFiberSphereCount = sendfibersphereCounts[0];
        int totalFiberBondCount = sendfiberbondCounts[0];
        int totalspherespherecontactCount = spherespherecontactCounts[0];
        int totalspherefibercontactCount = spherefibercontactCounts[0];
        int totalsphereplanewallcontactCount = sphereplanewallcontactCounts[0];
        int totalspherecylinderwallcontactCount = spherecylinderwallcontactCounts[0];
        int totalfiberfibercontactCount = fiberfibercontactCounts[0];
        int totalfiberplanewallcontactCount = fiberplanewallcontactCounts[0];
        int totalfibercylinderwallcontactCount = fibercylinderwallcontactCounts[0];
        for (int i = 1; i < sendsphereCounts.size(); ++i)
        {
            sendsphereDispls[i] = sendsphereDispls[i - 1] + sendsphereCounts[i - 1];
            sendfibersphereDispls[i] = sendfibersphereDispls[i - 1] + sendfibersphereCounts[i - 1];
            sendfiberbondDispls[i] = sendfiberbondDispls[i - 1] + sendfiberbondCounts[i - 1];
            spherespherecontactDispls[i] = spherespherecontactDispls[i - 1] + spherespherecontactCounts[i - 1];
            spherefibercontactDispls[i] = spherefibercontactDispls[i - 1] + spherefibercontactCounts[i - 1];
            sphereplanewallcontactDispls[i] = sphereplanewallcontactDispls[i - 1] + sphereplanewallcontactCounts[i - 1];
            spherecylinderwallcontactDispls[i] = spherecylinderwallcontactDispls[i - 1] + spherecylinderwallcontactCounts[i - 1];

            fiberfibercontactDispls[i] = fiberfibercontactDispls[i - 1] + fiberfibercontactCounts[i - 1];
            fiberplanewallcontactDispls[i] = fiberplanewallcontactDispls[i - 1] + fiberplanewallcontactCounts[i - 1];
            fibercylinderwallcontactDispls[i] = fibercylinderwallcontactDispls[i - 1] + fibercylinderwallcontactCounts[i - 1];

            totalSphereCount += sendsphereCounts[i];
            totalFiberSphereCount += sendfibersphereCounts[i];
            totalFiberBondCount += sendfiberbondCounts[i];
            totalspherespherecontactCount += spherespherecontactCounts[i];
            totalspherefibercontactCount += spherefibercontactCounts[i];
            totalsphereplanewallcontactCount += sphereplanewallcontactCounts[i];
            totalspherecylinderwallcontactCount += spherecylinderwallcontactCounts[i];
            totalfiberfibercontactCount += fiberfibercontactCounts[i];
            totalfiberplanewallcontactCount += fiberplanewallcontactCounts[i];
            totalfibercylinderwallcontactCount += fibercylinderwallcontactCounts[i];
        }
        receivedSphereData.resize(totalSphereCount);
        receivedFiberSphereData.resize(totalFiberSphereCount);
        receivedFiberBondData.resize(totalFiberBondCount);
        receivedspherespherecontactinformationData.resize(totalspherespherecontactCount);
        receivedspherefibercontactinformationData.resize(totalspherefibercontactCount);
        receivedsphereplanewallcontactinformationData.resize(totalsphereplanewallcontactCount);
        receivedspherecylinderwallcontactinformationData.resize(totalspherecylinderwallcontactCount);

        receivedfiberfibercontactinformationData.resize(totalfiberfibercontactCount);
        receivedfiberplanewallcontactinformationData.resize(totalfiberplanewallcontactCount);
        receivedfibercylinderwallcontactinformationData.resize(totalfibercylinderwallcontactCount);
    }

    MPI_Gatherv(sendsphereData.data(), sendSphereCount, mpi_sphereparticleDataType, receivedSphereData.data(), sendsphereCounts.data(), sendsphereDispls.data(), mpi_sphereparticleDataType, 0, MPI_COMM_WORLD);
    MPI_Gatherv(sendfibersphereData.data(), sendFiberSphereCount, mpi_sphereparticleDataType, receivedFiberSphereData.data(), sendfibersphereCounts.data(), sendfibersphereDispls.data(), mpi_sphereparticleDataType, 0, MPI_COMM_WORLD);
    MPI_Gatherv(sendfiberbondData.data(), sendFiberBondCount, mpi_fiberbondDataType, receivedFiberBondData.data(), sendfiberbondCounts.data(), sendfiberbondDispls.data(), mpi_fiberbondDataType, 0, MPI_COMM_WORLD);

    MPI_Gatherv(sendspherespherecontactinformationData.data(), sendspherespherecontactCount, mpi_contactinformationDataType,
                receivedspherespherecontactinformationData.data(), spherespherecontactCounts.data(),
                spherespherecontactDispls.data(), mpi_contactinformationDataType, 0, MPI_COMM_WORLD);
    MPI_Gatherv(sendspherefibercontactinformationData.data(), sendspherefibercontactCount, mpi_contactinformationDataType,
                receivedspherefibercontactinformationData.data(), spherefibercontactCounts.data(),
                spherefibercontactDispls.data(), mpi_contactinformationDataType, 0, MPI_COMM_WORLD);
    MPI_Gatherv(sendsphereplanewallcontactinformationData.data(), sendsphereplanewallcontactCount, mpi_contactinformationDataType,
                receivedsphereplanewallcontactinformationData.data(), sphereplanewallcontactCounts.data(),
                sphereplanewallcontactDispls.data(), mpi_contactinformationDataType, 0, MPI_COMM_WORLD);
    MPI_Gatherv(sendspherecylinderwallcontactinformationData.data(), sendspherecylinderwallcontactCount, mpi_contactinformationDataType,
                receivedspherecylinderwallcontactinformationData.data(), spherecylinderwallcontactCounts.data(),
                spherecylinderwallcontactDispls.data(), mpi_contactinformationDataType, 0, MPI_COMM_WORLD);
    MPI_Gatherv(sendfiberfibercontactinformationData.data(), sendfiberfibercontactCount, mpi_contactinformationDataType,
                receivedfiberfibercontactinformationData.data(), fiberfibercontactCounts.data(),
                fiberfibercontactDispls.data(), mpi_contactinformationDataType, 0, MPI_COMM_WORLD);
    MPI_Gatherv(sendfiberplanewallcontactinformationData.data(), sendfiberplanewallcontactCount, mpi_contactinformationDataType,
                receivedfiberplanewallcontactinformationData.data(), fiberplanewallcontactCounts.data(),
                fiberplanewallcontactDispls.data(), mpi_contactinformationDataType, 0, MPI_COMM_WORLD);
    MPI_Gatherv(sendfibercylinderwallcontactinformationData.data(), sendfibercylinderwallcontactCount, mpi_contactinformationDataType,
                receivedfibercylinderwallcontactinformationData.data(), fibercylinderwallcontactCounts.data(),
                fibercylinderwallcontactDispls.data(), mpi_contactinformationDataType, 0, MPI_COMM_WORLD);

    if (world_rank == 0)
    {
        std::unordered_map<int, std::unique_ptr<SphereParticle>> newsphereparticles;
        std::unordered_map<int, std::unique_ptr<SphereParticle>> newfibersphereparticles;
        std::unordered_map<int, std::unique_ptr<SphereCylinderBond>> newfiberbonds;
        for (const auto &data : receivedSphereData)
        {
            Eigen::Vector3d position(data.position[0], data.position[1], data.position[2]);
            Eigen::Vector3d velocity(data.velocity[0], data.velocity[1], data.velocity[2]);
            Eigen::Vector3d omega(data.omega[0], data.omega[1], data.omega[2]);
            Eigen::Vector3d force(data.force[0], data.force[1], data.force[2]);
            Eigen::Vector3d torque(data.torque[0], data.torque[1], data.torque[2]);

            PropertyTypeID type = PropertyTypeID(data.category, data.subType);
            auto particle = std::make_unique<SphereParticle>(data.id, type, data.state, DEMproperties->getParticleManager(), position, velocity, omega, force, torque);
            newsphereparticles.insert({particle->getId(), std::move(particle)});
        }
        for (const auto &data : receivedFiberSphereData)
        {
            Eigen::Vector3d position(data.position[0], data.position[1], data.position[2]);
            Eigen::Vector3d velocity(data.velocity[0], data.velocity[1], data.velocity[2]);
            Eigen::Vector3d omega(data.omega[0], data.omega[1], data.omega[2]);
            Eigen::Vector3d force(data.force[0], data.force[1], data.force[2]);
            Eigen::Vector3d torque(data.torque[0], data.torque[1], data.torque[2]);

            PropertyTypeID type = PropertyTypeID(data.category, data.subType);
            auto particle = std::make_unique<SphereParticle>(
                data.id, type, data.state, DEMproperties->getParticleManager(), position, velocity, omega, force, torque);
            newfibersphereparticles.insert({particle->getId(), std::move(particle)});
        }
        for (const auto &data : receivedFiberBondData)
        {
            Eigen::Vector3d tangentialforce(data.tangentialforce[0], data.tangentialforce[1], data.tangentialforce[2]);
            Eigen::Vector3d twisttorque(data.twisttorque[0], data.twisttorque[1], data.twisttorque[2]);
            Eigen::Vector3d bendtorque(data.bendtorque[0], data.bendtorque[1], data.bendtorque[2]);

            auto bond = std::make_unique<SphereCylinderBond>(
                data.id,
                PropertyTypeID(data.category, data.subType),
                DEMproperties->getParticleManager(),
                data.fiberid,
                data.node1,
                data.node2,
                data.neighborelement1,
                data.neighborelement2,
                data.energyDissipation);
            newfiberbonds.insert({bond->getId(), std::move(bond)});
        }
        DEMproperties->setsphereParticles(newsphereparticles);
        DEMproperties->setfibersphereParticles(newfibersphereparticles);
        DEMproperties->setfiberbonds(newfiberbonds);

        for (const auto &data : receivedspherespherecontactinformationData)
        {
            ContactInformation info;
            info.particleId1 = data.particleId1;
            info.particleId2 = data.particleId2;
            info.normalForce = Eigen::Vector3d(data.normalForce[0], data.normalForce[1], data.normalForce[2]);
            info.tangentialForce = Eigen::Vector3d(data.tangentialForce[0], data.tangentialForce[1], data.tangentialForce[2]);
            info.tangentialDisplacement = Eigen::Vector3d(data.tangentialDisplacement[0], data.tangentialDisplacement[1], data.tangentialDisplacement[2]);
            info.previousTangentialVelocity = Eigen::Vector3d(data.previousTangentialVelocity[0], data.previousTangentialVelocity[1], data.previousTangentialVelocity[2]);
            DEMproperties->getContactForce()->getSphereSphereContactInformationList()[info.particleId1][info.particleId2] = info;
        }
        for (const auto &data : receivedspherefibercontactinformationData)
        {
            ContactInformation info;
            info.particleId1 = data.particleId1;
            info.particleId2 = data.particleId2;
            info.normalForce = Eigen::Vector3d(data.normalForce[0], data.normalForce[1], data.normalForce[2]);
            info.tangentialForce = Eigen::Vector3d(data.tangentialForce[0], data.tangentialForce[1], data.tangentialForce[2]);
            info.tangentialDisplacement = Eigen::Vector3d(data.tangentialDisplacement[0], data.tangentialDisplacement[1], data.tangentialDisplacement[2]);
            info.previousTangentialVelocity = Eigen::Vector3d(data.previousTangentialVelocity[0], data.previousTangentialVelocity[1], data.previousTangentialVelocity[2]);
            DEMproperties->getContactForce()->getSphereFiberContactInformationList()[info.particleId1][info.particleId2] = info;
        }
        for (const auto &data : receivedsphereplanewallcontactinformationData)
        {
            ContactInformation info;
            info.particleId1 = data.particleId1;
            info.particleId2 = data.particleId2;
            info.normalForce = Eigen::Vector3d(data.normalForce[0], data.normalForce[1], data.normalForce[2]);
            info.tangentialForce = Eigen::Vector3d(data.tangentialForce[0], data.tangentialForce[1], data.tangentialForce[2]);
            info.tangentialDisplacement = Eigen::Vector3d(data.tangentialDisplacement[0], data.tangentialDisplacement[1], data.tangentialDisplacement[2]);
            info.previousTangentialVelocity = Eigen::Vector3d(data.previousTangentialVelocity[0], data.previousTangentialVelocity[1], data.previousTangentialVelocity[2]);
            DEMproperties->getContactForce()->getplanewallSphereContactInformationList()[info.particleId1][info.particleId2] = info;
        }
        for (const auto &data : receivedspherecylinderwallcontactinformationData)
        {
            ContactInformation info;
            info.particleId1 = data.particleId1;
            info.particleId2 = data.particleId2;
            info.normalForce = Eigen::Vector3d(data.normalForce[0], data.normalForce[1], data.normalForce[2]);
            info.tangentialForce = Eigen::Vector3d(data.tangentialForce[0], data.tangentialForce[1], data.tangentialForce[2]);
            info.tangentialDisplacement = Eigen::Vector3d(data.tangentialDisplacement[0], data.tangentialDisplacement[1], data.tangentialDisplacement[2]);
            info.previousTangentialVelocity = Eigen::Vector3d(data.previousTangentialVelocity[0], data.previousTangentialVelocity[1], data.previousTangentialVelocity[2]);
            DEMproperties->getContactForce()->getcylinderwallSphereContactInformationList()[info.particleId1][info.particleId2] = info;
        }

        for (const auto &data : receivedfiberfibercontactinformationData)
        {
            ContactInformation info;
            info.particleId1 = data.particleId1;
            info.particleId2 = data.particleId2;
            info.normalForce = Eigen::Vector3d(data.normalForce[0], data.normalForce[1], data.normalForce[2]);
            info.tangentialForce = Eigen::Vector3d(data.tangentialForce[0], data.tangentialForce[1], data.tangentialForce[2]);
            info.tangentialDisplacement = Eigen::Vector3d(data.tangentialDisplacement[0], data.tangentialDisplacement[1], data.tangentialDisplacement[2]);
            info.previousTangentialVelocity = Eigen::Vector3d(data.previousTangentialVelocity[0], data.previousTangentialVelocity[1], data.previousTangentialVelocity[2]);
            DEMproperties->getContactForce()->getFiberFiberContactInformationList()[info.particleId1][info.particleId2] = info;
        }
        for (const auto &data : receivedfiberplanewallcontactinformationData)
        {
            ContactInformation info;
            info.particleId1 = data.particleId1;
            info.particleId2 = data.particleId2;
            info.normalForce = Eigen::Vector3d(data.normalForce[0], data.normalForce[1], data.normalForce[2]);
            info.tangentialForce = Eigen::Vector3d(data.tangentialForce[0], data.tangentialForce[1], data.tangentialForce[2]);
            info.tangentialDisplacement = Eigen::Vector3d(data.tangentialDisplacement[0], data.tangentialDisplacement[1], data.tangentialDisplacement[2]);
            info.previousTangentialVelocity = Eigen::Vector3d(data.previousTangentialVelocity[0], data.previousTangentialVelocity[1], data.previousTangentialVelocity[2]);
            DEMproperties->getContactForce()->getplanewallFiberContactInformationList()[info.particleId1][info.particleId2] = info;
        }
        for (const auto &data : receivedfibercylinderwallcontactinformationData)
        {
            ContactInformation info;
            info.particleId1 = data.particleId1;
            info.particleId2 = data.particleId2;
            info.normalForce = Eigen::Vector3d(data.normalForce[0], data.normalForce[1], data.normalForce[2]);
            info.tangentialForce = Eigen::Vector3d(data.tangentialForce[0], data.tangentialForce[1], data.tangentialForce[2]);
            info.tangentialDisplacement = Eigen::Vector3d(data.tangentialDisplacement[0], data.tangentialDisplacement[1], data.tangentialDisplacement[2]);
            info.previousTangentialVelocity = Eigen::Vector3d(data.previousTangentialVelocity[0], data.previousTangentialVelocity[1], data.previousTangentialVelocity[2]);
            DEMproperties->getContactForce()->getcylinderwallFiberContactInformationList()[info.particleId1][info.particleId2] = info;
        }
    }
}

void ParallelDEMModel::distributePlanewallToProcesses(int id)
{
    // PlanewallData planewallData;
    PlanewallData pw;
    if (world_rank == 0)
    {
        const auto &planewall = DEMproperties->getPlaneWall().at(id);
        pw.id = planewall->getId();
        pw.category = planewall->getType().getCategory();
        pw.subType = planewall->getType().getSubType();
        pw.state = planewall->getState();
        pw.normal[0] = planewall->getNormal()[0];
        pw.normal[1] = planewall->getNormal()[1];
        pw.normal[2] = planewall->getNormal()[2];

        pw.corner1[0] = planewall->getCorner1()[0];
        pw.corner1[1] = planewall->getCorner1()[1];
        pw.corner1[2] = planewall->getCorner1()[2];

        pw.corner2[0] = planewall->getCorner2()[0];
        pw.corner2[1] = planewall->getCorner2()[1];
        pw.corner2[2] = planewall->getCorner2()[2];

        pw.corner3[0] = planewall->getCorner3()[0];
        pw.corner3[1] = planewall->getCorner3()[1];
        pw.corner3[2] = planewall->getCorner3()[2];

        pw.velocity[0] = planewall->getVelocity()[0];
        pw.velocity[1] = planewall->getVelocity()[1];
        pw.velocity[2] = planewall->getVelocity()[2];

        pw.force[0] = planewall->getForce()[0];
        pw.force[1] = planewall->getForce()[1];
        pw.force[2] = planewall->getForce()[2];
    }
    MPI_Bcast(&pw, 1, mpi_planewallDataType, 0, MPI_COMM_WORLD);
    if (world_rank != 0)
    {
        std::unique_ptr<PlaneWall> planewalls;
        Eigen::Vector3d normal(pw.normal[0], pw.normal[1], pw.normal[2]);
        Eigen::Vector3d corner1(pw.corner1[0], pw.corner1[1], pw.corner1[2]);
        Eigen::Vector3d corner2(pw.corner2[0], pw.corner2[1], pw.corner2[2]);
        Eigen::Vector3d corner3(pw.corner3[0], pw.corner3[1], pw.corner3[2]);
        Eigen::Vector3d velocity(pw.velocity[0], pw.velocity[1], pw.velocity[2]);
        Eigen::Vector3d force(pw.force[0], pw.force[1], pw.force[2]);

        PropertyTypeID type = PropertyTypeID(pw.category, pw.subType);
        auto planewall = std::make_unique<PlaneWall>(pw.id, type, pw.state, normal, corner1, corner2, corner3, velocity, force);
        DEMproperties->setplanewall(planewall);
    }
}
void ParallelDEMModel::collectPlanewallfromProcesses(int id)
{
    std::vector<double> forcexvec(world_size);
    std::vector<double> forceyvec(world_size);
    std::vector<double> forcezvec(world_size);
    double localforcex = DEMproperties->getPlaneWall().at(id)->getForce().x();
    double localforcey = DEMproperties->getPlaneWall().at(id)->getForce().y();
    double localforcez = DEMproperties->getPlaneWall().at(id)->getForce().z();

    forcexvec[world_rank] = localforcex;
    forceyvec[world_rank] = localforcey;
    forcezvec[world_rank] = localforcez;

    MPI_Gather(&localforcex, 1, MPI_DOUBLE, forcexvec.data(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&localforcey, 1, MPI_DOUBLE, forceyvec.data(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&localforcez, 1, MPI_DOUBLE, forcezvec.data(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (world_rank == 0)
    {

        Eigen::Vector3d force = Eigen::Vector3d::Zero();
        for (int i = 1; i < forcexvec.size(); ++i)
        {
            force += Eigen::Vector3d(forcexvec[i], forceyvec[i], forcezvec[i]);
        }
        DEMproperties->getPlaneWall().at(id)->addForce(force);
    }
}

void ParallelDEMModel::neighborCommunication()
{
    ghostsphereparticles.clear();
    ghostfibershpereparticles.clear();
    ghostfiberbonds.clear();
    std::map<int, std::array<MPI_Request, 2>> sendSphereParticlesRequests;
    std::map<int, std::array<MPI_Request, 2>> sendFiberSphereParticlesRequests;
    std::map<int, std::array<MPI_Request, 2>> sendFiberBondsRequests;

    std::map<int, std::array<MPI_Request, 1>> receiveSphereParticlesRequests;
    std::map<int, std::array<MPI_Request, 1>> receiveFiberSphereParticlesRequests;
    std::map<int, std::array<MPI_Request, 1>> receiveFiberBondsRequests;

    std::map<int, std::vector<SphereParticleData>> recievesphereparticles;
    std::map<int, std::vector<SphereParticleData>> recievefibersphereparticles;
    std::map<int, std::vector<FiberBondData>> recievefiberbonds;

    std::map<int, std::vector<SphereParticleData>> sendsphereparticles;
    std::map<int, std::vector<SphereParticleData>> sendfibersphereparticles;
    std::map<int, std::vector<FiberBondData>> sendfiberbonds;

    std::map<int, int> sendsphereNumbers;
    std::map<int, int> sendfibersphereNumbers;
    std::map<int, int> sendfiberbondNumbers;
    // Initialize requests to MPI_REQUEST_NULL
    for (auto rankId : neighborRankId)
    {
        sendSphereParticlesRequests[rankId] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL};
        sendFiberSphereParticlesRequests[rankId] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL};
        sendFiberBondsRequests[rankId] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL};

        receiveSphereParticlesRequests[rankId] = {MPI_REQUEST_NULL};
        receiveFiberSphereParticlesRequests[rankId] = {MPI_REQUEST_NULL};
        receiveFiberBondsRequests[rankId] = {MPI_REQUEST_NULL};
    }

    auto &fibersphereparticles = DEMproperties->getfibersphereParticles();
    auto &sphereparticles = DEMproperties->getsphereParticles();
    auto &fiberbonds = DEMproperties->getfiberbonds();
    auto &gridbasedcontactdetection = DEMproperties->getGridBasedContactDetection();
    // gridbasedcontactdetection->assignParticleToGrid(sphereparticles, fiberbonds, fibersphereparticles);
    auto &gridSaveSphereParticles = gridbasedcontactdetection->getGridSaveSphereParticles();
    auto &gridSaveFiberBonds = gridbasedcontactdetection->getGridSaveFiberBonds();

    for (auto rankId : neighborRankId)
    {

        auto &sendsphereparticle = sendsphereparticles[rankId];
        auto &sendfibersphereparticle = sendfibersphereparticles[rankId];
        auto &sendfiberbond = sendfiberbonds[rankId];
        std::unordered_set<int> nodes;
        for (auto gridId : ghostLayer[rankId])
        {
            if (gridSaveSphereParticles.find(gridId) != gridSaveSphereParticles.end())
            {
                for (auto i : gridSaveSphereParticles.at(gridId))
                {
                    SphereParticleData sp;
                    sp.id = sphereparticles.at(i)->getId();
                    sp.category = sphereparticles.at(i)->getType().getCategory();
                    sp.subType = sphereparticles.at(i)->getType().getSubType();
                    sp.state = sphereparticles.at(i)->getState();

                    sp.position[0] = sphereparticles.at(i)->getPosition()[0];
                    sp.position[1] = sphereparticles.at(i)->getPosition()[1];
                    sp.position[2] = sphereparticles.at(i)->getPosition()[2];

                    sp.velocity[0] = sphereparticles.at(i)->getVelocity()[0];
                    sp.velocity[1] = sphereparticles.at(i)->getVelocity()[1];
                    sp.velocity[2] = sphereparticles.at(i)->getVelocity()[2];

                    sp.omega[0] = sphereparticles.at(i)->getOmega()[0];
                    sp.omega[1] = sphereparticles.at(i)->getOmega()[1];
                    sp.omega[2] = sphereparticles.at(i)->getOmega()[2];

                    sp.force[0] = sphereparticles.at(i)->getForce()[0];
                    sp.force[1] = sphereparticles.at(i)->getForce()[1];
                    sp.force[2] = sphereparticles.at(i)->getForce()[2];

                    sp.torque[0] = sphereparticles.at(i)->getTorque()[0];
                    sp.torque[1] = sphereparticles.at(i)->getTorque()[1];
                    sp.torque[2] = sphereparticles.at(i)->getTorque()[2];
                    sendsphereparticle.push_back(sp);
                }
            }
            if (gridSaveFiberBonds.find(gridId) != gridSaveFiberBonds.end())
            {
                for (auto i : gridSaveFiberBonds.at(gridId))
                {
                    FiberBondData fb;
                    nodes.insert(fiberbonds.at(i)->getNode1());
                    nodes.insert(fiberbonds.at(i)->getNode2());
                    fb.id = fiberbonds.at(i)->getId();
                    fb.category = fiberbonds.at(i)->getType().getCategory();
                    fb.subType = fiberbonds.at(i)->getType().getSubType();
                    fb.fiberid = fiberbonds.at(i)->getFiberId();
                    fb.node1 = fiberbonds.at(i)->getNode1();
                    fb.node2 = fiberbonds.at(i)->getNode2();
                    fb.neighborelement1 = fiberbonds.at(i)->getNeighborelement1();
                    fb.neighborelement2 = fiberbonds.at(i)->getNeighborelement2();
                    fb.energyDissipation = fiberbonds.at(i)->getEnergydissipation();
                    fb.tangentialforce[0] = fiberbonds.at(i)->getTangentialforce()[0];
                    fb.tangentialforce[1] = fiberbonds.at(i)->getTangentialforce()[1];
                    fb.tangentialforce[2] = fiberbonds.at(i)->getTangentialforce()[2];
                    fb.twisttorque[0] = fiberbonds.at(i)->getTwisttorque()[0];
                    fb.twisttorque[1] = fiberbonds.at(i)->getTwisttorque()[1];
                    fb.twisttorque[2] = fiberbonds.at(i)->getTwisttorque()[2];

                    fb.bendtorque[0] = fiberbonds.at(i)->getBendtorque()[0];
                    fb.bendtorque[1] = fiberbonds.at(i)->getBendtorque()[1];
                    fb.bendtorque[2] = fiberbonds.at(i)->getBendtorque()[2];
                    sendfiberbond.push_back(fb);
                }
            }
        }
        for (auto i : nodes)
        {
            SphereParticleData sp;
            sp.id = fibersphereparticles.at(i)->getId();
            sp.category = fibersphereparticles.at(i)->getType().getCategory();
            sp.subType = fibersphereparticles.at(i)->getType().getSubType();
            sp.state = fibersphereparticles.at(i)->getState();

            sp.position[0] = fibersphereparticles.at(i)->getPosition()[0];
            sp.position[1] = fibersphereparticles.at(i)->getPosition()[1];
            sp.position[2] = fibersphereparticles.at(i)->getPosition()[2];

            sp.velocity[0] = fibersphereparticles.at(i)->getVelocity()[0];
            sp.velocity[1] = fibersphereparticles.at(i)->getVelocity()[1];
            sp.velocity[2] = fibersphereparticles.at(i)->getVelocity()[2];

            sp.omega[0] = fibersphereparticles.at(i)->getOmega()[0];
            sp.omega[1] = fibersphereparticles.at(i)->getOmega()[1];
            sp.omega[2] = fibersphereparticles.at(i)->getOmega()[2];

            sp.force[0] = fibersphereparticles.at(i)->getForce()[0];
            sp.force[1] = fibersphereparticles.at(i)->getForce()[1];
            sp.force[2] = fibersphereparticles.at(i)->getForce()[2];

            sp.torque[0] = fibersphereparticles.at(i)->getTorque()[0];
            sp.torque[1] = fibersphereparticles.at(i)->getTorque()[1];
            sp.torque[2] = fibersphereparticles.at(i)->getTorque()[2];

            sendfibersphereparticle.push_back(sp);
        }

        int &sendsphereNumber = sendsphereNumbers[rankId];
        int &sendfibersphereNumber = sendfibersphereNumbers[rankId];
        int &sendfiberbondNumber = sendfiberbondNumbers[rankId];

        sendsphereNumber = sendsphereparticle.size();
        sendfibersphereNumber = sendfibersphereparticle.size();
        sendfiberbondNumber = sendfiberbond.size();
        // std::cout << "S" << world_rank << " " << rankId << " " << sendsphereNumber << std::endl;
        auto &req1 = sendSphereParticlesRequests[rankId];
        MPI_Isend(&sendsphereNumber, 1, MPI_INT, rankId, world_rank, MPI_COMM_WORLD, &req1[0]);
        if (sendsphereNumber > 0)
        {
            MPI_Isend(sendsphereparticle.data(), sendsphereNumber, mpi_sphereparticleDataType, rankId, world_rank, MPI_COMM_WORLD, &req1[1]);
        }
        auto &req2 = sendFiberSphereParticlesRequests[rankId];
        MPI_Isend(&sendfibersphereNumber, 1, MPI_INT, rankId, world_rank + 1, MPI_COMM_WORLD, &req2[0]);
        if (sendfibersphereNumber > 0)
        {
            MPI_Isend(sendfibersphereparticle.data(), sendfibersphereNumber, mpi_sphereparticleDataType, rankId, world_rank + 1, MPI_COMM_WORLD, &req2[1]);
        }
        auto &req3 = sendFiberBondsRequests[rankId];
        MPI_Isend(&sendfiberbondNumber, 1, MPI_INT, rankId, 2, MPI_COMM_WORLD, &req3[0]);
        // MPI_Send(&sendfiberbondNumber, 1, MPI_INT, rankId, world_rank + 2, MPI_COMM_WORLD);
        if (sendfiberbondNumber > 0)
        {
            MPI_Isend(sendfiberbond.data(), sendfiberbondNumber, mpi_fiberbondDataType, rankId, world_rank + 2, MPI_COMM_WORLD, &req3[1]);
        }
    }
    for (auto rankId : neighborRankId)
    {
        int receivesphereNumber = 0;
        int receivefiberbondNumber = 0;
        int receivefibersphereNumber = 0;
        auto &recievesphereparticle = recievesphereparticles[rankId];
        auto &recievefibersphereparticle = recievefibersphereparticles[rankId];
        auto &recievefiberbond = recievefiberbonds[rankId];

        MPI_Recv(&receivesphereNumber, 1, MPI_INT, rankId, rankId, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (receivesphereNumber > 0)
        {
            recievesphereparticle.resize(receivesphereNumber);
            MPI_Irecv(recievesphereparticle.data(), receivesphereNumber, mpi_sphereparticleDataType, rankId, rankId, MPI_COMM_WORLD, &receiveSphereParticlesRequests[rankId][0]);
        }

        MPI_Recv(&receivefibersphereNumber, 1, MPI_INT, rankId, rankId + 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (receivefibersphereNumber > 0)
        {
            recievefibersphereparticle.resize(receivefibersphereNumber);
            MPI_Irecv(recievefibersphereparticle.data(), receivefibersphereNumber, mpi_sphereparticleDataType, rankId, rankId + 1, MPI_COMM_WORLD, &receiveFiberSphereParticlesRequests[rankId][0]);
        }
        MPI_Recv(&receivefiberbondNumber, 1, MPI_INT, rankId, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        if (receivefiberbondNumber > 0)
        {
            recievefiberbond.resize(receivefiberbondNumber);
            MPI_Irecv(recievefiberbond.data(), receivefiberbondNumber, mpi_fiberbondDataType, rankId, rankId + 2, MPI_COMM_WORLD, &receiveFiberBondsRequests[rankId][0]);
        }
    }
    // Wait to make sure you actually send it out
    for (auto &kv : sendSphereParticlesRequests)
    {
        // auto &target_rank = kv.first;
        auto &reqs = kv.second;
        MPI_Wait(&reqs[0], MPI_STATUS_IGNORE);
        if (reqs[1] != MPI_REQUEST_NULL)
        { // Check if the second send was initialized
            MPI_Wait(&reqs[1], MPI_STATUS_IGNORE);
        }
    }
    for (auto &kv : sendFiberSphereParticlesRequests)
    {
        // auto &target_rank = kv.first;
        auto &reqs = kv.second;
        MPI_Wait(&reqs[0], MPI_STATUS_IGNORE);
        if (reqs[1] != MPI_REQUEST_NULL)
        { // Check if the second send was initialized
            MPI_Wait(&reqs[1], MPI_STATUS_IGNORE);
        }
    }
    for (auto &kv : sendFiberBondsRequests)
    {
        // auto &target_rank = kv.first;
        auto &reqs = kv.second;
        MPI_Wait(&reqs[0], MPI_STATUS_IGNORE);
        if (reqs[1] != MPI_REQUEST_NULL)
        { // Check if the second send was initialized
            MPI_Wait(&reqs[1], MPI_STATUS_IGNORE);
        }
    }

    // Wait to make sure you actually receive it
    for (auto &kv : receiveSphereParticlesRequests)
    {
        int src_rank = kv.first;
        auto &reqs = kv.second;
        if (reqs[0] != MPI_REQUEST_NULL)
        {
            MPI_Wait(&reqs[0], MPI_STATUS_IGNORE);
        }
        auto &recievesphereparticle = recievesphereparticles[src_rank];
        // std::cout << "R" << world_rank << " " << src_rank << " " << recievesphereparticle.size() << std::endl;

        // Put the nodes to the ghost_nodes
        for (const auto &data : recievesphereparticle)
        {
            Eigen::Vector3d position(data.position[0], data.position[1], data.position[2]);
            Eigen::Vector3d velocity(data.velocity[0], data.velocity[1], data.velocity[2]);
            Eigen::Vector3d omega(data.omega[0], data.omega[1], data.omega[2]);
            Eigen::Vector3d force(data.force[0], data.force[1], data.force[2]);
            Eigen::Vector3d torque(data.torque[0], data.torque[1], data.torque[2]);

            PropertyTypeID type = PropertyTypeID(data.category, data.subType);
            auto particle = std::make_unique<SphereParticle>(data.id, type, data.state, DEMproperties->getParticleManager(), position, velocity, omega, force, torque);
            ghostsphereparticles.insert({particle->getId(), std::move(particle)});
        }
    }
    for (auto &kv : receiveFiberSphereParticlesRequests)
    {
        int src_rank = kv.first;
        auto &reqs = kv.second;
        if (reqs[0] != MPI_REQUEST_NULL)
        {
            MPI_Wait(&reqs[0], MPI_STATUS_IGNORE);
        }
        auto &recievefibersphereparticle = recievefibersphereparticles[src_rank];

        for (const auto &data : recievefibersphereparticle)
        {
            Eigen::Vector3d position(data.position[0], data.position[1], data.position[2]);
            Eigen::Vector3d velocity(data.velocity[0], data.velocity[1], data.velocity[2]);
            Eigen::Vector3d omega(data.omega[0], data.omega[1], data.omega[2]);
            Eigen::Vector3d force(data.force[0], data.force[1], data.force[2]);
            Eigen::Vector3d torque(data.torque[0], data.torque[1], data.torque[2]);

            PropertyTypeID type = PropertyTypeID(data.category, data.subType);
            auto particle = std::make_unique<SphereParticle>(
                data.id, type, data.state, DEMproperties->getParticleManager(), position, velocity, omega, force, torque);
            ghostfibershpereparticles.insert({particle->getId(), std::move(particle)});
        }
    }

    for (auto &kv : receiveFiberBondsRequests)
    {
        int src_rank = kv.first;
        auto &reqs = kv.second;
        if (reqs[0] != MPI_REQUEST_NULL)
        {
            MPI_Wait(&reqs[0], MPI_STATUS_IGNORE);
        }
        auto &recievefiberbond = recievefiberbonds[src_rank];
        for (const auto &data : recievefiberbond)
        {
            Eigen::Vector3d tangentialforce(data.tangentialforce[0], data.tangentialforce[1], data.tangentialforce[2]);
            Eigen::Vector3d twisttorque(data.twisttorque[0], data.twisttorque[1], data.twisttorque[2]);
            Eigen::Vector3d bendtorque(data.bendtorque[0], data.bendtorque[1], data.bendtorque[2]);

            auto bond = std::make_unique<SphereCylinderBond>(
                data.id,
                PropertyTypeID(data.category, data.subType),
                DEMproperties->getParticleManager(),
                data.fiberid,
                data.node1,
                data.node2,
                data.neighborelement1,
                data.neighborelement2,
                data.energyDissipation);
            ghostfiberbonds.insert({bond->getId(), std::move(bond)});
        }
    }
}

void ParallelDEMModel::handleCollisions()
{
    auto &gridbasedcontactdetection = DEMproperties->getGridBasedContactDetection();
    auto &contactforce = DEMproperties->getContactForce();

    double timestep = DEMproperties->getTimestep();
    auto &fibersphereparticles = DEMproperties->getfibersphereParticles();
    auto &sphereparticles = DEMproperties->getsphereParticles();
    auto &fiberbonds = DEMproperties->getfiberbonds();
    std::unordered_map<int, std::unordered_set<int>> sphere_sphere_contact_paris;
    std::unordered_map<int, std::unordered_set<int>> sphere_fiber_contact_paris;
    std::unordered_map<int, std::unordered_set<int>> fiber_fiber_contact_paris;

    gridbasedcontactdetection->ghostParticleBroadPhase(ghostsphereparticles, ghostfiberbonds, ghostfibershpereparticles, fiberbonds, sphere_sphere_contact_paris,
                                                       sphere_fiber_contact_paris, fiber_fiber_contact_paris, neighborbonds);
    for (const auto &contact_list : sphere_sphere_contact_paris)
    {
        auto id1 = contact_list.first;
        for (auto id2 : contact_list.second)
        {
            if (ghostsphereparticles.find(id1) != ghostsphereparticles.end())
            {
                auto &sp1 = ghostsphereparticles[id1];
                auto &sp2 = sphereparticles.at(id2);

                contactforce->computeSphereSphereForce(sp1, sp2, timestep);
            }
            else
            {

                auto &sp1 = sphereparticles.at(id1);
                auto &sp2 = ghostsphereparticles[id2];
                contactforce->computeSphereSphereForce(sp1, sp2, timestep);
            }
        }
    }

    for (const auto &contact_list : sphere_fiber_contact_paris)
    {
        auto id1 = contact_list.first;
        for (auto id2 : contact_list.second)
        {
            if (ghostsphereparticles.find(id1) != ghostsphereparticles.end())
            {
                auto &sp1 = ghostsphereparticles[id1];
                auto &sc2 = fiberbonds.at(id2);

                auto &node1 = fibersphereparticles.at(sc2->getNode1());
                auto &node2 = fibersphereparticles.at(sc2->getNode2());

                contactforce->computeSphereFiberForce(sp1, sc2, node1, node2, timestep);
            }
            else
            {
                auto &sp1 = sphereparticles.at(id1);
                auto &sc2 = ghostfiberbonds[id2];

                auto &node1 = ghostfibershpereparticles[sc2->getNode1()];
                auto &node2 = ghostfibershpereparticles[sc2->getNode2()];

                contactforce->computeSphereFiberForce(sp1, sc2, node1, node2, timestep);
            }
        }
    }

    for (const auto &contact_list : fiber_fiber_contact_paris)
    {
        auto id1 = contact_list.first;
        for (auto id2 : contact_list.second)
        {
            if (ghostfiberbonds.find(id1) != ghostfiberbonds.end())
            {
                const auto &sc1 = ghostfiberbonds[id1];
                const auto &sc2 = fiberbonds.at(id2);

                auto &sconenode1 = ghostfibershpereparticles[sc1->getNode1()];
                auto &sconenode2 = ghostfibershpereparticles[sc1->getNode2()];
                auto &sctwonode1 = fibersphereparticles.at(sc2->getNode1());
                auto &sctwonode2 = fibersphereparticles.at(sc2->getNode2());

                contactforce->computeFiberFiberForce(sc1, sconenode1, sconenode2, sc2, sctwonode1, sctwonode2, timestep);
            }
            else
            {
                const auto &sc1 = fiberbonds.at(id1);
                const auto &sc2 = ghostfiberbonds[id2];

                auto &sconenode1 = fibersphereparticles.at(sc1->getNode1());
                auto &sconenode2 = fibersphereparticles.at(sc1->getNode2());
                auto &sctwonode1 = ghostfibershpereparticles[sc2->getNode1()];
                auto &sctwonode2 = ghostfibershpereparticles[sc2->getNode2()];
                contactforce->computeFiberFiberForce(sc1, sconenode1, sconenode2, sc2, sctwonode1, sctwonode2, timestep);
            }
        }
    }
}

void ParallelDEMModel::exchangeSharednodes()
{
    auto &fibersphereparticles = DEMproperties->getfibersphereParticles();
    auto &fiberbonds = DEMproperties->getfiberbonds();
    std::map<int, std::vector<SphereParticleData>> sendsharedsphereparticles;
    std::map<int, std::vector<SphereParticleData>> recievesharedsphereparticles;
    std::map<int, std::array<MPI_Request, 2>> sendsharedsphereparticlesRequests;
    std::map<int, std::array<MPI_Request, 1>> receivesharedsphereparticlesRequests;
    for (auto rankId : neighborRankId)
    {
        sendsharedsphereparticlesRequests[rankId] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL};
        receivesharedsphereparticlesRequests[rankId] = {MPI_REQUEST_NULL};
    }
    for (auto &pair : neighborbonds)
    {
        int id1 = pair.first;
        for (auto id2 : pair.second)
        {
            int shardId = ghostfiberbonds[id1]->getNode1();
            if (fiberbonds.at(id2)->getNode1() != shardId && fiberbonds.at(id2)->getNode2() != shardId)
            {
                shardId = ghostfiberbonds[id1]->getNode2();
            }

            auto &node1 = ghostfibershpereparticles.at(ghostfiberbonds[id1]->getNode1());
            auto &node2 = ghostfibershpereparticles.at(ghostfiberbonds[id1]->getNode2());

            double centerX = 0.5 * (node1->getPosition()[0] + node2->getPosition()[0]);
            double centerY = 0.5 * (node1->getPosition()[1] + node2->getPosition()[1]);
            double centerZ = 0.5 * (node1->getPosition()[2] + node2->getPosition()[2]);

            int cellX = static_cast<int>(centerX / gridSizeX);
            int cellY = static_cast<int>(centerY / gridSizeY);
            int cellZ = static_cast<int>(centerZ / gridSizeZ);

            int cellId = cellX + cellZ * numberOfGridX + cellY * numberOfGridX * numberOfGridZ;
            for (auto rankId : neighborRankId)
            {

                if (cellId >= offset[rankId * 2] && cellId <= offset[rankId * 2 + 1])
                {

                    auto &sendsharedsphereparticle = sendsharedsphereparticles[rankId];
                    SphereParticleData sp;
                    sp.id = fibersphereparticles.at(shardId)->getId();
                    sp.category = fibersphereparticles.at(shardId)->getType().getCategory();
                    sp.subType = fibersphereparticles.at(shardId)->getType().getSubType();
                    sp.state = fibersphereparticles.at(shardId)->getState();

                    sp.position[0] = fibersphereparticles.at(shardId)->getPosition()[0];
                    sp.position[1] = fibersphereparticles.at(shardId)->getPosition()[1];
                    sp.position[2] = fibersphereparticles.at(shardId)->getPosition()[2];

                    sp.velocity[0] = fibersphereparticles.at(shardId)->getVelocity()[0];
                    sp.velocity[1] = fibersphereparticles.at(shardId)->getVelocity()[1];
                    sp.velocity[2] = fibersphereparticles.at(shardId)->getVelocity()[2];

                    sp.omega[0] = fibersphereparticles.at(shardId)->getOmega()[0];
                    sp.omega[1] = fibersphereparticles.at(shardId)->getOmega()[1];
                    sp.omega[2] = fibersphereparticles.at(shardId)->getOmega()[2];

                    sp.force[0] = fibersphereparticles.at(shardId)->getForce()[0];
                    sp.force[1] = fibersphereparticles.at(shardId)->getForce()[1];
                    sp.force[2] = fibersphereparticles.at(shardId)->getForce()[2];

                    sp.torque[0] = fibersphereparticles.at(shardId)->getTorque()[0];
                    sp.torque[1] = fibersphereparticles.at(shardId)->getTorque()[1];
                    sp.torque[2] = fibersphereparticles.at(shardId)->getTorque()[2];
                    sendsharedsphereparticle.push_back(sp);
                    break;
                }
            }
        }
    }
    std::map<int, int> sendsphereNumbers;
    for (auto rankId : neighborRankId)
    {
        auto &sendsharedsphereparticle = sendsharedsphereparticles[rankId];
        int &sendsphereNumber = sendsphereNumbers[rankId];
        sendsphereNumber = sendsharedsphereparticle.size();

        auto &req = sendsharedsphereparticlesRequests[rankId];
        MPI_Isend(&sendsphereNumber, 1, MPI_INT, rankId, world_rank, MPI_COMM_WORLD, &req[0]);
        if (sendsphereNumber > 0)
        {
            MPI_Isend(sendsharedsphereparticle.data(), sendsphereNumber, mpi_sphereparticleDataType, rankId, world_rank, MPI_COMM_WORLD, &req[1]);
        }
    }
    for (auto rankId : neighborRankId)
    {
        int receivesphereNumber = 0;
        auto &recievesharedsphereparticle = recievesharedsphereparticles[rankId];
        MPI_Recv(&receivesphereNumber, 1, MPI_INT, rankId, rankId, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (receivesphereNumber > 0)
        {
            recievesharedsphereparticle.resize(receivesphereNumber);
            MPI_Irecv(recievesharedsphereparticle.data(), receivesphereNumber, mpi_sphereparticleDataType, rankId, rankId, MPI_COMM_WORLD, &receivesharedsphereparticlesRequests[rankId][0]);
        }
    }
    for (auto &kv : sendsharedsphereparticlesRequests)
    {
        // auto &target_rank = kv.first;
        auto &reqs = kv.second;
        MPI_Wait(&reqs[0], MPI_STATUS_IGNORE);
        if (reqs[1] != MPI_REQUEST_NULL)
        { // Check if the second send was initialized
            MPI_Wait(&reqs[1], MPI_STATUS_IGNORE);
        }
    }
    for (auto &kv : receivesharedsphereparticlesRequests)
    {
        int src_rank = kv.first;
        auto &reqs = kv.second;
        if (reqs[0] != MPI_REQUEST_NULL)
        {
            MPI_Wait(&reqs[0], MPI_STATUS_IGNORE);
        }
        auto &recievesharedsphereparticle = recievesharedsphereparticles[src_rank];
        for (const auto &data : recievesharedsphereparticle)
        {

            Eigen::Vector3d force(data.force[0], data.force[1], data.force[2]);
            Eigen::Vector3d torque(data.torque[0], data.torque[1], data.torque[2]);

            fibersphereparticles.at(data.id)->addForce(force);
            fibersphereparticles.at(data.id)->addTorque(torque);
            // std::cout << fibersphereparticles.at(data.id)->getVelocity().x() << " " << data.velocity[0] << std::endl;
        }
    }
}

void ParallelDEMModel::motion()
{
    std::map<int, std::array<MPI_Request, 2>> sendSphereParticlesRequests;
    std::map<int, std::array<MPI_Request, 2>> sendFiberSphereParticlesRequests;
    std::map<int, std::array<MPI_Request, 2>> sendFiberBondsRequests;

    std::map<int, std::array<MPI_Request, 2>> sendSphereSphereContactRequests;
    std::map<int, std::array<MPI_Request, 2>> sendSphereFiberContactRequests;
    std::map<int, std::array<MPI_Request, 2>> sendSpherePlanewallContactRequests;
    std::map<int, std::array<MPI_Request, 2>> sendSphereCylinderwallContactRequests;
    std::map<int, std::array<MPI_Request, 2>> sendFiberFiberContactRequests;
    std::map<int, std::array<MPI_Request, 2>> sendFiberPlanewallContactRequests;
    std::map<int, std::array<MPI_Request, 2>> sendFiberCylinderwallContactRequests;

    std::map<int, std::array<MPI_Request, 1>> receiveSphereParticlesRequests;
    std::map<int, std::array<MPI_Request, 1>> receiveFiberSphereParticlesRequests;
    std::map<int, std::array<MPI_Request, 1>> receiveFiberBondsRequests;

    std::map<int, std::array<MPI_Request, 1>> receiveSphereSphereContactRequests;
    std::map<int, std::array<MPI_Request, 1>> receiveSphereFiberContactRequests;
    std::map<int, std::array<MPI_Request, 1>> receiveSpherePlanewallContactRequests;
    std::map<int, std::array<MPI_Request, 1>> receiveSphereCylinderwallContactRequests;
    std::map<int, std::array<MPI_Request, 1>> receiveFiberFiberContactRequests;
    std::map<int, std::array<MPI_Request, 1>> receiveFiberPlanewallContactRequests;
    std::map<int, std::array<MPI_Request, 1>> receiveFiberCylinderwallContactRequests;

    std::map<int, std::vector<SphereParticleData>> receivesphereparticles;
    std::map<int, std::vector<SphereParticleData>> receivefibersphereparticles;
    std::map<int, std::vector<FiberBondData>> receivefiberbonds;

    std::map<int, std::vector<ContactInformationData>> receivespherespherecontacts;
    std::map<int, std::vector<ContactInformationData>> receivespherefibercontacts;
    std::map<int, std::vector<ContactInformationData>> receivesphereplanewallcontacts;
    std::map<int, std::vector<ContactInformationData>> receivespherecylinderwallcontacts;
    std::map<int, std::vector<ContactInformationData>> receivefiberfibercontacts;
    std::map<int, std::vector<ContactInformationData>> receivefiberplanewallcontacts;
    std::map<int, std::vector<ContactInformationData>> receivefibercylinderwallcontacts;

    std::map<int, std::vector<SphereParticleData>> sendsphereparticles;
    std::map<int, std::vector<SphereParticleData>> sendfibersphereparticles;
    std::map<int, std::vector<FiberBondData>> sendfiberbonds;

    std::map<int, std::vector<ContactInformationData>> sendspherespherecontacts;
    std::map<int, std::vector<ContactInformationData>> sendspherefibercontacts;
    std::map<int, std::vector<ContactInformationData>> sendsphereplanewallcontacts;
    std::map<int, std::vector<ContactInformationData>> sendspherecylinderwallcontacts;
    std::map<int, std::vector<ContactInformationData>> sendfiberfibercontacts;
    std::map<int, std::vector<ContactInformationData>> sendfiberplanewallcontacts;
    std::map<int, std::vector<ContactInformationData>> sendfibercylinderwallcontacts;
    for (auto rankId : neighborRankId)
    {
        sendSphereParticlesRequests[rankId] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL};
        sendFiberSphereParticlesRequests[rankId] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL};
        sendFiberBondsRequests[rankId] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL};

        sendSphereSphereContactRequests[rankId] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL};
        sendSphereFiberContactRequests[rankId] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL};
        sendSpherePlanewallContactRequests[rankId] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL};
        sendSphereCylinderwallContactRequests[rankId] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL};
        sendFiberFiberContactRequests[rankId] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL};
        sendFiberPlanewallContactRequests[rankId] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL};
        sendFiberCylinderwallContactRequests[rankId] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL};

        receiveSphereParticlesRequests[rankId] = {MPI_REQUEST_NULL};
        receiveFiberSphereParticlesRequests[rankId] = {MPI_REQUEST_NULL};
        receiveFiberBondsRequests[rankId] = {MPI_REQUEST_NULL};

        receiveSphereSphereContactRequests[rankId] = {MPI_REQUEST_NULL};
        receiveSphereFiberContactRequests[rankId] = {MPI_REQUEST_NULL};
        receiveSpherePlanewallContactRequests[rankId] = {MPI_REQUEST_NULL};
        receiveSphereCylinderwallContactRequests[rankId] = {MPI_REQUEST_NULL};
        receiveFiberFiberContactRequests[rankId] = {MPI_REQUEST_NULL};
        receiveFiberPlanewallContactRequests[rankId] = {MPI_REQUEST_NULL};
        receiveFiberCylinderwallContactRequests[rankId] = {MPI_REQUEST_NULL};
    }

    auto &sphereparticles = DEMproperties->getsphereParticles();
    auto &fiberbonds = DEMproperties->getfiberbonds();
    auto &fibersphereparticles = DEMproperties->getfibersphereParticles();

    auto &spheresphereContactInformationList = DEMproperties->getContactForce()->getSphereSphereContactInformationList();
    auto &spherefiberContactInformationList = DEMproperties->getContactForce()->getSphereFiberContactInformationList();
    auto &sphereplanewallContactInformationList = DEMproperties->getContactForce()->getplanewallSphereContactInformationList();
    auto &spherecylinderwallContactInformationList = DEMproperties->getContactForce()->getcylinderwallSphereContactInformationList();

    auto &fiberfiberContactInformationList = DEMproperties->getContactForce()->getFiberFiberContactInformationList();
    auto &fiberplanewallContactInformationList = DEMproperties->getContactForce()->getplanewallFiberContactInformationList();
    auto &fibercylinderwallContactInformationList = DEMproperties->getContactForce()->getcylinderwallFiberContactInformationList();
    std::unordered_map<int, std::unordered_set<int>> nodes;

    std::unordered_set<int> spherekeysToRemove;
    std::unordered_set<int> fiberbondskeysToRemove;
    std::unordered_set<int> fibersphereskeysToRemove;

    for (auto &sphere : sphereparticles)
    {
        int cellX = static_cast<int>(sphere.second->getPosition()[0] / gridSizeX);
        int cellY = static_cast<int>(sphere.second->getPosition()[1] / gridSizeY);
        int cellZ = static_cast<int>(sphere.second->getPosition()[2] / gridSizeZ);
        int cellId = cellX + cellZ * numberOfGridX + cellY * numberOfGridX * numberOfGridZ;
        for (auto rankId : neighborRankId)
        {
            if (cellId >= offset[rankId * 2] && cellId <= offset[rankId * 2 + 1])
            {

                spherekeysToRemove.insert(sphere.second->getId());
                SphereParticleData sp;
                sp.id = sphere.second->getId();
                sp.category = sphere.second->getType().getCategory();
                sp.subType = sphere.second->getType().getSubType();
                sp.state = sphere.second->getState();

                sp.position[0] = sphere.second->getPosition()[0];
                sp.position[1] = sphere.second->getPosition()[1];
                sp.position[2] = sphere.second->getPosition()[2];

                sp.velocity[0] = sphere.second->getVelocity()[0];
                sp.velocity[1] = sphere.second->getVelocity()[1];
                sp.velocity[2] = sphere.second->getVelocity()[2];

                sp.omega[0] = sphere.second->getOmega()[0];
                sp.omega[1] = sphere.second->getOmega()[1];
                sp.omega[2] = sphere.second->getOmega()[2];

                sp.force[0] = sphere.second->getForce()[0];
                sp.force[1] = sphere.second->getForce()[1];
                sp.force[2] = sphere.second->getForce()[2];

                sp.torque[0] = sphere.second->getTorque()[0];
                sp.torque[1] = sphere.second->getTorque()[1];
                sp.torque[2] = sphere.second->getTorque()[2];
                sendsphereparticles[rankId].push_back(sp);

                for (auto &contactInfor : spheresphereContactInformationList[sphere.second->getId()])
                {
                    ContactInformationData ssc;
                    ssc.particleId1 = contactInfor.second.particleId1;
                    ssc.particleId2 = contactInfor.second.particleId2;

                    ssc.previousTangentialVelocity[0] = contactInfor.second.previousTangentialVelocity.x();
                    ssc.previousTangentialVelocity[1] = contactInfor.second.previousTangentialVelocity.y();
                    ssc.previousTangentialVelocity[2] = contactInfor.second.previousTangentialVelocity.z();

                    ssc.tangentialDisplacement[0] = contactInfor.second.tangentialDisplacement.x();
                    ssc.tangentialDisplacement[1] = contactInfor.second.tangentialDisplacement.y();
                    ssc.tangentialDisplacement[2] = contactInfor.second.tangentialDisplacement.z();

                    ssc.normalForce[0] = contactInfor.second.normalForce.x();
                    ssc.normalForce[1] = contactInfor.second.normalForce.y();
                    ssc.normalForce[2] = contactInfor.second.normalForce.z();

                    ssc.tangentialForce[0] = contactInfor.second.tangentialForce.x();
                    ssc.tangentialForce[1] = contactInfor.second.tangentialForce.y();
                    ssc.tangentialForce[2] = contactInfor.second.tangentialForce.z();
                    sendspherespherecontacts[rankId].push_back(ssc);
                }
                for (auto &spherespherecontacts : spheresphereContactInformationList)
                {
                    if (spherespherecontacts.second.find(sphere.second->getId()) != spherespherecontacts.second.end())
                    {
                        int i = sphere.second->getId();
                        ContactInformationData ssc;
                        ssc.particleId1 = spherespherecontacts.second[i].particleId1;
                        ssc.particleId2 = spherespherecontacts.second[i].particleId2;

                        ssc.previousTangentialVelocity[0] = spherespherecontacts.second[i].previousTangentialVelocity.x();
                        ssc.previousTangentialVelocity[1] = spherespherecontacts.second[i].previousTangentialVelocity.y();
                        ssc.previousTangentialVelocity[2] = spherespherecontacts.second[i].previousTangentialVelocity.z();

                        ssc.tangentialDisplacement[0] = spherespherecontacts.second[i].tangentialDisplacement.x();
                        ssc.tangentialDisplacement[1] = spherespherecontacts.second[i].tangentialDisplacement.y();
                        ssc.tangentialDisplacement[2] = spherespherecontacts.second[i].tangentialDisplacement.z();

                        ssc.normalForce[0] = spherespherecontacts.second[i].normalForce.x();
                        ssc.normalForce[1] = spherespherecontacts.second[i].normalForce.y();
                        ssc.normalForce[2] = spherespherecontacts.second[i].normalForce.z();

                        ssc.tangentialForce[0] = spherespherecontacts.second[i].tangentialForce.x();
                        ssc.tangentialForce[1] = spherespherecontacts.second[i].tangentialForce.y();
                        ssc.tangentialForce[2] = spherespherecontacts.second[i].tangentialForce.z();
                        sendspherespherecontacts[rankId].push_back(ssc);
                    }
                }

                for (auto &contactInfor : spherefiberContactInformationList[sphere.second->getId()])
                {
                    ContactInformationData ssc;
                    ssc.particleId1 = contactInfor.second.particleId1;
                    ssc.particleId2 = contactInfor.second.particleId2;

                    ssc.previousTangentialVelocity[0] = contactInfor.second.previousTangentialVelocity.x();
                    ssc.previousTangentialVelocity[1] = contactInfor.second.previousTangentialVelocity.y();
                    ssc.previousTangentialVelocity[2] = contactInfor.second.previousTangentialVelocity.z();

                    ssc.tangentialDisplacement[0] = contactInfor.second.tangentialDisplacement.x();
                    ssc.tangentialDisplacement[1] = contactInfor.second.tangentialDisplacement.y();
                    ssc.tangentialDisplacement[2] = contactInfor.second.tangentialDisplacement.z();

                    ssc.normalForce[0] = contactInfor.second.normalForce.x();
                    ssc.normalForce[1] = contactInfor.second.normalForce.y();
                    ssc.normalForce[2] = contactInfor.second.normalForce.z();

                    ssc.tangentialForce[0] = contactInfor.second.tangentialForce.x();
                    ssc.tangentialForce[1] = contactInfor.second.tangentialForce.y();
                    ssc.tangentialForce[2] = contactInfor.second.tangentialForce.z();
                    sendspherefibercontacts[rankId].push_back(ssc);
                }
                for (auto &sphereplanewllcontacts : sphereplanewallContactInformationList)
                {
                    if (sphereplanewllcontacts.second.find(sphere.second->getId()) != sphereplanewllcontacts.second.end())
                    {
                        int i = sphere.second->getId();
                        ContactInformationData ssc;
                        ssc.particleId1 = sphereplanewllcontacts.second[i].particleId1;
                        ssc.particleId2 = sphereplanewllcontacts.second[i].particleId2;

                        ssc.previousTangentialVelocity[0] = sphereplanewllcontacts.second[i].previousTangentialVelocity.x();
                        ssc.previousTangentialVelocity[1] = sphereplanewllcontacts.second[i].previousTangentialVelocity.y();
                        ssc.previousTangentialVelocity[2] = sphereplanewllcontacts.second[i].previousTangentialVelocity.z();

                        ssc.tangentialDisplacement[0] = sphereplanewllcontacts.second[i].tangentialDisplacement.x();
                        ssc.tangentialDisplacement[1] = sphereplanewllcontacts.second[i].tangentialDisplacement.y();
                        ssc.tangentialDisplacement[2] = sphereplanewllcontacts.second[i].tangentialDisplacement.z();

                        ssc.normalForce[0] = sphereplanewllcontacts.second[i].normalForce.x();
                        ssc.normalForce[1] = sphereplanewllcontacts.second[i].normalForce.y();
                        ssc.normalForce[2] = sphereplanewllcontacts.second[i].normalForce.z();

                        ssc.tangentialForce[0] = sphereplanewllcontacts.second[i].tangentialForce.x();
                        ssc.tangentialForce[1] = sphereplanewllcontacts.second[i].tangentialForce.y();
                        ssc.tangentialForce[2] = sphereplanewllcontacts.second[i].tangentialForce.z();
                        sendsphereplanewallcontacts[rankId].push_back(ssc);
                    }
                }
                for (auto &spherecylinderwllcontacts : spherecylinderwallContactInformationList)
                {
                    if (spherecylinderwllcontacts.second.find(sphere.second->getId()) != spherecylinderwllcontacts.second.end())
                    {
                        int i = sphere.second->getId();
                        ContactInformationData ssc;
                        ssc.particleId1 = spherecylinderwllcontacts.second[i].particleId1;
                        ssc.particleId2 = spherecylinderwllcontacts.second[i].particleId2;

                        ssc.previousTangentialVelocity[0] = spherecylinderwllcontacts.second[i].previousTangentialVelocity.x();
                        ssc.previousTangentialVelocity[1] = spherecylinderwllcontacts.second[i].previousTangentialVelocity.y();
                        ssc.previousTangentialVelocity[2] = spherecylinderwllcontacts.second[i].previousTangentialVelocity.z();

                        ssc.tangentialDisplacement[0] = spherecylinderwllcontacts.second[i].tangentialDisplacement.x();
                        ssc.tangentialDisplacement[1] = spherecylinderwllcontacts.second[i].tangentialDisplacement.y();
                        ssc.tangentialDisplacement[2] = spherecylinderwllcontacts.second[i].tangentialDisplacement.z();

                        ssc.normalForce[0] = spherecylinderwllcontacts.second[i].normalForce.x();
                        ssc.normalForce[1] = spherecylinderwllcontacts.second[i].normalForce.y();
                        ssc.normalForce[2] = spherecylinderwllcontacts.second[i].normalForce.z();

                        ssc.tangentialForce[0] = spherecylinderwllcontacts.second[i].tangentialForce.x();
                        ssc.tangentialForce[1] = spherecylinderwllcontacts.second[i].tangentialForce.y();
                        ssc.tangentialForce[2] = spherecylinderwllcontacts.second[i].tangentialForce.z();
                        sendspherecylinderwallcontacts[rankId].push_back(ssc);
                    }
                }

                break;
            }
        }
    }

    for (auto &bond : fiberbonds)
    {
        auto &node1 = fibersphereparticles.at(bond.second->getNode1());
        auto &node2 = fibersphereparticles.at(bond.second->getNode2());
        double centerX = 0.5 * (node1->getPosition()[0] + node2->getPosition()[0]);
        double centerY = 0.5 * (node1->getPosition()[1] + node2->getPosition()[1]);
        double centerZ = 0.5 * (node1->getPosition()[2] + node2->getPosition()[2]);

        int cellX = static_cast<int>(centerX / gridSizeX);
        int cellY = static_cast<int>(centerY / gridSizeY);
        int cellZ = static_cast<int>(centerZ / gridSizeZ);

        int cellId = cellX + cellZ * numberOfGridX + cellY * numberOfGridX * numberOfGridZ;
        for (auto rankId : neighborRankId)
        {

            if (cellId >= offset[rankId * 2] && cellId <= offset[rankId * 2 + 1])
            {
                fiberbondskeysToRemove.insert(bond.second->getId());
                FiberBondData fb;
                nodes[rankId].insert(bond.second->getNode1());
                nodes[rankId].insert(bond.second->getNode2());
                fb.id = bond.second->getId();
                fb.category = bond.second->getType().getCategory();
                fb.subType = bond.second->getType().getSubType();
                fb.fiberid = bond.second->getFiberId();
                fb.node1 = bond.second->getNode1();
                fb.node2 = bond.second->getNode2();
                fb.neighborelement1 = bond.second->getNeighborelement1();
                fb.neighborelement2 = bond.second->getNeighborelement2();
                fb.energyDissipation = bond.second->getEnergydissipation();
                fb.tangentialforce[0] = bond.second->getTangentialforce()[0];
                fb.tangentialforce[1] = bond.second->getTangentialforce()[1];
                fb.tangentialforce[2] = bond.second->getTangentialforce()[2];
                fb.twisttorque[0] = bond.second->getTwisttorque()[0];
                fb.twisttorque[1] = bond.second->getTwisttorque()[1];
                fb.twisttorque[2] = bond.second->getTwisttorque()[2];

                fb.bendtorque[0] = bond.second->getBendtorque()[0];
                fb.bendtorque[1] = bond.second->getBendtorque()[1];
                fb.bendtorque[2] = bond.second->getBendtorque()[2];
                sendfiberbonds[rankId].push_back(fb);
                for (auto &contactInfor : fiberfiberContactInformationList[bond.second->getId()])
                {
                    ContactInformationData ssc;
                    ssc.particleId1 = contactInfor.second.particleId1;
                    ssc.particleId2 = contactInfor.second.particleId2;

                    ssc.previousTangentialVelocity[0] = contactInfor.second.previousTangentialVelocity.x();
                    ssc.previousTangentialVelocity[1] = contactInfor.second.previousTangentialVelocity.y();
                    ssc.previousTangentialVelocity[2] = contactInfor.second.previousTangentialVelocity.z();

                    ssc.tangentialDisplacement[0] = contactInfor.second.tangentialDisplacement.x();
                    ssc.tangentialDisplacement[1] = contactInfor.second.tangentialDisplacement.y();
                    ssc.tangentialDisplacement[2] = contactInfor.second.tangentialDisplacement.z();

                    ssc.normalForce[0] = contactInfor.second.normalForce.x();
                    ssc.normalForce[1] = contactInfor.second.normalForce.y();
                    ssc.normalForce[2] = contactInfor.second.normalForce.z();

                    ssc.tangentialForce[0] = contactInfor.second.tangentialForce.x();
                    ssc.tangentialForce[1] = contactInfor.second.tangentialForce.y();
                    ssc.tangentialForce[2] = contactInfor.second.tangentialForce.z();
                    sendfiberfibercontacts[rankId].push_back(ssc);
                }
                for (auto &fiberfibercontacts : fiberfiberContactInformationList)
                {
                    if (fiberfibercontacts.second.find(bond.second->getId()) != fiberfibercontacts.second.end())
                    {
                        int i = bond.second->getId();
                        ContactInformationData ssc;
                        ssc.particleId1 = fiberfibercontacts.second[i].particleId1;
                        ssc.particleId2 = fiberfibercontacts.second[i].particleId2;

                        ssc.previousTangentialVelocity[0] = fiberfibercontacts.second[i].previousTangentialVelocity.x();
                        ssc.previousTangentialVelocity[1] = fiberfibercontacts.second[i].previousTangentialVelocity.y();
                        ssc.previousTangentialVelocity[2] = fiberfibercontacts.second[i].previousTangentialVelocity.z();

                        ssc.tangentialDisplacement[0] = fiberfibercontacts.second[i].tangentialDisplacement.x();
                        ssc.tangentialDisplacement[1] = fiberfibercontacts.second[i].tangentialDisplacement.y();
                        ssc.tangentialDisplacement[2] = fiberfibercontacts.second[i].tangentialDisplacement.z();

                        ssc.normalForce[0] = fiberfibercontacts.second[i].normalForce.x();
                        ssc.normalForce[1] = fiberfibercontacts.second[i].normalForce.y();
                        ssc.normalForce[2] = fiberfibercontacts.second[i].normalForce.z();

                        ssc.tangentialForce[0] = fiberfibercontacts.second[i].tangentialForce.x();
                        ssc.tangentialForce[1] = fiberfibercontacts.second[i].tangentialForce.y();
                        ssc.tangentialForce[2] = fiberfibercontacts.second[i].tangentialForce.z();
                        sendfiberfibercontacts[rankId].push_back(ssc);
                    }
                }
                for (auto &spherefibercontacts : spherefiberContactInformationList)
                {
                    if (spherefibercontacts.second.find(bond.second->getId()) != spherefibercontacts.second.end())
                    {
                        int i = bond.second->getId();
                        ContactInformationData ssc;
                        ssc.particleId1 = spherefibercontacts.second[i].particleId1;
                        ssc.particleId2 = spherefibercontacts.second[i].particleId2;

                        ssc.previousTangentialVelocity[0] = spherefibercontacts.second[i].previousTangentialVelocity.x();
                        ssc.previousTangentialVelocity[1] = spherefibercontacts.second[i].previousTangentialVelocity.y();
                        ssc.previousTangentialVelocity[2] = spherefibercontacts.second[i].previousTangentialVelocity.z();

                        ssc.tangentialDisplacement[0] = spherefibercontacts.second[i].tangentialDisplacement.x();
                        ssc.tangentialDisplacement[1] = spherefibercontacts.second[i].tangentialDisplacement.y();
                        ssc.tangentialDisplacement[2] = spherefibercontacts.second[i].tangentialDisplacement.z();

                        ssc.normalForce[0] = spherefibercontacts.second[i].normalForce.x();
                        ssc.normalForce[1] = spherefibercontacts.second[i].normalForce.y();
                        ssc.normalForce[2] = spherefibercontacts.second[i].normalForce.z();

                        ssc.tangentialForce[0] = spherefibercontacts.second[i].tangentialForce.x();
                        ssc.tangentialForce[1] = spherefibercontacts.second[i].tangentialForce.y();
                        ssc.tangentialForce[2] = spherefibercontacts.second[i].tangentialForce.z();
                        sendspherefibercontacts[rankId].push_back(ssc);
                    }
                }
                for (auto &fiberplanewllcontacts : fiberplanewallContactInformationList)
                {
                    if (fiberplanewllcontacts.second.find(bond.second->getId()) != fiberplanewllcontacts.second.end())
                    {
                        int i = bond.second->getId();
                        ContactInformationData ssc;
                        ssc.particleId1 = fiberplanewllcontacts.second[i].particleId1;
                        ssc.particleId2 = fiberplanewllcontacts.second[i].particleId2;

                        ssc.previousTangentialVelocity[0] = fiberplanewllcontacts.second[i].previousTangentialVelocity.x();
                        ssc.previousTangentialVelocity[1] = fiberplanewllcontacts.second[i].previousTangentialVelocity.y();
                        ssc.previousTangentialVelocity[2] = fiberplanewllcontacts.second[i].previousTangentialVelocity.z();

                        ssc.tangentialDisplacement[0] = fiberplanewllcontacts.second[i].tangentialDisplacement.x();
                        ssc.tangentialDisplacement[1] = fiberplanewllcontacts.second[i].tangentialDisplacement.y();
                        ssc.tangentialDisplacement[2] = fiberplanewllcontacts.second[i].tangentialDisplacement.z();

                        ssc.normalForce[0] = fiberplanewllcontacts.second[i].normalForce.x();
                        ssc.normalForce[1] = fiberplanewllcontacts.second[i].normalForce.y();
                        ssc.normalForce[2] = fiberplanewllcontacts.second[i].normalForce.z();

                        ssc.tangentialForce[0] = fiberplanewllcontacts.second[i].tangentialForce.x();
                        ssc.tangentialForce[1] = fiberplanewllcontacts.second[i].tangentialForce.y();
                        ssc.tangentialForce[2] = fiberplanewllcontacts.second[i].tangentialForce.z();
                        sendfiberplanewallcontacts[rankId].push_back(ssc);
                    }
                }
                for (auto &fibercylinderwllcontacts : fibercylinderwallContactInformationList)
                {
                    if (fibercylinderwllcontacts.second.find(bond.second->getId()) != fibercylinderwllcontacts.second.end())
                    {
                        int i = bond.second->getId();
                        ContactInformationData ssc;
                        ssc.particleId1 = fibercylinderwllcontacts.second[i].particleId1;
                        ssc.particleId2 = fibercylinderwllcontacts.second[i].particleId2;

                        ssc.previousTangentialVelocity[0] = fibercylinderwllcontacts.second[i].previousTangentialVelocity.x();
                        ssc.previousTangentialVelocity[1] = fibercylinderwllcontacts.second[i].previousTangentialVelocity.y();
                        ssc.previousTangentialVelocity[2] = fibercylinderwllcontacts.second[i].previousTangentialVelocity.z();

                        ssc.tangentialDisplacement[0] = fibercylinderwllcontacts.second[i].tangentialDisplacement.x();
                        ssc.tangentialDisplacement[1] = fibercylinderwllcontacts.second[i].tangentialDisplacement.y();
                        ssc.tangentialDisplacement[2] = fibercylinderwllcontacts.second[i].tangentialDisplacement.z();

                        ssc.normalForce[0] = fibercylinderwllcontacts.second[i].normalForce.x();
                        ssc.normalForce[1] = fibercylinderwllcontacts.second[i].normalForce.y();
                        ssc.normalForce[2] = fibercylinderwllcontacts.second[i].normalForce.z();

                        ssc.tangentialForce[0] = fibercylinderwllcontacts.second[i].tangentialForce.x();
                        ssc.tangentialForce[1] = fibercylinderwllcontacts.second[i].tangentialForce.y();
                        ssc.tangentialForce[2] = fibercylinderwllcontacts.second[i].tangentialForce.z();
                        sendfibercylinderwallcontacts[rankId].push_back(ssc);
                    }
                }
                break;
            }
        }
    }

    for (auto &node : nodes)
    {
        for (auto &i : node.second)
        {
            SphereParticleData sp;
            sp.id = fibersphereparticles.at(i)->getId();
            sp.category = fibersphereparticles.at(i)->getType().getCategory();
            sp.subType = fibersphereparticles.at(i)->getType().getSubType();
            sp.state = fibersphereparticles.at(i)->getState();

            sp.position[0] = fibersphereparticles.at(i)->getPosition()[0];
            sp.position[1] = fibersphereparticles.at(i)->getPosition()[1];
            sp.position[2] = fibersphereparticles.at(i)->getPosition()[2];

            sp.velocity[0] = fibersphereparticles.at(i)->getVelocity()[0];
            sp.velocity[1] = fibersphereparticles.at(i)->getVelocity()[1];
            sp.velocity[2] = fibersphereparticles.at(i)->getVelocity()[2];

            sp.omega[0] = fibersphereparticles.at(i)->getOmega()[0];
            sp.omega[1] = fibersphereparticles.at(i)->getOmega()[1];
            sp.omega[2] = fibersphereparticles.at(i)->getOmega()[2];

            sp.force[0] = fibersphereparticles.at(i)->getForce()[0];
            sp.force[1] = fibersphereparticles.at(i)->getForce()[1];
            sp.force[2] = fibersphereparticles.at(i)->getForce()[2];

            sp.torque[0] = fibersphereparticles.at(i)->getTorque()[0];
            sp.torque[1] = fibersphereparticles.at(i)->getTorque()[1];
            sp.torque[2] = fibersphereparticles.at(i)->getTorque()[2];

            sendfibersphereparticles[node.first].push_back(sp);
        }
    }

    std::map<int, int> sendsphereNumbers;
    std::map<int, int> sendfibersphereNumbers;
    std::map<int, int> sendfiberbondNumbers;

    std::map<int, int> sendspherespherecontactNumbers;
    std::map<int, int> sendspherefibercontactNumbers;
    std::map<int, int> sendsphereplanewallcontactNumbers;
    std::map<int, int> sendspherecylinderwallcontactNumbers;
    std::map<int, int> sendfiberfibercontactNumbers;
    std::map<int, int> sendfiberplanewallcontactNumbers;
    std::map<int, int> sendfibercylinderwallcontactNumbers;
    for (auto rankId : neighborRankId)
    {
        auto &sendsphereparticle = sendsphereparticles[rankId];
        auto &sendfibersphereparticle = sendfibersphereparticles[rankId];
        auto &sendfiberbond = sendfiberbonds[rankId];

        auto &sendspherespherecontact = sendspherespherecontacts[rankId];
        auto &sendspherefibercontact = sendspherefibercontacts[rankId];
        auto &sendsphereplanewallcontact = sendsphereplanewallcontacts[rankId];
        auto &sendspherecylinderwallcontact = sendspherecylinderwallcontacts[rankId];
        auto &sendfiberfibercontact = sendfiberfibercontacts[rankId];
        auto &sendfiberplanewallcontact = sendfiberplanewallcontacts[rankId];
        auto &sendfibercylinderwallcontact = sendfibercylinderwallcontacts[rankId];

        int &sendsphereNumber = sendsphereNumbers[rankId];
        int &sendfibersphereNumber = sendfibersphereNumbers[rankId];
        int &sendfiberbondNumber = sendfiberbondNumbers[rankId];

        int &sendspherespherecontactNumber = sendspherespherecontactNumbers[rankId];
        int &sendspherefibercontactNumber = sendspherefibercontactNumbers[rankId];
        int &sendsphereplanewallcontactNumber = sendsphereplanewallcontactNumbers[rankId];
        int &sendspherecylinderwallcontactNumber = sendspherecylinderwallcontactNumbers[rankId];
        int &sendfiberfibercontactNumber = sendfiberfibercontactNumbers[rankId];
        int &sendfiberplanewallcontactNumber = sendfiberplanewallcontactNumbers[rankId];
        int &sendfibercylinderwallcontactNumber = sendfibercylinderwallcontactNumbers[rankId];

        sendsphereNumber = sendsphereparticle.size();
        sendfibersphereNumber = sendfibersphereparticle.size();
        sendfiberbondNumber = sendfiberbond.size();

        sendspherespherecontactNumber = sendspherespherecontact.size();
        sendspherefibercontactNumber = sendspherefibercontact.size();
        sendsphereplanewallcontactNumber = sendsphereplanewallcontact.size();
        sendspherecylinderwallcontactNumber = sendspherecylinderwallcontact.size();
        sendfiberfibercontactNumber = sendfiberfibercontact.size();
        sendfiberplanewallcontactNumber = sendfiberplanewallcontact.size();
        sendfibercylinderwallcontactNumber = sendfibercylinderwallcontact.size();

        auto &req1 = sendSphereParticlesRequests[rankId];
        MPI_Isend(&sendsphereNumber, 1, MPI_INT, rankId, world_rank, MPI_COMM_WORLD, &req1[0]);
        if (sendsphereNumber > 0)
        {
            MPI_Isend(sendsphereparticle.data(), sendsphereNumber, mpi_sphereparticleDataType, rankId, world_rank, MPI_COMM_WORLD, &req1[1]);
        }
        auto &req2 = sendFiberSphereParticlesRequests[rankId];
        MPI_Isend(&sendfibersphereNumber, 1, MPI_INT, rankId, world_rank + 1, MPI_COMM_WORLD, &req2[0]);
        if (sendfibersphereNumber > 0)
        {
            MPI_Isend(sendfibersphereparticle.data(), sendfibersphereNumber, mpi_sphereparticleDataType, rankId, world_rank + 1, MPI_COMM_WORLD, &req2[1]);
        }
        auto &req3 = sendFiberBondsRequests[rankId];
        MPI_Isend(&sendfiberbondNumber, 1, MPI_INT, rankId, world_rank + 2, MPI_COMM_WORLD, &req3[0]);
        if (sendfiberbondNumber > 0)
        {
            MPI_Isend(sendfiberbond.data(), sendfiberbondNumber, mpi_fiberbondDataType, rankId, world_rank + 2, MPI_COMM_WORLD, &req3[1]);
        }

        auto &req4 = sendSphereSphereContactRequests[rankId];
        MPI_Isend(&sendspherespherecontactNumber, 1, MPI_INT, rankId, world_rank + 3, MPI_COMM_WORLD, &req4[0]);
        if (sendspherespherecontactNumber > 0)
        {
            MPI_Isend(sendspherespherecontact.data(), sendspherespherecontactNumber, mpi_contactinformationDataType, rankId, world_rank + 3, MPI_COMM_WORLD, &req4[1]);
        }

        auto &req5 = sendSphereFiberContactRequests[rankId];
        MPI_Isend(&sendspherefibercontactNumber, 1, MPI_INT, rankId, world_rank + 4, MPI_COMM_WORLD, &req5[0]);
        if (sendspherefibercontactNumber > 0)
        {
            MPI_Isend(sendspherefibercontact.data(), sendspherefibercontactNumber, mpi_contactinformationDataType, rankId, world_rank + 4, MPI_COMM_WORLD, &req5[1]);
        }
        auto &req6 = sendSpherePlanewallContactRequests[rankId];
        MPI_Isend(&sendsphereplanewallcontactNumber, 1, MPI_INT, rankId, world_rank + 5, MPI_COMM_WORLD, &req6[0]);
        if (sendsphereplanewallcontactNumber > 0)
        {
            MPI_Isend(sendsphereplanewallcontact.data(), sendsphereplanewallcontactNumber, mpi_contactinformationDataType, rankId, world_rank + 5, MPI_COMM_WORLD, &req6[1]);
        }
        auto &req7 = sendSphereCylinderwallContactRequests[rankId];
        MPI_Isend(&sendspherecylinderwallcontactNumber, 1, MPI_INT, rankId, world_rank + 6, MPI_COMM_WORLD, &req7[0]);
        if (sendspherecylinderwallcontactNumber > 0)
        {
            MPI_Isend(sendspherecylinderwallcontact.data(), sendspherecylinderwallcontactNumber, mpi_contactinformationDataType, rankId, world_rank + 6, MPI_COMM_WORLD, &req7[1]);
        }

        auto &req8 = sendFiberFiberContactRequests[rankId];
        MPI_Isend(&sendfiberfibercontactNumber, 1, MPI_INT, rankId, world_rank + 7, MPI_COMM_WORLD, &req8[0]);
        if (sendfiberfibercontactNumber > 0)
        {
            MPI_Isend(sendfiberfibercontact.data(), sendfiberfibercontactNumber, mpi_contactinformationDataType, rankId, world_rank + 7, MPI_COMM_WORLD, &req8[1]);
        }
        auto &req9 = sendFiberPlanewallContactRequests[rankId];
        MPI_Isend(&sendfiberplanewallcontactNumber, 1, MPI_INT, rankId, world_rank + 8, MPI_COMM_WORLD, &req9[0]);
        if (sendfiberplanewallcontactNumber > 0)
        {
            MPI_Isend(sendfiberplanewallcontact.data(), sendfiberplanewallcontactNumber, mpi_contactinformationDataType, rankId, world_rank + 8, MPI_COMM_WORLD, &req9[1]);
        }
        auto &req10 = sendFiberCylinderwallContactRequests[rankId];
        MPI_Isend(&sendfibercylinderwallcontactNumber, 1, MPI_INT, rankId, world_rank + 9, MPI_COMM_WORLD, &req10[0]);
        if (sendfibercylinderwallcontactNumber > 0)
        {
            MPI_Isend(sendfibercylinderwallcontact.data(), sendfibercylinderwallcontactNumber, mpi_contactinformationDataType, rankId, world_rank + 9, MPI_COMM_WORLD, &req10[1]);
        }
    }
    for (auto rankId : neighborRankId)
    {
        int receivesphereNumber = 0;
        int receivefiberbondNumber = 0;
        int receivefibersphereNumber = 0;
        int receivespherespherecontactNumber = 0;
        int receivespherefibercontactNumber = 0;
        int receivesphereplanewallcontactNumber = 0;
        int receivespherecylinderwallcontactNumber = 0;
        int receivefiberfibercontactNumber = 0;
        int receivefiberplanewallcontactNumber = 0;
        int receivefibercylinderwallcontactNumber = 0;

        auto &receivesphereparticle = receivesphereparticles[rankId];
        auto &receivefibersphereparticle = receivefibersphereparticles[rankId];
        auto &receivefiberbond = receivefiberbonds[rankId];

        auto &receivespherespherecontact = receivespherespherecontacts[rankId];
        auto &receivespherefibercontact = receivespherefibercontacts[rankId];
        auto &receivesphereplanewallcontact = receivesphereplanewallcontacts[rankId];
        auto &receivespherecylinderwallcontact = receivespherecylinderwallcontacts[rankId];
        auto &receivefiberfibercontact = receivefiberfibercontacts[rankId];
        auto &receivefiberplanewallcontact = receivefiberplanewallcontacts[rankId];
        auto &receivefibercylinderwallcontact = receivefibercylinderwallcontacts[rankId];

        MPI_Recv(&receivesphereNumber, 1, MPI_INT, rankId, rankId, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (receivesphereNumber > 0)
        {
            receivesphereparticle.resize(receivesphereNumber);
            MPI_Irecv(receivesphereparticle.data(), receivesphereNumber, mpi_sphereparticleDataType, rankId, rankId, MPI_COMM_WORLD, &receiveSphereParticlesRequests[rankId][0]);
        }
        MPI_Recv(&receivefibersphereNumber, 1, MPI_INT, rankId, rankId + 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (receivefibersphereNumber > 0)
        {

            receivefibersphereparticle.resize(receivefibersphereNumber);
            MPI_Irecv(receivefibersphereparticle.data(), receivefibersphereNumber, mpi_sphereparticleDataType, rankId, rankId + 1, MPI_COMM_WORLD, &receiveFiberSphereParticlesRequests[rankId][0]);
        }
        MPI_Recv(&receivefiberbondNumber, 1, MPI_INT, rankId, rankId + 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (receivefiberbondNumber > 0)
        {
            receivefiberbond.resize(receivefiberbondNumber);
            MPI_Irecv(receivefiberbond.data(), receivefiberbondNumber, mpi_fiberbondDataType, rankId, rankId + 2, MPI_COMM_WORLD, &receiveFiberBondsRequests[rankId][0]);
        }
        MPI_Recv(&receivespherespherecontactNumber, 1, MPI_INT, rankId, rankId + 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (receivespherespherecontactNumber > 0)
        {
            receivespherespherecontact.resize(receivespherespherecontactNumber);
            MPI_Irecv(receivespherespherecontact.data(), receivespherespherecontactNumber, mpi_contactinformationDataType, rankId,
                      rankId + 3, MPI_COMM_WORLD, &receiveSphereSphereContactRequests[rankId][0]);
        }
        MPI_Recv(&receivespherefibercontactNumber, 1, MPI_INT, rankId, rankId + 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (receivespherefibercontactNumber > 0)
        {
            receivespherefibercontact.resize(receivespherefibercontactNumber);
            MPI_Irecv(receivespherefibercontact.data(), receivespherefibercontactNumber, mpi_contactinformationDataType, rankId,
                      rankId + 4, MPI_COMM_WORLD, &receiveSphereFiberContactRequests[rankId][0]);
        }
        MPI_Recv(&receivesphereplanewallcontactNumber, 1, MPI_INT, rankId, rankId + 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (receivesphereplanewallcontactNumber > 0)
        {
            receivesphereplanewallcontact.resize(receivesphereplanewallcontactNumber);
            MPI_Irecv(receivesphereplanewallcontact.data(), receivesphereplanewallcontactNumber, mpi_contactinformationDataType, rankId,
                      rankId + 5, MPI_COMM_WORLD, &receiveSpherePlanewallContactRequests[rankId][0]);
        }
        MPI_Recv(&receivespherecylinderwallcontactNumber, 1, MPI_INT, rankId, rankId + 6, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (receivespherecylinderwallcontactNumber > 0)
        {
            receivespherecylinderwallcontact.resize(receivespherecylinderwallcontactNumber);
            MPI_Irecv(receivespherecylinderwallcontact.data(), receivespherecylinderwallcontactNumber, mpi_contactinformationDataType, rankId,
                      rankId + 6, MPI_COMM_WORLD, &receiveSphereCylinderwallContactRequests[rankId][0]);
        }
        MPI_Recv(&receivefiberfibercontactNumber, 1, MPI_INT, rankId, rankId + 7, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (receivefiberfibercontactNumber > 0)
        {
            receivefiberfibercontact.resize(receivefiberfibercontactNumber);
            MPI_Irecv(receivefiberfibercontact.data(), receivefiberfibercontactNumber, mpi_contactinformationDataType, rankId,
                      rankId + 7, MPI_COMM_WORLD, &receiveFiberFiberContactRequests[rankId][0]);
        }
        MPI_Recv(&receivefiberplanewallcontactNumber, 1, MPI_INT, rankId, rankId + 8, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (receivefiberplanewallcontactNumber > 0)
        {
            receivefiberplanewallcontact.resize(receivefiberplanewallcontactNumber);
            MPI_Irecv(receivefiberplanewallcontact.data(), receivefiberplanewallcontactNumber, mpi_contactinformationDataType, rankId,
                      rankId + 8, MPI_COMM_WORLD, &receiveFiberPlanewallContactRequests[rankId][0]);
        }
        MPI_Recv(&receivefibercylinderwallcontactNumber, 1, MPI_INT, rankId, rankId + 9, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (receivefibercylinderwallcontactNumber > 0)
        {
            receivefibercylinderwallcontact.resize(receivefibercylinderwallcontactNumber);
            MPI_Irecv(receivefibercylinderwallcontact.data(), receivefibercylinderwallcontactNumber, mpi_contactinformationDataType, rankId,
                      rankId + 9, MPI_COMM_WORLD, &receiveFiberCylinderwallContactRequests[rankId][0]);
        }
    }

    for (int key : spherekeysToRemove)
    {
        DEMproperties->removesphereParticles(key);
    }

    // Wait to make sure you actually send it out
    for (auto &kv : sendSphereParticlesRequests)
    {
        // auto &target_rank = kv.first;
        auto &reqs = kv.second;
        MPI_Wait(&reqs[0], MPI_STATUS_IGNORE);
        if (reqs[1] != MPI_REQUEST_NULL)
        { // Check if the second send was initialized
            MPI_Wait(&reqs[1], MPI_STATUS_IGNORE);
        }
    }
    for (auto &kv : sendFiberSphereParticlesRequests)
    {
        // auto &target_rank = kv.first;
        auto &reqs = kv.second;
        MPI_Wait(&reqs[0], MPI_STATUS_IGNORE);
        if (reqs[1] != MPI_REQUEST_NULL)
        { // Check if the second send was initialized
            MPI_Wait(&reqs[1], MPI_STATUS_IGNORE);
        }
    }
    for (auto &kv : sendFiberBondsRequests)
    {
        // auto &target_rank = kv.first;
        auto &reqs = kv.second;
        MPI_Wait(&reqs[0], MPI_STATUS_IGNORE);
        if (reqs[1] != MPI_REQUEST_NULL)
        { // Check if the second send was initialized
            MPI_Wait(&reqs[1], MPI_STATUS_IGNORE);
        }
    }
    for (auto &kv : sendSphereSphereContactRequests)
    {
        // auto &target_rank = kv.first;
        auto &reqs = kv.second;
        MPI_Wait(&reqs[0], MPI_STATUS_IGNORE);
        if (reqs[1] != MPI_REQUEST_NULL)
        { // Check if the second send was initialized
            MPI_Wait(&reqs[1], MPI_STATUS_IGNORE);
        }
    }
    for (auto &kv : sendSphereFiberContactRequests)
    {
        // auto &target_rank = kv.first;
        auto &reqs = kv.second;
        MPI_Wait(&reqs[0], MPI_STATUS_IGNORE);
        if (reqs[1] != MPI_REQUEST_NULL)
        { // Check if the second send was initialized
            MPI_Wait(&reqs[1], MPI_STATUS_IGNORE);
        }
    }
    for (auto &kv : sendSpherePlanewallContactRequests)
    {
        // auto &target_rank = kv.first;
        auto &reqs = kv.second;
        MPI_Wait(&reqs[0], MPI_STATUS_IGNORE);
        if (reqs[1] != MPI_REQUEST_NULL)
        { // Check if the second send was initialized
            MPI_Wait(&reqs[1], MPI_STATUS_IGNORE);
        }
    }
    for (auto &kv : sendSphereCylinderwallContactRequests)
    {
        // auto &target_rank = kv.first;
        auto &reqs = kv.second;
        MPI_Wait(&reqs[0], MPI_STATUS_IGNORE);
        if (reqs[1] != MPI_REQUEST_NULL)
        { // Check if the second send was initialized
            MPI_Wait(&reqs[1], MPI_STATUS_IGNORE);
        }
    }
    for (auto &kv : sendFiberFiberContactRequests)
    {
        // auto &target_rank = kv.first;
        auto &reqs = kv.second;
        MPI_Wait(&reqs[0], MPI_STATUS_IGNORE);
        if (reqs[1] != MPI_REQUEST_NULL)
        { // Check if the second send was initialized
            MPI_Wait(&reqs[1], MPI_STATUS_IGNORE);
        }
    }
    for (auto &kv : sendFiberPlanewallContactRequests)
    {
        // auto &target_rank = kv.first;
        auto &reqs = kv.second;
        MPI_Wait(&reqs[0], MPI_STATUS_IGNORE);
        if (reqs[1] != MPI_REQUEST_NULL)
        { // Check if the second send was initialized
            MPI_Wait(&reqs[1], MPI_STATUS_IGNORE);
        }
    }
    for (auto &kv : sendFiberCylinderwallContactRequests)
    {
        // auto &target_rank = kv.first;
        auto &reqs = kv.second;
        MPI_Wait(&reqs[0], MPI_STATUS_IGNORE);
        if (reqs[1] != MPI_REQUEST_NULL)
        { // Check if the second send was initialized
            MPI_Wait(&reqs[1], MPI_STATUS_IGNORE);
        }
    }

    // Wait to make sure you actually receive it
    for (auto &kv : receiveSphereParticlesRequests)
    {
        int src_rank = kv.first;
        auto &reqs = kv.second;
        if (reqs[0] != MPI_REQUEST_NULL)
        {
            MPI_Wait(&reqs[0], MPI_STATUS_IGNORE);
        }
        auto &receivesphereparticle = receivesphereparticles[src_rank];

        // Put the nodes to the ghost_nodes
        for (const auto &data : receivesphereparticle)
        {
            Eigen::Vector3d position(data.position[0], data.position[1], data.position[2]);
            Eigen::Vector3d velocity(data.velocity[0], data.velocity[1], data.velocity[2]);
            Eigen::Vector3d omega(data.omega[0], data.omega[1], data.omega[2]);
            Eigen::Vector3d force(data.force[0], data.force[1], data.force[2]);
            Eigen::Vector3d torque(data.torque[0], data.torque[1], data.torque[2]);

            PropertyTypeID type = PropertyTypeID(data.category, data.subType);
            auto particle = std::make_unique<SphereParticle>(data.id, type, data.state, DEMproperties->getParticleManager(), position, velocity, omega, force, torque);
            DEMproperties->addsphereParticles(particle);
        }
    }
    for (auto &kv : receiveFiberSphereParticlesRequests)
    {
        int src_rank = kv.first;
        auto &reqs = kv.second;
        if (reqs[0] != MPI_REQUEST_NULL)
        {
            MPI_Wait(&reqs[0], MPI_STATUS_IGNORE);
        }
        auto &recievefibersphereparticle = receivefibersphereparticles[src_rank];

        for (const auto &data : recievefibersphereparticle)
        {
            Eigen::Vector3d position(data.position[0], data.position[1], data.position[2]);
            Eigen::Vector3d velocity(data.velocity[0], data.velocity[1], data.velocity[2]);
            Eigen::Vector3d omega(data.omega[0], data.omega[1], data.omega[2]);
            Eigen::Vector3d force(data.force[0], data.force[1], data.force[2]);
            Eigen::Vector3d torque(data.torque[0], data.torque[1], data.torque[2]);

            PropertyTypeID type = PropertyTypeID(data.category, data.subType);
            auto particle = std::make_unique<SphereParticle>(
                data.id, type, data.state, DEMproperties->getParticleManager(), position, velocity, omega, force, torque);
            DEMproperties->addfibersphereParticles(particle);
        }
    }
    for (auto &kv : receiveFiberBondsRequests)
    {
        int src_rank = kv.first;
        auto &reqs = kv.second;
        if (reqs[0] != MPI_REQUEST_NULL)
        {
            MPI_Wait(&reqs[0], MPI_STATUS_IGNORE);
        }
        auto &receivefiberbond = receivefiberbonds[src_rank];
        for (const auto &data : receivefiberbond)
        {
            Eigen::Vector3d tangentialforce(data.tangentialforce[0], data.tangentialforce[1], data.tangentialforce[2]);
            Eigen::Vector3d twisttorque(data.twisttorque[0], data.twisttorque[1], data.twisttorque[2]);
            Eigen::Vector3d bendtorque(data.bendtorque[0], data.bendtorque[1], data.bendtorque[2]);

            auto bond = std::make_unique<SphereCylinderBond>(
                data.id,
                PropertyTypeID(data.category, data.subType),
                DEMproperties->getParticleManager(),
                data.fiberid,
                data.node1,
                data.node2,
                data.neighborelement1,
                data.neighborelement2,
                data.energyDissipation);
            DEMproperties->addfiberbonds(bond);
        }
    }
    for (auto &kv : receiveSphereSphereContactRequests)
    {
        int src_rank = kv.first;
        auto &reqs = kv.second;
        if (reqs[0] != MPI_REQUEST_NULL)
        {
            MPI_Wait(&reqs[0], MPI_STATUS_IGNORE);
        }
        auto &receivespherespherecontact = receivespherespherecontacts[src_rank];
        for (const auto &data : receivespherespherecontact)
        {
            ContactInformation info;
            info.particleId1 = data.particleId1;
            info.particleId2 = data.particleId2;
            info.normalForce = Eigen::Vector3d(data.normalForce[0], data.normalForce[1], data.normalForce[2]);
            info.tangentialForce = Eigen::Vector3d(data.tangentialForce[0], data.tangentialForce[1], data.tangentialForce[2]);
            info.tangentialDisplacement = Eigen::Vector3d(data.tangentialDisplacement[0], data.tangentialDisplacement[1], data.tangentialDisplacement[2]);
            info.previousTangentialVelocity = Eigen::Vector3d(data.previousTangentialVelocity[0], data.previousTangentialVelocity[1], data.previousTangentialVelocity[2]);
            DEMproperties->getContactForce()->getSphereSphereContactInformationList()[info.particleId1][info.particleId2] = info;
        }
    }
    for (auto &kv : receiveSphereFiberContactRequests)
    {
        int src_rank = kv.first;
        auto &reqs = kv.second;
        if (reqs[0] != MPI_REQUEST_NULL)
        {
            MPI_Wait(&reqs[0], MPI_STATUS_IGNORE);
        }
        auto &receivespherefibercontact = receivespherefibercontacts[src_rank];
        for (const auto &data : receivespherefibercontact)
        {
            ContactInformation info;
            info.particleId1 = data.particleId1;
            info.particleId2 = data.particleId2;
            info.normalForce = Eigen::Vector3d(data.normalForce[0], data.normalForce[1], data.normalForce[2]);
            info.tangentialForce = Eigen::Vector3d(data.tangentialForce[0], data.tangentialForce[1], data.tangentialForce[2]);
            info.tangentialDisplacement = Eigen::Vector3d(data.tangentialDisplacement[0], data.tangentialDisplacement[1], data.tangentialDisplacement[2]);
            info.previousTangentialVelocity = Eigen::Vector3d(data.previousTangentialVelocity[0], data.previousTangentialVelocity[1], data.previousTangentialVelocity[2]);
            DEMproperties->getContactForce()->getSphereFiberContactInformationList()[info.particleId1][info.particleId2] = info;
        }
    }
    for (auto &kv : receiveSpherePlanewallContactRequests)
    {
        int src_rank = kv.first;
        auto &reqs = kv.second;
        if (reqs[0] != MPI_REQUEST_NULL)
        {
            MPI_Wait(&reqs[0], MPI_STATUS_IGNORE);
        }
        auto &receivesphereplanewallcontact = receivesphereplanewallcontacts[src_rank];
        for (const auto &data : receivesphereplanewallcontact)
        {
            ContactInformation info;
            info.particleId1 = data.particleId1;
            info.particleId2 = data.particleId2;
            info.normalForce = Eigen::Vector3d(data.normalForce[0], data.normalForce[1], data.normalForce[2]);
            info.tangentialForce = Eigen::Vector3d(data.tangentialForce[0], data.tangentialForce[1], data.tangentialForce[2]);
            info.tangentialDisplacement = Eigen::Vector3d(data.tangentialDisplacement[0], data.tangentialDisplacement[1], data.tangentialDisplacement[2]);
            info.previousTangentialVelocity = Eigen::Vector3d(data.previousTangentialVelocity[0], data.previousTangentialVelocity[1], data.previousTangentialVelocity[2]);
            DEMproperties->getContactForce()->getplanewallSphereContactInformationList()[info.particleId1][info.particleId2] = info;
        }
    }
    for (auto &kv : receiveSphereCylinderwallContactRequests)
    {
        int src_rank = kv.first;
        auto &reqs = kv.second;
        if (reqs[0] != MPI_REQUEST_NULL)
        {
            MPI_Wait(&reqs[0], MPI_STATUS_IGNORE);
        }
        auto &receivespherecylinderwallcontact = receivespherecylinderwallcontacts[src_rank];
        for (const auto &data : receivespherecylinderwallcontact)
        {
            ContactInformation info;
            info.particleId1 = data.particleId1;
            info.particleId2 = data.particleId2;
            info.normalForce = Eigen::Vector3d(data.normalForce[0], data.normalForce[1], data.normalForce[2]);
            info.tangentialForce = Eigen::Vector3d(data.tangentialForce[0], data.tangentialForce[1], data.tangentialForce[2]);
            info.tangentialDisplacement = Eigen::Vector3d(data.tangentialDisplacement[0], data.tangentialDisplacement[1], data.tangentialDisplacement[2]);
            info.previousTangentialVelocity = Eigen::Vector3d(data.previousTangentialVelocity[0], data.previousTangentialVelocity[1], data.previousTangentialVelocity[2]);
            DEMproperties->getContactForce()->getcylinderwallSphereContactInformationList()[info.particleId1][info.particleId2] = info;
        }
    }

    for (auto &kv : receiveFiberFiberContactRequests)
    {
        int src_rank = kv.first;
        auto &reqs = kv.second;
        if (reqs[0] != MPI_REQUEST_NULL)
        {
            MPI_Wait(&reqs[0], MPI_STATUS_IGNORE);
        }
        auto &receivefiberfibercontact = receivefiberfibercontacts[src_rank];
        for (const auto &data : receivefiberfibercontact)
        {
            ContactInformation info;
            info.particleId1 = data.particleId1;
            info.particleId2 = data.particleId2;
            info.normalForce = Eigen::Vector3d(data.normalForce[0], data.normalForce[1], data.normalForce[2]);
            info.tangentialForce = Eigen::Vector3d(data.tangentialForce[0], data.tangentialForce[1], data.tangentialForce[2]);
            info.tangentialDisplacement = Eigen::Vector3d(data.tangentialDisplacement[0], data.tangentialDisplacement[1], data.tangentialDisplacement[2]);
            info.previousTangentialVelocity = Eigen::Vector3d(data.previousTangentialVelocity[0], data.previousTangentialVelocity[1], data.previousTangentialVelocity[2]);
            DEMproperties->getContactForce()->getFiberFiberContactInformationList()[info.particleId1][info.particleId2] = info;
        }
    }
    for (auto &kv : receiveFiberPlanewallContactRequests)
    {
        int src_rank = kv.first;
        auto &reqs = kv.second;
        if (reqs[0] != MPI_REQUEST_NULL)
        {
            MPI_Wait(&reqs[0], MPI_STATUS_IGNORE);
        }
        auto &receivefiberplanewallcontact = receivefiberplanewallcontacts[src_rank];
        for (const auto &data : receivefiberplanewallcontact)
        {
            ContactInformation info;
            info.particleId1 = data.particleId1;
            info.particleId2 = data.particleId2;
            info.normalForce = Eigen::Vector3d(data.normalForce[0], data.normalForce[1], data.normalForce[2]);
            info.tangentialForce = Eigen::Vector3d(data.tangentialForce[0], data.tangentialForce[1], data.tangentialForce[2]);
            info.tangentialDisplacement = Eigen::Vector3d(data.tangentialDisplacement[0], data.tangentialDisplacement[1], data.tangentialDisplacement[2]);
            info.previousTangentialVelocity = Eigen::Vector3d(data.previousTangentialVelocity[0], data.previousTangentialVelocity[1], data.previousTangentialVelocity[2]);
            DEMproperties->getContactForce()->getplanewallFiberContactInformationList()[info.particleId1][info.particleId2] = info;
        }
    }
    for (auto &kv : receiveFiberCylinderwallContactRequests)
    {
        int src_rank = kv.first;
        auto &reqs = kv.second;
        if (reqs[0] != MPI_REQUEST_NULL)
        {
            MPI_Wait(&reqs[0], MPI_STATUS_IGNORE);
        }
        auto &receivefibercylinderwallcontact = receivefibercylinderwallcontacts[src_rank];
        for (const auto &data : receivefibercylinderwallcontact)
        {
            ContactInformation info;
            info.particleId1 = data.particleId1;
            info.particleId2 = data.particleId2;
            info.normalForce = Eigen::Vector3d(data.normalForce[0], data.normalForce[1], data.normalForce[2]);
            info.tangentialForce = Eigen::Vector3d(data.tangentialForce[0], data.tangentialForce[1], data.tangentialForce[2]);
            info.tangentialDisplacement = Eigen::Vector3d(data.tangentialDisplacement[0], data.tangentialDisplacement[1], data.tangentialDisplacement[2]);
            info.previousTangentialVelocity = Eigen::Vector3d(data.previousTangentialVelocity[0], data.previousTangentialVelocity[1], data.previousTangentialVelocity[2]);
            DEMproperties->getContactForce()->getcylinderwallFiberContactInformationList()[info.particleId1][info.particleId2] = info;
        }
    }

    for (int key : fiberbondskeysToRemove)
    {
        if (fiberbonds.at(key)->getNeighborelement1() == -1 ||
            fiberbonds.find(fiberbonds.at(key)->getNeighborelement1()) == fiberbonds.end())
        {
            fibersphereskeysToRemove.insert(fiberbonds.at(key)->getNode1());
        }

        if (fiberbonds.at(key)->getNeighborelement2() == -1 ||
            fiberbonds.find(fiberbonds.at(key)->getNeighborelement2()) == fiberbonds.end())
        {
            fibersphereskeysToRemove.insert(fiberbonds.at(key)->getNode2());
        }
    }
    for (int key : fiberbondskeysToRemove)
    {
        DEMproperties->removefiberbonds(key);
    }

    for (int key : fibersphereskeysToRemove)
    {
        DEMproperties->removefibersphereParticles(key);
    }
}

void ParallelDEMModel::runSimulation()
{

    MPI_Barrier(MPI_COMM_WORLD);

    double currentTime = DEMproperties->getCurrentTime();
    double totalTime = DEMproperties->getTotalTime();
    double timeStep = DEMproperties->getTimestep();
    double criticalSpeed = DEMproperties->getCriticalSpeed();
    int taskShowInterval = DEMproperties->getTaskShowInterval();
    int showInterval = DEMproperties->getShowInterval();

    if (simParams.iternum < 0)
    {
        // generate particles
        int generate_number = 0;
        while (!DEMproperties->isGenerateComplete())
        {
            // DEMproperties->applyExternalForces();
            DEMproperties->handleCollisions();
            neighborCommunication();
            DEMproperties->bondforce();
            handleCollisions();
            MPI_Barrier(MPI_COMM_WORLD);
            exchangeSharednodes();

            DEMproperties->motion();
            motion();
            if (generate_number % showInterval == 0)
            {
                collectParticlesFromProcesses();
                if (world_rank == 0)
                {
                    DEMproperties->generateRemainingParticles();
                    vis->Update();
                    simParams.gernerateFiberFlag = DEMproperties->getgernerateFiberFlag();
                    simParams.gernerateSphereFlag = DEMproperties->getGernerateSphereFlag();
                }
                MPI_Bcast(&simParams, 1, mpi_simulationParamsType, 0, MPI_COMM_WORLD);
                if (world_rank != 0)
                {
                    DEMproperties->setGernerateSphereFlag(simParams.gernerateSphereFlag);
                    DEMproperties->setGernerateFiberFlag(simParams.gernerateFiberFlag);
                }

                distributeParticlesToProcesses();
            }
            generate_number++;
        }
        // reach a quasi-static state

        generate_number = 0;
        double averageVelocity = 0;
        while (averageVelocity > criticalSpeed || generate_number < 100000)
        {
            // DEMproperties->applyExternalForces();
            DEMproperties->handleCollisions();
            neighborCommunication();
            DEMproperties->bondforce();
            handleCollisions();
            exchangeSharednodes();
            DEMproperties->motion();
            motion();

            if (generate_number % showInterval == 0)
            {
                collectParticlesFromProcesses();
                if (world_rank == 0)
                {
                    vis->Update();
                    averageVelocity = DEMproperties->getAverageVelocity();
                    std::cout << generate_number << " average velocity is " << averageVelocity << std::endl;
                }
                MPI_Bcast(&averageVelocity, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

                distributeParticlesToProcesses();
            }
            generate_number++;
        }
    }
    // Check if the directory exists
    DEMproperties->initial_task();
    timeStep = DEMproperties->getTimestep();
    std::string file_folder = "DEMProperties";

    if (world_rank == 0)
    {
        std::filesystem::path dirPath(file_folder); // Directory path
        if (!std::filesystem::exists(dirPath))
        {
            // Create the directory if it does not exist
            std::filesystem::create_directories(dirPath);
        }
    }
    while (currentTime < totalTime)
    {
        // DEMproperties->applyExternalForces();
        DEMproperties->handleCollisions();
        neighborCommunication();
        DEMproperties->bondforce();
        handleCollisions();
        MPI_Barrier(MPI_COMM_WORLD);
        exchangeSharednodes();

        std::vector<int> planewallIndex;
        std::vector<int> cylinderwallIndex;
        int planewallnumber = 0;
        int cylinderwallnumber = 0;
        if (world_rank == 0)
        {
            auto &connectedgeometrys = DEMproperties->getConnectedgeometrys();
            for (const auto &connectedgeometry : connectedgeometrys)
            {
                for (auto &geometry : connectedgeometry.second->getGeometrys())
                {
                    if (std::get<0>(geometry) == "PLANEWALL")
                    {
                        int id = std::get<1>(geometry);
                        planewallIndex.push_back(id);
                    }
                    else if (std::get<0>(geometry) == "CYLINDERCONTAINER")
                    {
                        int id = std::get<1>(geometry);
                        cylinderwallIndex.push_back(id);
                    }
                }
            }
            planewallnumber = static_cast<int>(planewallIndex.size());
            cylinderwallnumber = static_cast<int>(cylinderwallIndex.size());
        }
        MPI_Bcast(&planewallnumber, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&cylinderwallnumber, 1, MPI_INT, 0, MPI_COMM_WORLD);
        
        if (world_rank != 0)
        {
            planewallIndex.resize(planewallnumber);
            cylinderwallIndex.resize(cylinderwallnumber);
        }
        if (planewallnumber > 0)
        {
            MPI_Bcast(planewallIndex.data(), planewallnumber, MPI_INT, 0, MPI_COMM_WORLD);
        }
        if (cylinderwallnumber > 0)
        {
            MPI_Bcast(cylinderwallIndex.data(), cylinderwallnumber, MPI_INT, 0, MPI_COMM_WORLD);
        }
        for (int i = 0; i < planewallIndex.size(); ++i)
        {
            collectPlanewallfromProcesses(planewallIndex[i]);
        }
        for (int i = 0; i < cylinderwallIndex.size(); ++i)
        {
            // collectPlanewallfromProcesses(cylinderwallIndex[i]);
        }
        if (world_rank == 0)
        {
            DEMproperties->dotask();
        }
        MPI_Barrier(MPI_COMM_WORLD);
        for (int i = 0; i < planewallIndex.size(); ++i)
        {
            distributePlanewallToProcesses(planewallIndex[i]);
        }
        DEMproperties->motion();
        motion();
        // Update visualization at specified intervals
        if (simParams.iternum % taskShowInterval == 0)
        {
            collectParticlesFromProcesses();
            if (world_rank == 0)
            {
                vis->Update();
                std::string filename = file_folder + "/DEMProperties" + std::to_string(simParams.iternum) + ".dat";
                DEMproperties->saveToFile(filename);
            }
            distributeParticlesToProcesses();
        }
        currentTime += timeStep;
        DEMproperties->setCurrentTime(currentTime);
        simParams.iternum++;
    }
}

void ParallelDEMModel::createMPITypeForSimulationParameters()
{
    const int nitems = 13;                                                                                                                                                            // SimulationParameters中的项目数量
    int blocklengths[nitems] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3};                                                                                                               // 每个项目的长度
    MPI_Datatype types[nitems] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_DOUBLE, MPI_INT, MPI_INT, MPI_C_BOOL, MPI_C_BOOL, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE}; // 每个项目的类型
    MPI_Aint offsets[nitems];                                                                                                                                                         // 每个项目的偏移量

    offsets[0] = offsetof(SimulationParameters, timestep);
    offsets[1] = offsetof(SimulationParameters, taskTimestep);
    offsets[2] = offsetof(SimulationParameters, currentTime);
    offsets[3] = offsetof(SimulationParameters, totalTime);
    offsets[4] = offsetof(SimulationParameters, showInterval);
    offsets[5] = offsetof(SimulationParameters, criticalSpeed);
    offsets[6] = offsetof(SimulationParameters, taskShowInterval);
    offsets[7] = offsetof(SimulationParameters, iternum);
    offsets[8] = offsetof(SimulationParameters, gernerateSphereFlag);
    offsets[9] = offsetof(SimulationParameters, gernerateFiberFlag);

    offsets[10] = offsetof(SimulationParameters, gravity);
    offsets[11] = offsetof(SimulationParameters, globalforces);
    offsets[12] = offsetof(SimulationParameters, simulationdimensions);

    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_simulationParamsType);
    MPI_Type_commit(&mpi_simulationParamsType);
}
void ParallelDEMModel::createMPITypeForSphereParticleData()
{
    const int nitems = 9;
    int blocklengths[nitems] = {1, 1, 1, 1, 3, 3, 3, 3, 3};
    MPI_Datatype types[nitems] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    MPI_Aint offsets[nitems];

    offsets[0] = offsetof(SphereParticleData, id);
    offsets[1] = offsetof(SphereParticleData, category);
    offsets[2] = offsetof(SphereParticleData, subType);
    offsets[3] = offsetof(SphereParticleData, state);
    offsets[4] = offsetof(SphereParticleData, position);
    offsets[5] = offsetof(SphereParticleData, velocity);
    offsets[6] = offsetof(SphereParticleData, omega);
    offsets[7] = offsetof(SphereParticleData, force);
    offsets[8] = offsetof(SphereParticleData, torque);

    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_sphereparticleDataType);
    MPI_Type_commit(&mpi_sphereparticleDataType);
}
void ParallelDEMModel::createMPITypeForFiberBondData()
{
    const int nitems = 12;
    int blocklengths[nitems] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3};
    MPI_Datatype types[nitems] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    MPI_Aint offsets[nitems];

    offsets[0] = offsetof(FiberBondData, id);
    offsets[1] = offsetof(FiberBondData, category);
    offsets[2] = offsetof(FiberBondData, subType);
    offsets[3] = offsetof(FiberBondData, fiberid);
    offsets[4] = offsetof(FiberBondData, node1);
    offsets[5] = offsetof(FiberBondData, node2);
    offsets[6] = offsetof(FiberBondData, neighborelement1);
    offsets[7] = offsetof(FiberBondData, neighborelement2);
    offsets[8] = offsetof(FiberBondData, energyDissipation);
    offsets[9] = offsetof(FiberBondData, tangentialforce);
    offsets[10] = offsetof(FiberBondData, twisttorque);
    offsets[11] = offsetof(FiberBondData, bendtorque);

    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_fiberbondDataType);
    MPI_Type_commit(&mpi_fiberbondDataType);
}
void ParallelDEMModel::createMPITypeForSphereparticleProperties()
{
    const int nitems = 11;                                                                                                                                       // Number of items in SphereParticlePropertiesData
    int blocklengths[nitems] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};                                                                                                // All items are single values
    MPI_Datatype types[nitems] = {MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE}; // Types of items
    MPI_Aint offsets[nitems];                                                                                                                                    // Offsets of items

    // Calculate offsets
    offsets[0] = offsetof(SphereParticlePropertiesData, category);
    offsets[1] = offsetof(SphereParticlePropertiesData, subType);
    offsets[2] = offsetof(SphereParticlePropertiesData, density);
    offsets[3] = offsetof(SphereParticlePropertiesData, rolling_friction_coefficient);
    offsets[4] = offsetof(SphereParticlePropertiesData, slide_friction_coefficient);
    offsets[5] = offsetof(SphereParticlePropertiesData, Young_modulus);
    offsets[6] = offsetof(SphereParticlePropertiesData, restitution);
    offsets[7] = offsetof(SphereParticlePropertiesData, poisson_ratio);
    offsets[8] = offsetof(SphereParticlePropertiesData, momentofinertia);
    offsets[9] = offsetof(SphereParticlePropertiesData, radius);
    offsets[10] = offsetof(SphereParticlePropertiesData, mass);
    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_sphereparticlepropertiesDataType);
    MPI_Type_commit(&mpi_sphereparticlepropertiesDataType);
}
void ParallelDEMModel::createMPITypeForFiberProperties()
{
    const int nitems = 19;                                                                // Number of items in FiberPropertiesData
    int blocklengths[nitems] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}; // All items are single values
    MPI_Datatype types[nitems] = {MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
                                  MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
                                  MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_DOUBLE}; // Types of items
    MPI_Aint offsets[nitems];                                                               // Offsets of items

    // Calculate offsets
    offsets[0] = offsetof(FiberPropertiesData, category);
    offsets[1] = offsetof(FiberPropertiesData, subType);
    offsets[2] = offsetof(FiberPropertiesData, density);
    offsets[3] = offsetof(FiberPropertiesData, rolling_friction_coefficient);
    offsets[4] = offsetof(FiberPropertiesData, slide_friction_coefficient);
    offsets[5] = offsetof(FiberPropertiesData, Young_modulus);
    offsets[6] = offsetof(FiberPropertiesData, restitution);
    offsets[7] = offsetof(FiberPropertiesData, poisson_ratio);
    offsets[8] = offsetof(FiberPropertiesData, aspectratio);
    offsets[9] = offsetof(FiberPropertiesData, nodemomentofinertia);
    offsets[10] = offsetof(FiberPropertiesData, radius);
    offsets[11] = offsetof(FiberPropertiesData, nodemass);
    offsets[12] = offsetof(FiberPropertiesData, elementlength);
    offsets[13] = offsetof(FiberPropertiesData, normalmodulus);
    offsets[14] = offsetof(FiberPropertiesData, shearmodulus);
    offsets[15] = offsetof(FiberPropertiesData, twistmodulus);
    offsets[16] = offsetof(FiberPropertiesData, bendingmodulus);
    offsets[17] = offsetof(FiberPropertiesData, nodenumber);
    offsets[18] = offsetof(FiberPropertiesData, bonddampingcoefficient);

    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_fiberpropertiesDataType);
    MPI_Type_commit(&mpi_fiberpropertiesDataType);
}
void ParallelDEMModel::createMPITypeForPlanewallProperties()
{
    const int nitems = 9;                                                                                                                // Number of items in PlanewallPropertiesData
    int blocklengths[nitems] = {1, 1, 1, 1, 1, 1, 1, 1, 1};                                                                              // All items are single values
    MPI_Datatype types[nitems] = {MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE}; // Types of items
    MPI_Aint offsets[nitems];                                                                                                            // Offsets of items

    // Calculate offsets
    offsets[0] = offsetof(PlanewallPropertiesData, category);
    offsets[1] = offsetof(PlanewallPropertiesData, subType);
    offsets[2] = offsetof(PlanewallPropertiesData, density);
    offsets[3] = offsetof(PlanewallPropertiesData, rolling_friction_coefficient);
    offsets[4] = offsetof(PlanewallPropertiesData, slide_friction_coefficient);
    offsets[5] = offsetof(PlanewallPropertiesData, Young_modulus);
    offsets[6] = offsetof(PlanewallPropertiesData, restitution);
    offsets[7] = offsetof(PlanewallPropertiesData, poisson_ratio);
    offsets[8] = offsetof(PlanewallPropertiesData, thickness);

    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_planewallpropertiesDataType);
    MPI_Type_commit(&mpi_planewallpropertiesDataType);
}
void ParallelDEMModel::createMPITypeForContactInformation()
{
    const int nitems = 6;
    int blocklengths[nitems] = {1, 1, 3, 3, 3, 3};
    MPI_Datatype types[nitems] = {MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE}; // Types of items
    MPI_Aint offsets[nitems];                                                                        // Offsets of items

    // Calculate offsets
    offsets[0] = offsetof(ContactInformationData, particleId1);
    offsets[1] = offsetof(ContactInformationData, particleId2);
    offsets[2] = offsetof(ContactInformationData, tangentialDisplacement);
    offsets[3] = offsetof(ContactInformationData, previousTangentialVelocity);
    offsets[4] = offsetof(ContactInformationData, normalForce);
    offsets[5] = offsetof(ContactInformationData, tangentialForce);

    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_contactinformationDataType);
    MPI_Type_commit(&mpi_contactinformationDataType);
}
void ParallelDEMModel::createMPITypeForPlanewall()
{
    const int nitems = 10;
    int blocklengths[nitems] = {1, 1, 1, 1, 3, 3, 3, 3, 3, 3};
    MPI_Datatype types[nitems] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE}; // Types of items
    MPI_Aint offsets[nitems];                                                                                                                  // Offsets of items

    // Calculate offsets
    offsets[0] = offsetof(PlanewallData, id);
    offsets[1] = offsetof(PlanewallData, category);
    offsets[2] = offsetof(PlanewallData, subType);
    offsets[3] = offsetof(PlanewallData, state);

    offsets[4] = offsetof(PlanewallData, normal);
    offsets[5] = offsetof(PlanewallData, corner1);
    offsets[6] = offsetof(PlanewallData, corner2);
    offsets[7] = offsetof(PlanewallData, corner3);
    offsets[8] = offsetof(PlanewallData, velocity);

    offsets[9] = offsetof(PlanewallData, force);

    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_planewallDataType);
    MPI_Type_commit(&mpi_planewallDataType);
}
void ParallelDEMModel::createMPITypeForCylinderwall()
{
    const int nitems = 9;
    int blocklengths[nitems] = {1, 1, 1, 1, 1, 3, 3, 3, 3};
    MPI_Datatype types[nitems] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE}; // Types of items
    MPI_Aint offsets[nitems];                                                                                                      // Offsets of items

    // Calculate offsets
    offsets[0] = offsetof(CylinderwallData, id);
    offsets[1] = offsetof(CylinderwallData, category);
    offsets[2] = offsetof(CylinderwallData, subType);
    offsets[3] = offsetof(CylinderwallData, state);

    offsets[4] = offsetof(CylinderwallData, radius);

    offsets[5] = offsetof(CylinderwallData, startpoint);
    offsets[6] = offsetof(CylinderwallData, endpoint);

    offsets[7] = offsetof(CylinderwallData, velocity);
    offsets[8] = offsetof(CylinderwallData, force);

    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_cylinderwallDataType);
    MPI_Type_commit(&mpi_cylinderwallDataType);
}

void ParallelDEMModel::createMPITypeForCylinderwallProperties()
{
    const int nitems = 10;                                                                                                                           // Number of items in PlanewallPropertiesData
    int blocklengths[nitems] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};                                                                                       // All items are single values
    MPI_Datatype types[nitems] = {MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_C_BOOL}; // Types of items
    MPI_Aint offsets[nitems];                                                                                                                        // Offsets of items

    // Calculate offsets
    offsets[0] = offsetof(CylinderwallPropertiesData, category);
    offsets[1] = offsetof(CylinderwallPropertiesData, subType);
    offsets[2] = offsetof(CylinderwallPropertiesData, density);
    offsets[3] = offsetof(CylinderwallPropertiesData, rolling_friction_coefficient);
    offsets[4] = offsetof(CylinderwallPropertiesData, slide_friction_coefficient);
    offsets[5] = offsetof(CylinderwallPropertiesData, Young_modulus);
    offsets[6] = offsetof(CylinderwallPropertiesData, restitution);
    offsets[7] = offsetof(CylinderwallPropertiesData, poisson_ratio);
    offsets[8] = offsetof(CylinderwallPropertiesData, thickness);
    offsets[9] = offsetof(CylinderwallPropertiesData, hollow);

    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_cylinderwallpropertiesDataType);
    MPI_Type_commit(&mpi_cylinderwallpropertiesDataType);
}

void ParallelDEMModel::finalizeMPIType()
{
    // 释放MPI数据类型
    MPI_Type_free(&mpi_simulationParamsType);

    MPI_Type_free(&mpi_sphereparticleDataType);
    MPI_Type_free(&mpi_fiberbondDataType);

    MPI_Type_free(&mpi_sphereparticlepropertiesDataType);
    MPI_Type_free(&mpi_fiberpropertiesDataType);
    MPI_Type_free(&mpi_planewallpropertiesDataType);
    MPI_Type_free(&mpi_cylinderwallpropertiesDataType);

    MPI_Type_free(&mpi_contactinformationDataType);

    MPI_Type_free(&mpi_cylinderwallDataType);
    MPI_Type_free(&mpi_planewallDataType);
}
