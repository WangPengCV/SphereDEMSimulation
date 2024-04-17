#include "GridBasedContactDetection.h"
GridBasedContactDetection::GridBasedContactDetection()
{
}
void GridBasedContactDetection::initial(double domain_x, double domain_y, double domain_z, double girdsize)
{
    domainX = domain_x;
    domainY = domain_y;
    domainZ = domain_z;
    numberOfGridX = static_cast<int>(std::ceil(domainX / girdsize));
    numberOfGridY = static_cast<int>(std::ceil(domainY / girdsize));
    numberOfGridZ = static_cast<int>(std::ceil(domainZ / girdsize));

    gridSizeX = domainX / numberOfGridX;
    gridSizeY = domainY / numberOfGridY;
    gridSizeZ = domainZ / numberOfGridZ;

    gridsize = std::min(std::min(gridSizeX, gridSizeY), gridSizeZ);
}

int GridBasedContactDetection::getGridIndex(double x, double y, double z)
{
    int cellX = static_cast<int>(x / gridSizeX);
    if (cellX >= numberOfGridX)
        cellX = numberOfGridX - 1;
    if (cellX < 0)
        cellX = 0;
    int cellY = static_cast<int>(y / gridSizeY);
    if (cellY >= numberOfGridY)
        cellY = numberOfGridY - 1;
    if (cellY < 0)
        cellY = 0;
    int cellZ = static_cast<int>(z / gridSizeZ);
    if (cellZ >= numberOfGridZ)
        cellZ = numberOfGridY - 1;
    if (cellZ < 0)
        cellZ = 0;
    return cellX + cellZ * numberOfGridX + cellY * numberOfGridX * numberOfGridZ;
}

void GridBasedContactDetection::assignParticleToGrid(const std::unordered_map<int, std::unique_ptr<SphereParticle>> &sphereparticles,
                                                     const std::unordered_map<int, std::unique_ptr<SphereCylinderBond>> &fiberbonds,
                                                     const std::unordered_map<int, std::unique_ptr<SphereParticle>> &fibersphereparticles)
{
    gridSaveSphereParticles.clear(); // Clear previous data
    gridSaveFiberBonds.clear();
    for (const auto &particle : sphereparticles)
    {
        int category = particle.second->getType().getCategory();
        auto manager = particle.second->getParticlePropertyManager();

        auto radius = manager->getSphereProperties(particle.second->getType())->getRadius();
        int minCellX = static_cast<int>((particle.second->getPosition().x() - radius) / gridSizeX);
        if (minCellX >= numberOfGridX)
            minCellX = numberOfGridX - 1;
        if (minCellX < 0)
            minCellX = 0;

        int minCellY = static_cast<int>((particle.second->getPosition().y() - radius) / gridSizeY);
        if (minCellY >= numberOfGridY)
            minCellY = numberOfGridY - 1;
        if (minCellY < 0)
            minCellY = 0;

        int minCellZ = static_cast<int>((particle.second->getPosition().z() - radius) / gridSizeZ);
        if (minCellZ >= numberOfGridZ)
            minCellZ = numberOfGridZ - 1;
        if (minCellZ < 0)
            minCellZ = 0;

        int maxCellX = static_cast<int>((particle.second->getPosition().x() + radius) / gridSizeX);
        if (maxCellX >= numberOfGridX)
            maxCellX = numberOfGridX - 1;
        if (maxCellX < 0)
            maxCellX = 0;

        int maxCellY = static_cast<int>((particle.second->getPosition().y() + radius) / gridSizeY);
        if (maxCellY >= numberOfGridY)
            maxCellY = numberOfGridY - 1;
        if (maxCellY < 0)
            maxCellY = 0;

        int maxCellZ = static_cast<int>((particle.second->getPosition().z() + radius) / gridSizeZ);
        if (maxCellZ >= numberOfGridZ)
            maxCellZ = numberOfGridZ - 1;
        if (maxCellZ < 0)
            maxCellZ = 0;
        // Iterate over the range of grid cells and assign the particle's ID
        for (int x = minCellX; x <= maxCellX; x++)
        {
            for (int z = minCellZ; z <= maxCellZ; z++)
            {
                for (int y = minCellY; y <= maxCellY; y++)
                {

                    int gridIndex = x + z * numberOfGridX + y * numberOfGridX * numberOfGridZ;
                    gridSaveSphereParticles[gridIndex].insert(particle.second->getId());
                }
            }
        }
    }
    for (const auto &particle : fiberbonds)
    {
        int category = particle.second->getType().getCategory();
        auto manager = particle.second->getParticlePropertyManager();

        auto radius = manager->getFiberProperties(particle.second->getType())->getRadius();

        int node1Id = particle.second->getNode1();
        int node2Id = particle.second->getNode2();

        Eigen::Vector3d startposition = fibersphereparticles.at(node1Id)->getPosition();
        Eigen::Vector3d endposition = fibersphereparticles.at(node2Id)->getPosition();

        int numSegments = static_cast<int>(std::ceil(manager->getFiberProperties(particle.second->getType())->getElementlength() / (2 * radius)));
        Eigen::Vector3d segmentVector = (endposition - startposition) / numSegments;
        for (int i = 0; i <= numSegments; ++i)
        {
            Eigen::Vector3d center;
            if (i < numSegments)
            {
                center = startposition + i * segmentVector;
            }
            else
            {
                // For the last segment, set the center to the endposition
                center = endposition;
            }

            int minCellX = static_cast<int>((center.x() - 1.42 * radius) / gridSizeX);
            if (minCellX >= numberOfGridX)
                minCellX = numberOfGridX - 1;
            if (minCellX < 0)
                minCellX = 0;

            int minCellY = static_cast<int>((center.y() - 1.42 * radius) / gridSizeY);
            if (minCellY >= numberOfGridY)
                minCellY = numberOfGridY - 1;
            if (minCellY < 0)
                minCellY = 0;

            int minCellZ = static_cast<int>((center.z() - 1.42 * radius) / gridSizeZ);
            if (minCellZ >= numberOfGridZ)
                minCellZ = numberOfGridZ - 1;
            if (minCellZ < 0)
                minCellZ = 0;

            int maxCellX = static_cast<int>((center.x() + 1.42 * radius) / gridSizeX);
            if (maxCellX >= numberOfGridX)
                maxCellX = numberOfGridX - 1;
            if (maxCellX < 0)
                maxCellX = 0;

            int maxCellY = static_cast<int>((center.y() + 1.42 * radius) / gridSizeY);
            if (maxCellY >= numberOfGridY)
                maxCellY = numberOfGridY - 1;
            if (maxCellY < 0)
                maxCellY = 0;

            int maxCellZ = static_cast<int>((center.z() + 1.42 * radius) / gridSizeZ);
            if (maxCellZ >= numberOfGridZ)
                maxCellZ = numberOfGridZ - 1;
            if (maxCellZ < 0)
                maxCellZ = 0;

            // Iterate over the range of grid cells and assign the particle's ID
            for (int x = minCellX; x <= maxCellX; x++)
            {
                for (int z = minCellZ; z <= maxCellZ; z++)
                {
                    for (int y = minCellY; y <= maxCellY; y++)
                    {

                        int gridIndex = x + z * numberOfGridX + y * numberOfGridX * numberOfGridZ;
                        gridSaveFiberBonds[gridIndex].insert(particle.second->getId());
                    }
                }
            }
        }
    }
}

void GridBasedContactDetection::ParticleBroadPhase(const std::unordered_map<int, std::unique_ptr<SphereParticle>> &sphereparticles,
                                                   const std::unordered_map<int, std::unique_ptr<SphereCylinderBond>> &fiberbonds,
                                                   const std::unordered_map<int, std::unique_ptr<SphereParticle>> &fibersphereparticles,
                                                   std::unordered_map<int, std::unordered_set<int>> &SScontactparis,
                                                   std::unordered_map<int, std::unordered_set<int>> &SFcontactparis,
                                                   std::unordered_map<int, std::unordered_set<int>> &FFcontactparis)
{
    assignParticleToGrid(sphereparticles, fiberbonds, fibersphereparticles);
    for (const auto &grid_entry : gridSaveSphereParticles)
    {
        const auto &sphere_in_grid = grid_entry.second;

        for (auto it1 = sphere_in_grid.begin(); it1 != sphere_in_grid.end(); ++it1)
        {

            for (auto it2 = std::next(it1); it2 != sphere_in_grid.end(); ++it2)
            {
                int id1 = *it1;
                int id2 = *it2;
                SScontactparis[std::min(id1, id2)].insert(std::max(id1, id2));
            }
            if (gridSaveFiberBonds.find(grid_entry.first) != gridSaveFiberBonds.end())
            {
                for (auto it3 = gridSaveFiberBonds[grid_entry.first].begin(); it3 != gridSaveFiberBonds[grid_entry.first].end(); ++it3)
                {
                    int sphere_id = *it1;
                    int fiber_id = *it3;
                    SFcontactparis[sphere_id].insert(fiber_id);
                }
            }
        }
    }
    for (const auto &grid_entry : gridSaveFiberBonds)
    {
        const auto &fibers_in_grid = grid_entry.second;

        for (auto it1 = fibers_in_grid.begin(); it1 != fibers_in_grid.end(); ++it1)
        {
            auto it2 = it1;
            ++it2; // Start the second iterator ahead of the first
            for (; it2 != fibers_in_grid.end(); ++it2)
            {
                int id1 = *it1;
                int id2 = *it2;
                if (fiberbonds.at(id1)->getFiberId() != fiberbonds.at(id2)->getFiberId())
                {
                    FFcontactparis[std::min(id1, id2)].insert(std::max(id1, id2));
                }
            }
        }
    }
}

void GridBasedContactDetection::ghostParticleBroadPhase(const std::unordered_map<int, std::unique_ptr<SphereParticle>> &ghostsphereparticles,
                                                        const std::unordered_map<int, std::unique_ptr<SphereCylinderBond>> &ghostfiberbonds,
                                                        const std::unordered_map<int, std::unique_ptr<SphereParticle>> &ghostfibershpereparticles,
                                                        const std::unordered_map<int, std::unique_ptr<SphereCylinderBond>> &fiberbonds,
                                                        std::unordered_map<int, std::unordered_set<int>> &SScontactparis,
                                                        std::unordered_map<int, std::unordered_set<int>> &SFcontactparis,
                                                        std::unordered_map<int, std::unordered_set<int>> &FFcontactparis,
                                                        std::unordered_map<int, std::unordered_set<int>> &neighborbonds)
{
    neighborbonds.clear();
    std::unordered_map<int, std::unordered_set<int>> gridSaveghostSphereParticles;
    std::unordered_map<int, std::unordered_set<int>> gridSaveghostFiberBonds;
    for (const auto &particle : ghostsphereparticles)
    {
        int category = particle.second->getType().getCategory();
        auto manager = particle.second->getParticlePropertyManager();

        auto radius = manager->getSphereProperties(particle.second->getType())->getRadius();
        int minCellX = static_cast<int>((particle.second->getPosition().x() - radius) / gridSizeX);
        if (minCellX >= numberOfGridX)
            minCellX = numberOfGridX - 1;
        if (minCellX < 0)
            minCellX = 0;

        int minCellY = static_cast<int>((particle.second->getPosition().y() - radius) / gridSizeY);
        if (minCellY >= numberOfGridY)
            minCellY = numberOfGridY - 1;
        if (minCellY < 0)
            minCellY = 0;

        int minCellZ = static_cast<int>((particle.second->getPosition().z() - radius) / gridSizeZ);
        if (minCellZ >= numberOfGridZ)
            minCellZ = numberOfGridZ - 1;
        if (minCellZ < 0)
            minCellZ = 0;

        int maxCellX = static_cast<int>((particle.second->getPosition().x() + radius) / gridSizeX);
        if (maxCellX >= numberOfGridX)
            maxCellX = numberOfGridX - 1;
        if (maxCellX < 0)
            maxCellX = 0;

        int maxCellY = static_cast<int>((particle.second->getPosition().y() + radius) / gridSizeY);
        if (maxCellY >= numberOfGridY)
            maxCellY = numberOfGridY - 1;
        if (maxCellY < 0)
            maxCellY = 0;

        int maxCellZ = static_cast<int>((particle.second->getPosition().z() + radius) / gridSizeZ);
        if (maxCellZ >= numberOfGridZ)
            maxCellZ = numberOfGridZ - 1;
        if (maxCellZ < 0)
            maxCellZ = 0;
        // Iterate over the range of grid cells and assign the particle's ID
        for (int x = minCellX; x <= maxCellX; x++)
        {
            for (int z = minCellZ; z <= maxCellZ; z++)
            {
                for (int y = minCellY; y <= maxCellY; y++)
                {

                    int gridIndex = x + z * numberOfGridX + y * numberOfGridX * numberOfGridZ;
                    gridSaveghostSphereParticles[gridIndex].insert(particle.second->getId());
                }
            }
        }
    }
    for (const auto &particle : ghostfiberbonds)
    {
        int category = particle.second->getType().getCategory();
        auto manager = particle.second->getParticlePropertyManager();

        auto radius = manager->getFiberProperties(particle.second->getType())->getRadius();

        int node1Id = particle.second->getNode1();
        int node2Id = particle.second->getNode2();

        Eigen::Vector3d startposition = ghostfibershpereparticles.at(node1Id)->getPosition();
        Eigen::Vector3d endposition = ghostfibershpereparticles.at(node2Id)->getPosition();

        int numSegments = static_cast<int>(std::ceil(manager->getFiberProperties(particle.second->getType())->getElementlength() / (2 * radius)));
        Eigen::Vector3d segmentVector = (endposition - startposition) / numSegments;
        for (int i = 0; i <= numSegments; ++i)
        {
            Eigen::Vector3d center;
            if (i < numSegments)
            {
                center = startposition + i * segmentVector;
            }
            else
            {
                // For the last segment, set the center to the endposition
                center = endposition;
            }

            int minCellX = static_cast<int>((center.x() - 1.42 * radius) / gridSizeX);
            if (minCellX >= numberOfGridX)
                minCellX = numberOfGridX - 1;
            if (minCellX < 0)
                minCellX = 0;

            int minCellY = static_cast<int>((center.y() - 1.42 * radius) / gridSizeY);
            if (minCellY >= numberOfGridY)
                minCellY = numberOfGridY - 1;
            if (minCellY < 0)
                minCellY = 0;

            int minCellZ = static_cast<int>((center.z() - 1.42 * radius) / gridSizeZ);
            if (minCellZ >= numberOfGridZ)
                minCellZ = numberOfGridZ - 1;
            if (minCellZ < 0)
                minCellZ = 0;

            int maxCellX = static_cast<int>((center.x() + 1.42 * radius) / gridSizeX);
            if (maxCellX >= numberOfGridX)
                maxCellX = numberOfGridX - 1;
            if (maxCellX < 0)
                maxCellX = 0;

            int maxCellY = static_cast<int>((center.y() + 1.42 * radius) / gridSizeY);
            if (maxCellY >= numberOfGridY)
                maxCellY = numberOfGridY - 1;
            if (maxCellY < 0)
                maxCellY = 0;

            int maxCellZ = static_cast<int>((center.z() + 1.42 * radius) / gridSizeZ);
            if (maxCellZ >= numberOfGridZ)
                maxCellZ = numberOfGridZ - 1;
            if (maxCellZ < 0)
                maxCellZ = 0;

            // Iterate over the range of grid cells and assign the particle's ID
            for (int x = minCellX; x <= maxCellX; x++)
            {
                for (int z = minCellZ; z <= maxCellZ; z++)
                {
                    for (int y = minCellY; y <= maxCellY; y++)
                    {

                        int gridIndex = x + z * numberOfGridX + y * numberOfGridX * numberOfGridZ;
                        gridSaveghostFiberBonds[gridIndex].insert(particle.second->getId());
                    }
                }
            }
        }
    }

    for (const auto &grid_entry : gridSaveghostSphereParticles)
    {
        int gridId = grid_entry.first;
        const auto &sphere_in_grid = grid_entry.second;

        if (gridSaveSphereParticles.find(gridId) != gridSaveSphereParticles.end())
        {
            for (auto it1 = sphere_in_grid.begin(); it1 != sphere_in_grid.end(); ++it1)
            {
                for (auto it2 = gridSaveSphereParticles[gridId].begin(); it2 != gridSaveSphereParticles[gridId].end(); ++it2)
                {
                    int id1 = *it1;
                    int id2 = *it2;
                    SScontactparis[std::min(id1, id2)].insert(std::max(id1, id2));
                }
            }
        }
        if (gridSaveFiberBonds.find(gridId) != gridSaveFiberBonds.end())
        {
            for (auto it1 = sphere_in_grid.begin(); it1 != sphere_in_grid.end(); ++it1)
            {
                for (auto it3 = gridSaveFiberBonds[gridId].begin(); it3 != gridSaveFiberBonds[gridId].end(); ++it3)
                {
                    int sphere_id = *it1;
                    int fiber_id = *it3;
                    SFcontactparis[sphere_id].insert(fiber_id);
                }
            }
        }
    }

    for (const auto &grid_entry : gridSaveghostFiberBonds)
    {
        int gridId = grid_entry.first;
        const auto &fiber_in_grid = grid_entry.second;

        if (gridSaveSphereParticles.find(gridId) != gridSaveSphereParticles.end())
        {
            for (auto it1 = fiber_in_grid.begin(); it1 != fiber_in_grid.end(); ++it1)
            {
                for (auto it2 = gridSaveSphereParticles[gridId].begin(); it2 != gridSaveSphereParticles[gridId].end(); ++it2)
                {
                    int fiber_id = *it1;
                    int sphere_id = *it2;
                    SFcontactparis[sphere_id].insert(fiber_id);
                }
            }
        }
        if (gridSaveFiberBonds.find(gridId) != gridSaveFiberBonds.end())
        {
            for (auto it1 = fiber_in_grid.begin(); it1 != fiber_in_grid.end(); ++it1)
            {
                for (auto it3 = gridSaveFiberBonds[gridId].begin(); it3 != gridSaveFiberBonds[gridId].end(); ++it3)
                {
                    int id1 = *it1;
                    int id2 = *it3;
                    if (ghostfiberbonds.at(id1)->getFiberId() != fiberbonds.at(id2)->getFiberId())
                    {
                        FFcontactparis[std::min(id1, id2)].insert(std::max(id1, id2));
                    }
                    else
                    {
                        if (ghostfiberbonds.at(id1)->getNeighborelement1() == id2 || ghostfiberbonds.at(id1)->getNeighborelement2() == id2)
                        {
                            neighborbonds[id1].insert(id2);
                        }
                    }
                }
            }
        }
    }
}
void GridBasedContactDetection::planewallBroadPhase(const std::unordered_map<int, std::unique_ptr<PlaneWall>> &planewalls,
                                                    std::unordered_map<int, std::unordered_set<int>> &PScontactpairs,
                                                    std::unordered_map<int, std::unordered_set<int>> &PFcontactpairs)
{

    for (auto &plane : planewalls)
    {
        if (plane.second->getState() == 1 || planeSaveGrid.find(plane.second->getId()) == planeSaveGrid.end())
        {
            plane.second->generateMesh(gridsize);
            std::vector<Eigen::Vector3d> meshvertices = plane.second->getMeshVertices();
            planeSaveGrid[plane.second->getId()].clear();
            std::unordered_set<int> &gridIndices = planeSaveGrid[plane.second->getId()];

            for (auto &vertice : meshvertices)
            {

                gridIndices.insert(getGridIndex(vertice.x(), vertice.y(), vertice.z()));
            }
        }
    }

    for (const auto &grid_entry : planeSaveGrid)
    {
        auto wall_id = grid_entry.first;
        const auto &grid_indices = grid_entry.second;
        auto &pscontactSet = PScontactpairs[wall_id];
        auto &pfcontactSet = PFcontactpairs[wall_id];

        for (int grid_id : grid_indices)
        {
            if (gridSaveSphereParticles.find(grid_id) != gridSaveSphereParticles.end())
            {
                const auto &particlesInGrid = gridSaveSphereParticles[grid_id];

                pscontactSet.insert(particlesInGrid.begin(), particlesInGrid.end());
            }

            if (gridSaveFiberBonds.find(grid_id) != gridSaveFiberBonds.end())
            {
                const auto &particlesInGrid = gridSaveFiberBonds[grid_id];

                pfcontactSet.insert(particlesInGrid.begin(), particlesInGrid.end());
            }
        }
    }
}
void GridBasedContactDetection::RectangularContainerBroadPhase(const std::unique_ptr<RectangularContainer> &rectangularcontainer,
                                                               std::unordered_map<int, std::unordered_set<int>> &PScontactpairs,
                                                               std::unordered_map<int, std::unordered_set<int>> &PFcontactpairs)
{

    for (auto &plane : rectangularcontainer->getPlaneWall())
    {
        // if the plane can move freely or rectangularcontainerSaveGrid have no gird
        if (plane->getState() == 1 || rectangularcontainerSaveGrid.find(plane->getId()) == rectangularcontainerSaveGrid.end())
        {
            plane->generateMesh(gridsize);
            std::vector<Eigen::Vector3d> meshvertices = plane->getMeshVertices();
            rectangularcontainerSaveGrid[plane->getId()].clear();
            std::unordered_set<int> &gridIndices = rectangularcontainerSaveGrid[plane->getId()];

            for (auto &vertice : meshvertices)
            {
                gridIndices.insert(getGridIndex(vertice.x(), vertice.y(), vertice.z()));
            }
        }
    }

    for (const auto &grid_entry : rectangularcontainerSaveGrid)
    {
        auto wall_id = grid_entry.first;
        const auto &grid_indices = grid_entry.second;
        auto &pscontactSet = PScontactpairs[wall_id];
        auto &pfcontactSet = PFcontactpairs[wall_id];

        for (int grid_id : grid_indices)
        {
            if (gridSaveSphereParticles.find(grid_id) != gridSaveSphereParticles.end())
            {
                const auto &particlesInGrid = gridSaveSphereParticles[grid_id];

                pscontactSet.insert(particlesInGrid.begin(), particlesInGrid.end());
            }

            if (gridSaveFiberBonds.find(grid_id) != gridSaveFiberBonds.end())
            {
                const auto &particlesInGrid = gridSaveFiberBonds[grid_id];

                pfcontactSet.insert(particlesInGrid.begin(), particlesInGrid.end());
            }
        }
    }
}
void GridBasedContactDetection::cylinderwallBroadPhase(const std::unordered_map<int, std::unique_ptr<CylinderContainer>> &cylinderwalls,
                                                       std::unordered_map<int, std::unordered_set<int>> &CScontactparis,
                                                       std::unordered_map<int, std::unordered_set<int>> &CFcontactpairs)
{
    for (auto &cylinderwall : cylinderwalls)
    {
        if (cylinderwall.second->getState() == 1 || cylindercontainerSaveGrid.find(cylinderwall.second->getId()) == cylindercontainerSaveGrid.end())
        {
            cylinderwall.second->generateMesh(gridsize);
            const std::vector<Eigen::Vector3d> &meshvertices = cylinderwall.second->getMeshVertices();
            cylindercontainerSaveGrid[cylinderwall.second->getId()].clear();
            std::unordered_set<int> &gridIndices = cylindercontainerSaveGrid[cylinderwall.second->getId()];

            for (auto &vertice : meshvertices)
            {

                gridIndices.insert(getGridIndex(vertice.x(), vertice.y(), vertice.z()));
            }
        }
    }

    for (const auto &grid_entry : cylindercontainerSaveGrid)
    {
        auto wall_id = grid_entry.first;
        const auto &grid_indices = grid_entry.second;
        auto &cscontactSet = CScontactparis[wall_id];
        auto &cfcontactSet = CFcontactpairs[wall_id];

        for (int grid_id : grid_indices)
        {
            if (gridSaveSphereParticles.find(grid_id) != gridSaveSphereParticles.end())
            {
                const auto &particlesInGrid = gridSaveSphereParticles[grid_id];

                cscontactSet.insert(particlesInGrid.begin(), particlesInGrid.end());
            }

            if (gridSaveFiberBonds.find(grid_id) != gridSaveFiberBonds.end())
            {
                const auto &particlesInGrid = gridSaveFiberBonds[grid_id];

                cfcontactSet.insert(particlesInGrid.begin(), particlesInGrid.end());
            }
        }
    }
    // for (const auto &cylinder : cylinderwalls)
    // {
    //     auto wall_id = cylinder.first;

    //     auto &cscontactSet = CScontactparis[wall_id];
    //     auto &cfcontactSet = CFcontactpairs[wall_id];
    //     for(const auto& particlesInGrid :  gridSaveSphereParticles)
    //     {
    //         cscontactSet.insert(particlesInGrid.second.begin(), particlesInGrid.second.end());
    //     }
    // }
}