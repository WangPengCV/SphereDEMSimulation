#include "DEMProperties.h"
DEMProperties::DEMProperties()
{
    manager = std::make_shared<ParticlePropertyManager>();
    contactforce = std::make_shared<ContactForce>();
    gridbasedcontactdetection = std::make_shared<GridBasedContactDetection>();

    globalforces = Eigen::Vector3d::Zero();
    gravity = Eigen::Vector3d::Zero();
    gernerateSphereFlag = true;
    gernerateFiberFlag = true;
    contactEnergydissipation = 0;
    bondEnergydissipation = 0;
}
void DEMProperties::loadFromFile(const std::string &filename)
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        throw std::runtime_error("Unable to open file: " + filename);
    }
    else
    {

        std::string line;
        while (std::getline(file, line))
        {
            line_process(line);
            parseLine(line);
        }
    }
}
bool DEMProperties::parseVector3d(std::istringstream &iss, Eigen::Vector3d &vec)
{
    for (int i = 0; i < 3; ++i)
    {
        if (!(iss >> vec[i]))
        {
            return false;
        }
    }
    return true;
}
void DEMProperties::parseLine(const std::string &line)
{
    // Parse the line and set properties accordingly
    // Example: "GRAVITY 0 -9.81 0" sets the gravity vector
    // Example: "PARTICLE SPHERE 1 0 0 0" defines a sphere particle
    // and so on...
    if (line.empty() || line[0] == '#')
        return; // Skip empty lines and comments

    std::istringstream iss(line);
    std::string entryType, token;
    iss >> entryType;

    if (entryType == "SPHERE_PROPERTIES")
    {
        parseSphereProperties(iss);
    }
    else if (entryType == "CYLINDER_PROPERTIES")
    {
    }
    else if (entryType == "PLANEWALL_PROPERTIES")
    {
        parsePlanewallProperties(iss);
    }
    else if (entryType == "FIBER_PROPERTIES")
    {
        parseFiberProperties(iss);
    }
    else if (entryType == "CYLINDERWALL_PROPERTIES")
    {
        parseCylinderwallProperties(iss);
    }
    else if (entryType == "BOUNDARY")
    {
        std::string boundaryType;
        iss >> boundaryType;

        if (boundaryType == "PLANEWALL")
        {
            parsePlaneWall(iss);
        }
        else if (boundaryType == "CYLINDERCONTAINER")
        {
            parseCylinderWall(iss);
        }
        // Handle other boundary types similarly
    }
    else if (entryType == "CONNECTEDGEOMETRY")
    {
        parseConnectedGeometry(iss);
    }
    else if (entryType == "RANDOM_PARTICLE")
    {
        parseRandomParticle(iss);
    }
    else if (entryType == "PARTICLE")
    {
        parseSpecificParticle(iss);
    }
    else if (entryType == "Fiber")
    {
        parseFiberParticle(iss);
    }
    else if (entryType == "TIMESTEP")
    {
        iss >> timestep;
    }
    else if (entryType == "TASKTIMESTEP")
    {
        iss >> taskTimestep;
    }
    else if (entryType == "TOTALTIME")
    {
        iss >> totalTime;
    }
    else if (entryType == "CURRENTTIME")
    {
        iss >> currentTime;
    }
    else if (entryType == "SHOWINTERVAL")
    {
        iss >> showInterval;
    }
    else if (entryType == "TASKSHOWINTERVAL")
    {
        iss >> taskShowInterval;
    }
    else if (entryType == "CRITICALSPEED")
    {
        iss >> criticalSpeed;
    }
    else if (entryType == "GRAVITY")
    {
        iss >> gravity.x() >> gravity.y() >> gravity.z();
    }
    else if (entryType == "FORCE")
    {
        std::string forceType;
        iss >> forceType;
        if (forceType == "SPECIFIC")
        {
            Eigen::Vector3d force;
            int particle_id;
            iss >> particle_id >> force.x() >> force.y() >> force.z();
            specificforces[particle_id] = force;
        }
        else if (forceType == "GLOBAL")
        {
            iss >> globalforces.x() >> globalforces.y() >> globalforces.z();
        }
    }
    else if (entryType == "DIMENSIONS")
    {
        iss >> simulationdimensions.x() >> simulationdimensions.y() >> simulationdimensions.z();
    }
    else if (entryType == "contactEnergydissipation")
    {
        iss >> contactEnergydissipation;
    }
    else if (entryType == "bondEnergydissipation")
    {
        iss >> bondEnergydissipation;
    }
    else if (entryType == "TASK")
    {
        int id;
        std::string boundarytype, movetype;

        iss >> id >> movetype;
        if (movetype == "SPRINGOSCILLATOR")
        {
            double k, c, p;
            iss >> k >> c >> p;
            std::shared_ptr<DynamicBoundaryMover> dynamicMover = std::make_shared<SpringOscillator>(k, c, p);
            dynamicBoundaryTask.insert({id, dynamicMover});
        }
        else if (movetype == "CONSTANTMOVE")
        {
            Eigen::Vector3d velocity;
            double period;
            int cycles;
            iss >> velocity.x() >> velocity.y() >> velocity.z() >> period;
            std::shared_ptr<DynamicBoundaryMover> dynamicMover = std::make_shared<ConstantMove>(velocity,period);
            dynamicBoundaryTask.insert({id, dynamicMover});
        }
    }
    else if (entryType == "SphereToSphereContact")
    {
        ContactInformation info;
        iss >> info.particleId1 >> info.particleId2;
        iss >> info.normalForce.x() >> info.normalForce.y() >> info.normalForce.z();
        iss >> info.tangentialForce.x() >> info.tangentialForce.y() >> info.tangentialForce.z();
        iss >> info.tangentialDisplacement.x() >> info.tangentialDisplacement.y() >> info.tangentialDisplacement.z();
        iss >> info.previousTangentialVelocity.x() >> info.previousTangentialVelocity.y() >> info.previousTangentialVelocity.z();

        contactforce->getSphereSphereContactInformationList()[info.particleId1][info.particleId2] = info;
    }
    else if (entryType == "SphereToPlanewallContact")
    {
        ContactInformation info;
        iss >> info.particleId1 >> info.particleId2;
        iss >> info.normalForce.x() >> info.normalForce.y() >> info.normalForce.z();
        iss >> info.tangentialForce.x() >> info.tangentialForce.y() >> info.tangentialForce.z();
        iss >> info.tangentialDisplacement.x() >> info.tangentialDisplacement.y() >> info.tangentialDisplacement.z();
        iss >> info.previousTangentialVelocity.x() >> info.previousTangentialVelocity.y() >> info.previousTangentialVelocity.z();

        contactforce->getplanewallSphereContactInformationList()[info.particleId1][info.particleId2] = info;
    }
    else if (entryType == "SphereToFiberContact")
    {
        ContactInformation info;
        iss >> info.particleId1 >> info.particleId2;
        iss >> info.normalForce.x() >> info.normalForce.y() >> info.normalForce.z();
        iss >> info.tangentialForce.x() >> info.tangentialForce.y() >> info.tangentialForce.z();
        iss >> info.tangentialDisplacement.x() >> info.tangentialDisplacement.y() >> info.tangentialDisplacement.z();
        iss >> info.previousTangentialVelocity.x() >> info.previousTangentialVelocity.y() >> info.previousTangentialVelocity.z();
        contactforce->getSphereFiberContactInformationList()[info.particleId1][info.particleId2] = info;
    }
    else if (entryType == "FiberToFiberContact")
    {
        ContactInformation info;
        iss >> info.particleId1 >> info.particleId2;
        iss >> info.normalForce.x() >> info.normalForce.y() >> info.normalForce.z();
        iss >> info.tangentialForce.x() >> info.tangentialForce.y() >> info.tangentialForce.z();
        iss >> info.tangentialDisplacement.x() >> info.tangentialDisplacement.y() >> info.tangentialDisplacement.z();
        iss >> info.previousTangentialVelocity.x() >> info.previousTangentialVelocity.y() >> info.previousTangentialVelocity.z();

        // Assuming a method to store fiber to fiber contact information
        contactforce->getFiberFiberContactInformationList()[info.particleId1][info.particleId2] = info;
    }
    else if (entryType == "FiberToPlanewallContact")
    {
        ContactInformation info;
        iss >> info.particleId1 >> info.particleId2;
        iss >> info.normalForce.x() >> info.normalForce.y() >> info.normalForce.z();
        iss >> info.tangentialForce.x() >> info.tangentialForce.y() >> info.tangentialForce.z();
        iss >> info.tangentialDisplacement.x() >> info.tangentialDisplacement.y() >> info.tangentialDisplacement.z();
        iss >> info.previousTangentialVelocity.x() >> info.previousTangentialVelocity.y() >> info.previousTangentialVelocity.z();

        // Store the contact information
        contactforce->getplanewallFiberContactInformationList()[info.particleId1][info.particleId2] = info;
    }
    else if (entryType == "FiberToCylinderwallContact")
    {
        ContactInformation info;
        iss >> info.particleId1 >> info.particleId2;
        iss >> info.normalForce.x() >> info.normalForce.y() >> info.normalForce.z();
        iss >> info.tangentialForce.x() >> info.tangentialForce.y() >> info.tangentialForce.z();
        iss >> info.tangentialDisplacement.x() >> info.tangentialDisplacement.y() >> info.tangentialDisplacement.z();
        iss >> info.previousTangentialVelocity.x() >> info.previousTangentialVelocity.y() >> info.previousTangentialVelocity.z();

        // Store the contact information
        contactforce->getcylinderwallFiberContactInformationList()[info.particleId1][info.particleId2] = info;
    }
    else if (entryType == "SphereToCylinderwallContact")
    {
        ContactInformation info;
        iss >> info.particleId1 >> info.particleId2;
        iss >> info.normalForce.x() >> info.normalForce.y() >> info.normalForce.z();
        iss >> info.tangentialForce.x() >> info.tangentialForce.y() >> info.tangentialForce.z();
        iss >> info.tangentialDisplacement.x() >> info.tangentialDisplacement.y() >> info.tangentialDisplacement.z();
        iss >> info.previousTangentialVelocity.x() >> info.previousTangentialVelocity.y() >> info.previousTangentialVelocity.z();

        // Store the contact information
        contactforce->getcylinderwallSphereContactInformationList()[info.particleId1][info.particleId2] = info;
    }
}
void DEMProperties::line_process(std::string &line)
{
    for (char &c : line)
    {
        if (c == '\t' || c == ',' || c == ';' || c == '\r' || c == '\n')
            c = ' ';
    }
    line.erase(0, line.find_first_not_of(" "));
    line.erase(line.find_last_not_of(" ") + 1);
}
void DEMProperties::parseSphereProperties(std::istringstream &iss)
{
    int category_id, subtype;
    double density, radius, rollingFriction, slidingFriction, youngModulus, restitution, poissonRatio;

    try
    {
        iss >> category_id >> subtype >> density >> radius >> rollingFriction >> slidingFriction >> youngModulus >> restitution >> poissonRatio;

        // Use the parsed values...
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error parsing SPHERE_PROPERTIES: " << e.what() << std::endl;
        return;
    }
    // Calculations based on the properties
    double mass = 4.0 / 3.0 * PI * std::pow(radius, 3) * density;
    double moment_of_inertia = 2.0 / 5.0 * mass * std::pow(radius, 2);

    // Create SphereProperties object and store or process as needed
    auto sphereProperties = std::make_shared<SphereProperties>(density, mass, radius, rollingFriction, slidingFriction,
                                                               youngModulus, restitution, poissonRatio, moment_of_inertia);
    // Assuming you have a mechanism to add these properties to your manager
    manager->addSphereProperties(PropertyTypeID(category_id, subtype), sphereProperties);
}
void DEMProperties::parsePlanewallProperties(std::istringstream &iss)
{
    int category_id, subtype;
    double density, thickness, rollingFriction, slidingFriction, youngModulus, restitution, poissonRatio;

    try
    {
        iss >> category_id >> subtype >> density >> thickness >> rollingFriction >> slidingFriction >> youngModulus >> restitution >> poissonRatio;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error parsing Planewall_Propertie: " << e.what() << std::endl;
        return;
    }
    auto properties = std::make_shared<PlanewallProperties>(density, thickness, rollingFriction, slidingFriction,
                                                            youngModulus, restitution, poissonRatio);
    manager->addPlanewallProperties(PropertyTypeID(category_id, subtype), properties);
}
void DEMProperties::parseFiberProperties(std::istringstream &iss)
{
    int category_id, subtype, nodenumber;
    double density, radius, rollingFriction, slidingFriction, youngModulus, restitution,
        poissonRatio, normalmodulus, shearmodulus, twistmodulus, bendingmodulus, aspectratio, bonddampingcoefficient;

    try
    {
        iss >> category_id >> subtype >> density >> radius >> rollingFriction >> slidingFriction >> youngModulus >> restitution >> poissonRatio >> normalmodulus >> shearmodulus >> twistmodulus >> bendingmodulus >> nodenumber >> aspectratio >> bonddampingcoefficient;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error parsing Planewall_Propertie: " << e.what() << std::endl;
        return;
    }
    double elementlength = (aspectratio - 1) * 2 * radius / (nodenumber - 1);
    auto properties = std::make_shared<FiberProperties>(density, rollingFriction, slidingFriction,
                                                        youngModulus, restitution, poissonRatio, radius, elementlength, normalmodulus,
                                                        shearmodulus, twistmodulus, bendingmodulus, nodenumber, aspectratio, bonddampingcoefficient);
    manager->addFiberProperties(PropertyTypeID(category_id, subtype), properties);
}
void DEMProperties::parsePlaneWall(std::istringstream &iss)
{
    int id, category_id, subtype, state;
    Eigen::Vector3d normal, corner1, corner2, corner3, velocity;
    try
    {
        iss >> id >> category_id >> subtype >> state;
        if (!parseVector3d(iss, normal) || !parseVector3d(iss, corner1) || !parseVector3d(iss, corner2) || !parseVector3d(iss, corner3) || !parseVector3d(iss, velocity))
        {
            throw std::runtime_error("Error parsing PLANEWALL boundary vectors.");
        }
        auto pw = std::make_unique<PlaneWall>(id, PropertyTypeID(category_id, subtype), state, normal, corner1, corner2, corner3, velocity);
        planewalls.insert({pw->getId(), std::move(pw)});
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error parsing PlaneWall: " << e.what() << std::endl;
        return;
    }
}
void DEMProperties::parseCylinderwallProperties(std::istringstream &iss)
{
    int category_id, subtype;
    double density, thickness, rollingFriction, slidingFriction, youngModulus, restitution, poissonRatio;
    bool hollow;
    try
    {
        iss >> category_id >> subtype >> density >> thickness >> rollingFriction >> slidingFriction >> youngModulus >> restitution >> poissonRatio >> hollow;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error parsing Cylinderwall_Propertie: " << e.what() << std::endl;
        return;
    }
    auto properties = std::make_shared<CylinderwallProperties>(density, thickness, rollingFriction, slidingFriction,
                                                               youngModulus, restitution, poissonRatio, hollow);
    manager->addCylinerwallProperties(PropertyTypeID(category_id, subtype), properties);
}
void DEMProperties::parseCylinderWall(std::istringstream &iss)
{
    int id, category_id, subtype, state;
    double radius;
    Eigen::Vector3d startpoint, endpoint, velocity;
    try
    {
        iss >> id >> category_id >> subtype >> state >> radius;
        if (!parseVector3d(iss, startpoint) || !parseVector3d(iss, endpoint) || !parseVector3d(iss, velocity))
        {
            throw std::runtime_error("Error parsing CYLINDERWALL boundary vectors.");
        }
        auto cw = std::make_unique<CylinderContainer>(id, PropertyTypeID(category_id, subtype), state, radius, startpoint, endpoint, velocity);
        cylinderwalls.insert({cw->getId(), std::move(cw)});
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error parsing CylinderContainer: " << e.what() << std::endl;
        return;
    }
}
void DEMProperties::parseConnectedGeometry(std::istringstream &iss)
{
    int id, numWalls;
    std::string geometryType, dof;
    int geoId;
    int centerIndex = 0;
    Eigen::Vector3d center = Eigen::Vector3d::Zero();
    Eigen::Vector3d velocity;

    std::vector<std::tuple<std::string, int>> geometrys;
    try
    {
        iss >> id >> numWalls;
        for (int i = 0; i < numWalls; ++i)
        {
            iss >> geometryType >> geoId;
            geometrys.push_back(std::make_tuple(geometryType, geoId));

            if (geometryType == "SPHERE")
            {
                center += sphereparticles[geoId]->getPosition();
                centerIndex++;
                velocity = sphereparticles[geoId]->getVelocity();
                // mass += manager->getSphereProperties(sphereparticles[geoId]->getType())->getMass();
            }
            else if (geometryType == "Fiber")
            {
                for (int j = 0; j < fibers[geoId]->getSphereCylinderBond().size(); ++j)
                {
                    if (j == 0)
                    {
                        int node1Id = fiberbonds[fibers[geoId]->getSphereCylinderBond()[j]]->getNode1();
                        center += fibershpereparticles[node1Id]->getPosition();
                        centerIndex++;
                        velocity = fibershpereparticles[node1Id]->getVelocity();
                    }
                    int node2Id = fiberbonds[fibers[geoId]->getSphereCylinderBond()[j]]->getNode2();
                    center += fibershpereparticles[node2Id]->getPosition();
                    centerIndex++;
                }
                // mass += manager->getFiberProperties(fibers[geoId]->getType())->getNodenumber() * manager->getFiberProperties(fibers[geoId]->getType())->getNodemass();
            }
            else if (geometryType == "CYLINDERCONTAINER")
            {
                center += cylinderwalls[geoId]->getStartpoint();
                centerIndex++;
                center += cylinderwalls[geoId]->getEndpoint();
                centerIndex++;
                velocity = cylinderwalls[geoId]->getVelocity();
                // mass += cylinderwalls[geoId]->getMass();
            }
            else if (geometryType == "PLANEWALL")
            {
                center += planewalls[geoId]->getCorner1();
                centerIndex++;
                center += planewalls[geoId]->getCorner2();
                centerIndex++;
                center += planewalls[geoId]->getCorner3();
                centerIndex++;
                velocity = planewalls[geoId]->getVelocity();
                // mass += planewalls[geoId]->getMass();
            }
        }
        center /= centerIndex;
        auto cg = std::make_unique<ConnectedGeometry>(id, center, velocity);
        cg->setGeometry(geometrys);
        // cg->setMass(mass);
        iss >> dof;
        if (dof == "X")
        {
            cg->setDegreeOfFreedom(ConnectedGeometry::DegreeOfFreedom::Y, false); // 开启 Y 方向的自由度
            cg->setDegreeOfFreedom(ConnectedGeometry::DegreeOfFreedom::Z, false); // 关闭 Z 方向的自由度
        }
        else if (dof == "Y")
        {
            cg->setDegreeOfFreedom(ConnectedGeometry::DegreeOfFreedom::X, false); // 开启 Y 方向的自由度
            cg->setDegreeOfFreedom(ConnectedGeometry::DegreeOfFreedom::Z, false); // 关闭 Z 方向的自由度
        }
        else if (dof == "Z")
        {
            cg->setDegreeOfFreedom(ConnectedGeometry::DegreeOfFreedom::X, false); // 开启 Y 方向的自由度
            cg->setDegreeOfFreedom(ConnectedGeometry::DegreeOfFreedom::Y, false); // 关闭 Z 方向的自由度
        }

        connectedgeometrys.insert({id, std::move(cg)});
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error parsing connected geometrys: " << e.what() << std::endl;
        return;
    }
}
void DEMProperties::setCurrentTime(double currentTime)
{
    this->currentTime = currentTime;
}
void DEMProperties::setTotalTime(double totalTime)
{
    this->totalTime = totalTime;
}
void DEMProperties::setTimestep(double timestep)
{
    this->timestep = timestep;
}
void DEMProperties::setTaskTimestep(double taskTimestep)
{
    this->taskTimestep = taskTimestep;
}
void DEMProperties::setShowInterval(int showInterval)
{
    this->showInterval = showInterval;
}
void DEMProperties::setTaskShowInterval(double taskShowInterval)
{
    this->taskShowInterval = taskShowInterval;
}
void DEMProperties::setCriticalSpeed(double criticalSpeed)
{
    this->criticalSpeed = criticalSpeed;
}

void DEMProperties::setGravity(const Eigen::Vector3d &gravity)
{
    this->gravity = gravity;
}

void DEMProperties::setGlobalForces(const Eigen::Vector3d &globalforces)
{
    this->globalforces = globalforces;
}
void DEMProperties::setSimulationDimensions(const Eigen::Vector3d &simulationdimensions)
{
    this->simulationdimensions = simulationdimensions;
}

void DEMProperties::setsphereParticles(std::unordered_map<int, std::unique_ptr<SphereParticle>> &sphereParticles)
{
    this->sphereparticles.clear();
    for (auto &pair : sphereParticles)
    {
        this->sphereparticles.emplace(pair.first, std::move(pair.second));
    }
}
void DEMProperties::removesphereParticles(int id)
{
    sphereparticles.erase(id);
}
void DEMProperties::addsphereParticles(std::unique_ptr<SphereParticle> &sp)
{
    sphereparticles.insert({sp->getId(), std::move(sp)});
}

void DEMProperties::setfibersphereParticles(std::unordered_map<int, std::unique_ptr<SphereParticle>> &fiberSphereParticles)
{
    this->fibershpereparticles.clear();
    for (auto &pair : fiberSphereParticles)
    {
        this->fibershpereparticles.emplace(pair.first, std::move(pair.second));
    }
}
void DEMProperties::removefibersphereParticles(int id)
{
    fibershpereparticles.erase(id);
}

void DEMProperties::addfibersphereParticles(std::unique_ptr<SphereParticle> &sp)
{
    fibershpereparticles.insert({sp->getId(), std::move(sp)});
}
void DEMProperties::setfiberbonds(std::unordered_map<int, std::unique_ptr<SphereCylinderBond>> &fiberBonds)
{
    this->fiberbonds.clear();

    for (auto &pair : fiberBonds)
    {
        this->fiberbonds.emplace(pair.first, std::move(pair.second));
    }
}
void DEMProperties::removefiberbonds(int id)
{
    fiberbonds.erase(id);
}
void DEMProperties::addfiberbonds(std::unique_ptr<SphereCylinderBond> &sc)
{
    fiberbonds.insert({sc->getId(), std::move(sc)});
}

void DEMProperties::setplanewalls(std::unordered_map<int, std::unique_ptr<PlaneWall>> &planeWalls)
{
    this->planewalls.clear();
    for (auto &pair : planeWalls)
    {
        this->planewalls.emplace(pair.first, std::move(pair.second));
    }
}
void DEMProperties::setplanewall(std::unique_ptr<PlaneWall> &planeWall)
{

    this->planewalls[planeWall->getId()] = std::move(planeWall);
}

void DEMProperties::setcylinderwalls(std::unordered_map<int, std::unique_ptr<CylinderContainer>> &cylinderWalls)
{
    this->cylinderwalls.clear();
    for (auto &pair : cylinderWalls)
    {
        this->cylinderwalls.emplace(pair.first, std::move(pair.second));
    }
}
void DEMProperties::setcylinderwall(std::unique_ptr<CylinderContainer> &cylinderWall)
{
    this->cylinderwalls[cylinderWall->getId()] = std::move(cylinderWall);
}
void DEMProperties::setGernerateFiberFlag(bool gernerateFiberFlag)
{
    this->gernerateFiberFlag = gernerateFiberFlag;
}

void DEMProperties::setGernerateSphereFlag(bool gernerateSphereFlag)
{
    this->gernerateSphereFlag = gernerateSphereFlag;
}

void DEMProperties::saveToFile(const std::string &filename) const
{

    std::ofstream file(filename);
    if (!file.is_open())
    {
        throw std::runtime_error("Unable to open file for saving: " + filename);
    }
    file.precision(std::numeric_limits<double>::digits10 + 1);

    file << "DIMENSIONS,  " << simulationdimensions.x() << ", " << simulationdimensions.y() << ", " << simulationdimensions.z() << "\n";
    file << "TIMESTEP, " << timestep << "\n";
    file << "TASKTIMESTEP, " << taskTimestep << "\n";
    file << "CURRENTTIME, " << currentTime << "\n";
    file << "TOTALTIME, " << totalTime << "\n";
    file << "SHOWINTERVAL, " << showInterval << "\n";
    file << "TASKSHOWINTERVAL, " << taskShowInterval << "\n";
    file << "CRITICALSPEED, " << criticalSpeed << "\n";

    file << "GRAVITY, " << gravity.x() << ", " << gravity.y() << ", " << gravity.z() << "\n";

    for (const auto &pair : specificforces)
    {
        const auto &id = pair.first;
        const auto &force = pair.second;

        file << "FORCE, "
             << "SPECIFIC, " << id << "," << force.x() << ", " << force.y() << ", " << force.z() << "\n";
    }

    file << "FORCE, "
         << "GLOBAL, " << globalforces.x() << ", " << globalforces.y() << ", " << globalforces.z() << "\n";
    file << "contactEnergydissipation, " << contactEnergydissipation << "\n";
    file << "bondEnergydissipation, " << bondEnergydissipation << "\n";
    for (const auto &pair : manager->getParticleProperties())
    {
        const auto &id = pair.first;
        const auto &prop = pair.second;

        if (auto sphereProerties = std::dynamic_pointer_cast<SphereProperties>(prop))
        {
            file << "SPHERE_PROPERTIES, " << id.getCategory() << ", " << id.getSubType() << ", " << prop->save_tostring() << "\n";
        }
        else if (auto planewallProerties = std::dynamic_pointer_cast<PlanewallProperties>(prop))
        {
            file << "PLANEWALL_PROPERTIES, " << id.getCategory() << ", " << id.getSubType() << ", " << prop->save_tostring() << "\n";
        }
        else if (auto planewallProerties = std::dynamic_pointer_cast<FiberProperties>(prop))
        {
            file << "FIBER_PROPERTIES, " << id.getCategory() << ", " << id.getSubType() << ", " << prop->save_tostring() << "\n";
        }
        else if (auto cylinderwallProerties = std::dynamic_pointer_cast<CylinderwallProperties>(prop))
        {
            file << "CYLINDERWALL_PROPERTIES, " << id.getCategory() << ", " << id.getSubType() << ", " << prop->save_tostring() << "\n";
        }
        // Add other property types as needed
    }

    // Write boundary conditions
    for (const auto &boundary : planewalls)
    {
        file << "BOUNDARY, " << boundary.second->save_tostring() << "\n";
    }

    for (const auto &boundary : cylinderwalls)
    {
        file << "BOUNDARY, " << boundary.second->save_tostring() << "\n";
    }

    for (const auto &geometry : connectedgeometrys)
    {
        file << "CONNECTEDGEOMETRY, " << geometry.second->save_tostring() << "\n";
    }
    for (auto &subtask : dynamicBoundaryTask)
    {
        file << "TASK, ";
        file << subtask.first << ", ";
        if (auto &springoscillator = std::dynamic_pointer_cast<SpringOscillator>(subtask.second))
        {
            file << "SPRINGOSCILLATOR, ";
            file << springoscillator->save_tostring() << std::endl;
        }
        else if (auto &constantmove = std::dynamic_pointer_cast<ConstantMove>(subtask.second))
        {
            file << "CONSTANTMOVE, ";
            file << constantmove->save_tostring() << std::endl;
        }
    }

    for (const auto &sphere : sphereparticles)
    {
        file << sphere.second->save_tostring() << std::endl;
    }

    for (const auto &fiber : fibers)
    {
        file << fiber->save_tostring() << std::endl;
    }
    for (const auto &bond : fiberbonds)
    {
        file << bond.second->save_tostring() << std::endl;
    }
    for (const auto &sphere : fibershpereparticles)
    {
        file << sphere.second->save_fibertostring() << std::endl;
    }

    // Saving sphere-to-sphere contact information
    for (const auto &pair1 : contactforce->getSphereSphereContactInformationList())
    {
        for (const auto &pair2 : pair1.second)
        {
            const ContactInformation &info = pair2.second;
            file << "SphereToSphereContact, "
                 << info.particleId1 << ", "
                 << info.particleId2 << ", "
                 << info.normalForce.x() << ", " << info.normalForce.y() << ", " << info.normalForce.z() << ", "
                 << info.tangentialForce.x() << ", " << info.tangentialForce.y() << ", " << info.tangentialForce.z() << ", "
                 << info.tangentialDisplacement.x() << ", " << info.tangentialDisplacement.y() << ", " << info.tangentialDisplacement.z() << ", "
                 << info.previousTangentialVelocity.x() << ", " << info.previousTangentialVelocity.y() << ", " << info.previousTangentialVelocity.z() << "\n";
        }
    }
    // Saving sphere-to-planewall contact information
    for (const auto &pair1 : contactforce->getplanewallSphereContactInformationList())
    {
        for (const auto &pair2 : pair1.second)
        {
            const ContactInformation &info = pair2.second;
            file << "SphereToPlanewallContact, "
                 << info.particleId1 << ", "
                 << info.particleId2 << ", "
                 << info.normalForce.x() << ", " << info.normalForce.y() << ", " << info.normalForce.z() << ", "
                 << info.tangentialForce.x() << ", " << info.tangentialForce.y() << ", " << info.tangentialForce.z() << ", "
                 << info.tangentialDisplacement.x() << ", " << info.tangentialDisplacement.y() << ", " << info.tangentialDisplacement.z() << ", "
                 << info.previousTangentialVelocity.x() << ", " << info.previousTangentialVelocity.y() << ", " << info.previousTangentialVelocity.z() << "\n";
        }
    }
    // Saving sphere-to-fiber contact information
    for (const auto &pair1 : contactforce->getSphereFiberContactInformationList())
    {
        for (const auto &pair2 : pair1.second)
        {
            const ContactInformation &info = pair2.second;
            file << "SphereToFiberContact, "
                 << info.particleId1 << ", "
                 << info.particleId2 << ", "
                 << info.normalForce.x() << ", " << info.normalForce.y() << ", " << info.normalForce.z() << ", "
                 << info.tangentialForce.x() << ", " << info.tangentialForce.y() << ", " << info.tangentialForce.z() << ", "
                 << info.tangentialDisplacement.x() << ", " << info.tangentialDisplacement.y() << ", " << info.tangentialDisplacement.z() << ", "
                 << info.previousTangentialVelocity.x() << ", " << info.previousTangentialVelocity.y() << ", " << info.previousTangentialVelocity.z() << "\n";
        }
    }

    // Saving fiber-to-fiber contact information
    for (const auto &pair1 : contactforce->getFiberFiberContactInformationList())
    {
        for (const auto &pair2 : pair1.second)
        {
            const ContactInformation &info = pair2.second;
            file << "FiberToFiberContact, "
                 << info.particleId1 << ", "
                 << info.particleId2 << ", "
                 << info.normalForce.x() << ", " << info.normalForce.y() << ", " << info.normalForce.z() << ", "
                 << info.tangentialForce.x() << ", " << info.tangentialForce.y() << ", " << info.tangentialForce.z() << ", "
                 << info.tangentialDisplacement.x() << ", " << info.tangentialDisplacement.y() << ", " << info.tangentialDisplacement.z() << ", "
                 << info.previousTangentialVelocity.x() << ", " << info.previousTangentialVelocity.y() << ", " << info.previousTangentialVelocity.z() << "\n";
        }
    }

    // Saving planewall-to-fiber contact information
    for (const auto &pair1 : contactforce->getplanewallFiberContactInformationList())
    {
        for (const auto &pair2 : pair1.second)
        {
            const ContactInformation &info = pair2.second;
            file << "FiberToPlanewallContact, "
                 << info.particleId1 << ", "
                 << info.particleId2 << ", "
                 << info.normalForce.x() << ", " << info.normalForce.y() << ", " << info.normalForce.z() << ", "
                 << info.tangentialForce.x() << ", " << info.tangentialForce.y() << ", " << info.tangentialForce.z() << ", "
                 << info.tangentialDisplacement.x() << ", " << info.tangentialDisplacement.y() << ", " << info.tangentialDisplacement.z() << ", "
                 << info.previousTangentialVelocity.x() << ", " << info.previousTangentialVelocity.y() << ", " << info.previousTangentialVelocity.z() << "\n";
        }
    }
    // Saving cylinderwall-to-sphere contact information
    for (const auto &pair1 : contactforce->getcylinderwallSphereContactInformationList())
    {
        for (const auto &pair2 : pair1.second)
        {
            const ContactInformation &info = pair2.second;
            file << "SphereToCylinderwallContact, "
                 << info.particleId1 << ", "
                 << info.particleId2 << ", "
                 << info.normalForce.x() << ", " << info.normalForce.y() << ", " << info.normalForce.z() << ", "
                 << info.tangentialForce.x() << ", " << info.tangentialForce.y() << ", " << info.tangentialForce.z() << ", "
                 << info.tangentialDisplacement.x() << ", " << info.tangentialDisplacement.y() << ", " << info.tangentialDisplacement.z() << ", "
                 << info.previousTangentialVelocity.x() << ", " << info.previousTangentialVelocity.y() << ", " << info.previousTangentialVelocity.z() << "\n";
        }
    }
    // Saving cylinderwall-to-fiber contact information
    for (const auto &pair1 : contactforce->getcylinderwallFiberContactInformationList())
    {
        for (const auto &pair2 : pair1.second)
        {
            const ContactInformation &info = pair2.second;
            file << "FiberToCylinderwallContact, "
                 << info.particleId1 << ", "
                 << info.particleId2 << ", "
                 << info.normalForce.x() << ", " << info.normalForce.y() << ", " << info.normalForce.z() << ", "
                 << info.tangentialForce.x() << ", " << info.tangentialForce.y() << ", " << info.tangentialForce.z() << ", "
                 << info.tangentialDisplacement.x() << ", " << info.tangentialDisplacement.y() << ", " << info.tangentialDisplacement.z() << ", "
                 << info.previousTangentialVelocity.x() << ", " << info.previousTangentialVelocity.y() << ", " << info.previousTangentialVelocity.z() << "\n";
        }
    }

    file.close();
}
void DEMProperties::parseRandomParticle(std::istringstream &iss)
{

    std::string type;
    iss >> type;
    if (type == "SPHERE")
    {
        int categoryId, subTypeId, state, count;
        double xmin, xmax, ymin, ymax, zmin, zmax;
        iss >> categoryId >> subTypeId >> state >> count >> xmin >> xmax >> ymin >> ymax >> zmin >> zmax;

        int id_index = sphereparticles.size();

        gernerateSphereFlag = true;
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> disX(xmin, xmax);
        std::uniform_real_distribution<> disY(ymin, ymax);
        std::uniform_real_distribution<> disZ(zmin, zmax);
        int failedSpheres = 0;
        double sphereRadius = manager->getSphereProperties(PropertyTypeID(categoryId, subTypeId))->getRadius();

        // Generate 'count' random particles within the specified ranges
        for (int i = 0; i < count; ++i)
        {
            double x, y, z;
            bool validPosition = false;
            int count_number = 0;
            do
            {
                count_number++;
                validPosition = true;
                x = disX(gen);
                y = disY(gen);
                z = disZ(gen);
                Eigen::Vector3d center = Eigen::Vector3d(x, y, z);
                for (const auto &plane : planewalls)
                {
                    // Retrieve the plane wall's point and normal
                    Eigen::Vector3d planePoint = plane.second->getCorner1();

                    Eigen::Vector3d planeNormal = plane.second->getNormal();

                    // Compute the vector from a point on the plane to the sphere's center
                    Eigen::Vector3d vecToSphereCenter = center - planePoint;

                    // Compute the distance from the sphere's center to the plane
                    double distanceToPlane = vecToSphereCenter.dot(planeNormal);

                    // Calculate the overlap (penetration depth)
                    double overlap_pw = sphereRadius - std::abs(distanceToPlane);
                    if (overlap_pw > 0)
                    {
                        validPosition = false;
                        break;
                    }
                }
                if (validPosition)
                {
                    for (const auto &cylinder : cylinderwalls)
                    {
                        Eigen::Vector3d cylinderStart = cylinder.second->getStartpoint();
                        Eigen::Vector3d cylinderEnd = cylinder.second->getEndpoint();
                        // Calculate the axis direction of the cylinder
                        Eigen::Vector3d cylinderAxis = (cylinderEnd - cylinderStart).normalized();

                        // Calculate the vector from the start point of the cylinder to the sphere center
                        Eigen::Vector3d startToCenter = center - cylinderStart;

                        // Project the vector onto the cylinder axis to find the closest point on the axis to the sphere center
                        double projection = startToCenter.dot(cylinderAxis);
                        projection = std::clamp(projection, 0.0, 1.0);
                        Eigen::Vector3d closestPoint = cylinderStart + projection * cylinderAxis;
                        double distanceToaxis = (closestPoint - center).norm();
                        double overlap = 0;
                        if (distanceToaxis < cylinder.second->getRadius())
                        {
                            if (!manager->getCylinderwallProperties(cylinder.second->getType())->getHollow())
                            {
                                validPosition = false;
                                break;
                            }
                            else
                            {
                                overlap = distanceToaxis - (cylinder.second->getRadius() - sphereRadius);
                            }
                        }
                        else
                        {
                            overlap = cylinder.second->getRadius() + sphereRadius - distanceToaxis;
                        }
                        if (overlap > 0)
                        {
                            validPosition = false;
                            break;
                        }
                    }
                }

                if (validPosition)
                {
                    for (const auto &existparticle : sphereparticles)
                    {
                        double distance = (existparticle.second->getPosition() - center).norm();
                        double sphere1Radius = manager->getSphereProperties(existparticle.second->getType())->getRadius();
                        double overlap_pp = (sphere1Radius + sphereRadius) - distance;

                        if (overlap_pp > 0)
                        {
                            validPosition = false;
                            break;
                        }
                    }
                }
                if (validPosition)
                {
                    for (const auto &bond : fiberbonds)
                    {
                        Eigen::Vector3d onenode1position = fibershpereparticles[bond.second->getNode1()]->getPosition();
                        Eigen::Vector3d onenode2position = fibershpereparticles[bond.second->getNode2()]->getPosition();

                        double fiberRadius = manager->getFiberProperties(bond.second->getType())->getRadius();
                        double t;
                        Eigen::Vector3d projection;
                        SphereCylinderBond::computeOverlap(onenode1position, onenode2position, center, t, projection);
                        double distance = (center - projection).norm();

                        double overlap = (sphereRadius + fiberRadius) - distance;
                        if (overlap > 0)
                        {
                            validPosition = false;
                            break;
                        }

                        if (!validPosition)
                            break;
                    }
                }
            } while (!validPosition && count_number < 1000);
            if (validPosition)
            {
                auto sp = std::make_unique<SphereParticle>(id_index, PropertyTypeID(categoryId, subTypeId), state, manager, Eigen::Vector3d(x, y, z));
                sphereparticles.insert({sp->getId(), std::move(sp)});
                id_index++;
            }
            else
            {
                gernerateSphereFlag = false;
                failedSpheres = count - i;
                // Use an ostringstream to format the data as a string
                std::ostringstream oss;
                oss << type << " "
                    << categoryId << " "
                    << subTypeId << " "
                    << state << " "
                    << failedSpheres << " "
                    << xmin << " "
                    << xmax << " "
                    << ymin << " "
                    << ymax << " "
                    << zmin << " "
                    << zmax;

                gernerateInfor[PropertyTypeID(categoryId, subTypeId)] = oss.str();

                break;
            }
        }
        if (gernerateSphereFlag)
        {
            // Use an ostringstream to format the data as a string
            std::ostringstream oss;
            oss << type << " "
                << categoryId << " "
                << subTypeId << " "
                << state << " "
                << failedSpheres << " "
                << xmin << " "
                << xmax << " "
                << ymin << " "
                << ymax << " "
                << zmin << " "
                << zmax;

            gernerateInfor[PropertyTypeID(categoryId, subTypeId)] = oss.str();
        }
        else
        {
            std::cout << "Can't generate all spheres, number of remaining spheres: " << failedSpheres << std::endl;
        }
    }
    else if (type == "FIBER")
    {
        int categoryId, subTypeId, startstate, endstate, count;
        double xmin, xmax, ymin, ymax, zmin, zmax;
        iss >> categoryId >> subTypeId >> startstate >> endstate >> count >> xmin >> xmax >> ymin >> ymax >> zmin >> zmax;

        const auto &fiberproperties = manager->getFiberProperties(PropertyTypeID(categoryId, subTypeId));
        // FiberProperties copyfiberproperties = *fiberproperties;

        int bond_id_index = fiberbonds.size();
        int sphere_id_index = fibershpereparticles.size();
        int fiber_id_index = fibers.size();
        gernerateFiberFlag = true;
        int failedFibers = 0;
        double elementlength = fiberproperties->getElementlength();
        double fiberRadius = fiberproperties->getRadius();
        std::random_device rd;
        // unsigned int seed = 300; // 固定种子值
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> disX(xmin, xmax);
        std::uniform_real_distribution<> disY(ymin, ymax);
        std::uniform_real_distribution<> disZ(zmin, zmax);

        std::uniform_real_distribution<> oriX(-1, 1);
        std::uniform_real_distribution<> oriY(-1, 1);
        std::uniform_real_distribution<> oriZ(-1, 1);

        for (int i = 0; i < count; ++i)
        {
            double x, y, z, ox, oy, oz;
            bool validPosition = false;
            int count_number = 0;
            std::vector<Eigen::Vector3d> tempsphereposition;
            tempsphereposition.resize(fiberproperties->getNodenumber());
            do
            {
                count_number++;
                validPosition = true;
                // Random starting position of the fiber
                x = disX(gen);
                y = disY(gen);
                z = disZ(gen);
                ox = oriX(gen);
                oy = oriY(gen);
                oz = oriZ(gen);

                Eigen::Vector3d startPos(x, y, z);
                // Random orientation of the fiber
                Eigen::Vector3d orientation(ox, oy, oz);
                orientation.normalize();
                for (int k = 0; k < fiberproperties->getNodenumber(); ++k)
                {

                    Eigen::Vector3d position = startPos + k * elementlength * orientation;

                    tempsphereposition[k] = position;
                }

                for (const auto &plane : planewalls)
                {
                    // Retrieve the plane wall's point and normal
                    const Eigen::Vector3d &planePoint = plane.second->getCorner1();
                    const Eigen::Vector3d &planeNormal = plane.second->getNormal();
                    for (int j = 0; j < tempsphereposition.size() - 1; ++j)
                    {
                        double distanceToPlane;
                        SphereCylinderBond::computeOverlap(tempsphereposition[j], tempsphereposition[j + 1], plane.second, distanceToPlane);

                        double overlap_pw = fiberRadius - distanceToPlane;
                        if (overlap_pw > 0)
                        {
                            validPosition = false;
                            break;
                        }
                    }
                    if (!validPosition)
                        break;
                }
                if (validPosition)
                {
                    for (const auto &cylinder : cylinderwalls)
                    {
                        Eigen::Vector3d cylinderStart = cylinder.second->getStartpoint();
                        Eigen::Vector3d cylinderEnd = cylinder.second->getEndpoint();
                        // Calculate the axis direction of the cylinder
                        Eigen::Vector3d cylinderAxis = (cylinderEnd - cylinderStart).normalized();
                        for (int j = 0; j < tempsphereposition.size() - 1; ++j)
                        {
                            // Calculate the vector from the start point of the cylinder to the sphere center
                            Eigen::Vector3d startToCenter = tempsphereposition[j] - cylinderStart;
                            Eigen::Vector3d endToCenter = tempsphereposition[j + 1] - cylinderStart;

                            // Project the vector onto the cylinder axis to find the closest point on the axis to the sphere center
                            double projection1 = startToCenter.dot(cylinderAxis);
                            projection1 = std::clamp(projection1, 0.0, 1.0);
                            Eigen::Vector3d closestPoint1 = cylinderStart + projection1 * cylinderAxis;
                            double distanceToaxis1 = (closestPoint1 - tempsphereposition[j]).norm();

                            double projection2 = endToCenter.dot(cylinderAxis);
                            projection2 = std::clamp(projection2, 0.0, 1.0);
                            Eigen::Vector3d closestPoint2 = cylinderStart + projection2 * cylinderAxis;
                            double distanceToaxis2 = (closestPoint2 - tempsphereposition[j + 1]).norm();

                            double overlap = 0;
                            if (distanceToaxis1 < cylinder.second->getRadius() || distanceToaxis2 < cylinder.second->getRadius())
                            {
                                if (!manager->getCylinderwallProperties(cylinder.second->getType())->getHollow())
                                {
                                    validPosition = false;
                                    break;
                                }
                                else
                                {
                                    double overlap1 = distanceToaxis1 - (cylinder.second->getRadius() - fiberRadius);
                                    double overlap2 = distanceToaxis2 - (cylinder.second->getRadius() - fiberRadius);
                                    overlap = overlap1 > overlap2 ? overlap1 : overlap2;
                                }
                            }
                            else
                            {
                                double overlap1 = cylinder.second->getRadius() + fiberRadius - distanceToaxis1;
                                double overlap2 = cylinder.second->getRadius() + fiberRadius - distanceToaxis2;
                                overlap = overlap1 > overlap2 ? overlap1 : overlap2;
                            }
                            if (overlap > 0)
                            {
                                validPosition = false;
                                break;
                            }
                        }
                        if (!validPosition)
                            break;
                    }
                }

                // Check overlap with existing spheres
                if (validPosition)
                {
                    for (const auto &sphere : sphereparticles)
                    {
                        // Retrieve the sphere's center and radius
                        Eigen::Vector3d sphereCenter = sphere.second->getPosition();
                        double sphereRadius = manager->getSphereProperties(sphere.second->getType())->getRadius();
                        for (int j = 0; j < tempsphereposition.size() - 1; ++j)
                        {
                            Eigen::Vector3d projection;
                            double t;
                            SphereCylinderBond::computeOverlap(tempsphereposition[j], tempsphereposition[j + 1], sphereCenter, t, projection);
                            double distance = (projection - sphereCenter).norm();

                            double overlap = (fiberRadius + sphereRadius) - distance;
                            if (overlap > 0)
                            {
                                validPosition = false;
                                break;
                            }
                        }
                        if (!validPosition)
                            break;
                    }
                }

                // Check overlap with existing fibers
                if (validPosition)
                {
                    for (const auto &bond : fiberbonds)
                    {
                        Eigen::Vector3d onenode1position = fibershpereparticles[bond.second->getNode1()]->getPosition();
                        Eigen::Vector3d onenode2position = fibershpereparticles[bond.second->getNode2()]->getPosition();

                        double fiber1Radius = manager->getFiberProperties(bond.second->getType())->getRadius();
                        for (int j = 0; j < tempsphereposition.size() - 1; ++j)
                        {
                            double s = 0, t = 0;
                            SphereCylinderBond::computeOverlap(tempsphereposition[j], tempsphereposition[j + 1], s,
                                                               onenode1position, onenode2position, t);
                            Eigen::Vector3d projection1 = (1 - s) * tempsphereposition[j] + s * tempsphereposition[j + 1];
                            Eigen::Vector3d projection2 = (1 - t) * onenode1position + t * onenode2position;
                            double distance = (projection1 - projection2).norm();

                            double overlap = (fiberRadius + fiber1Radius) - distance;
                            if (overlap > 0)
                            {
                                validPosition = false;
                                break;
                            }
                        }
                        if (!validPosition)
                            break;
                    }
                }

            } while (!validPosition && count_number < 50000);
            if (validPosition)
            {
                auto fiber = std::make_unique<Fiber>(fiber_id_index, PropertyTypeID(categoryId, subTypeId), startstate, endstate, manager);
                std::vector<int> spherecylinderbondID;
                for (int k = 0; k < fiberproperties->getNodenumber(); ++k)
                {
                    auto createSphereParticle = [&](int state)
                    {
                        Eigen::Vector3d position = tempsphereposition[k];

                        return std::make_unique<SphereParticle>(
                            sphere_id_index + k, PropertyTypeID(categoryId, subTypeId), state, manager, position);
                    };

                    // Determine the state of the SphereParticle
                    int state = (k == 0) ? startstate : (k == fiberproperties->getNodenumber() - 1) ? endstate
                                                                                                    : 1;
                    auto sp = createSphereParticle(state);
                    fibershpereparticles.insert({sp->getId(), std::move(sp)});

                    // Determine the SphereCylinderBond neighbors
                    int neighborElement1 = (k == 0) ? -1 : bond_id_index + k - 1;
                    int neighborElement2 = (k == fiberproperties->getNodenumber() - 2) ? -1 : bond_id_index + k + 1;
                    // Create SphereCylinderBond except for the last element
                    if (k != fiberproperties->getNodenumber() - 1)
                    {
                        auto bond = std::make_unique<SphereCylinderBond>(
                            bond_id_index + k,
                            PropertyTypeID(categoryId, subTypeId),
                            manager,
                            fiber_id_index,
                            sphere_id_index + k,
                            sphere_id_index + k + 1,
                            neighborElement1,
                            neighborElement2);
                        fiberbonds.insert({bond->getId(), std::move(bond)});
                        spherecylinderbondID.push_back(bond_id_index + k);
                    }
                }
                fiber->setSphereCylinderBond(spherecylinderbondID);
                fibers.push_back(std::move(fiber));
                sphere_id_index += fiberproperties->getNodenumber();
                fiber_id_index++;
                bond_id_index += fiberproperties->getNodenumber() - 1;
            }
            else
            {
                gernerateFiberFlag = false;
                failedFibers = count - i;
                // Use an ostringstream to format the data as a string
                std::ostringstream oss;
                oss << type << " "
                    << categoryId << " "
                    << subTypeId << " "
                    << startstate << " "
                    << endstate << " "
                    << failedFibers << " "
                    << xmin << " "
                    << xmax << " "
                    << ymin << " "
                    << ymax << " "
                    << zmin << " "
                    << zmax;

                gernerateInfor[PropertyTypeID(categoryId, subTypeId)] = oss.str();

                break;
            }
        }
        if (gernerateFiberFlag)
        {
            // Use an ostringstream to format the data as a string
            std::ostringstream oss;
            oss << type << " "
                << categoryId << " "
                << subTypeId << " "
                << startstate << " "
                << endstate << " "
                << failedFibers << " "
                << xmin << " "
                << xmax << " "
                << ymin << " "
                << ymax << " "
                << zmin << " "
                << zmax;

            gernerateInfor[PropertyTypeID(categoryId, subTypeId)] = oss.str();
        }
        else
        {
            std::cout << "Can't generate all fibers, number of remaining fibers: " << failedFibers << std::endl;
        }
    }
}
void DEMProperties::parseSpecificParticle(std::istringstream &iss)
{
    std::string type;
    iss >> type;
    if (type == "SPHERE")
    {
        int id, categoryId, subTypeId, state;
        double x, y, z, vx, vy, vz;
        iss >> id >> categoryId >> subTypeId >> state >> x >> y >> z >> vx >> vy >> vz;
        // Create and add the specific particle to your simulation
        auto sp = std::make_unique<SphereParticle>(id, PropertyTypeID(categoryId, subTypeId), state, manager, Eigen::Vector3d(x, y, z), Eigen::Vector3d(vx, vy, vz));

        sphereparticles.insert({sp->getId(), std::move(sp)});
    }
    else if (type == "FIBER")
    {
        int id, categoryId, subTypeId, startstate, endstate;
        double startx, starty, startz, vx, vy, vz;
        iss >> id >> categoryId >> subTypeId >> startstate >> endstate >> startx >> starty >> startz >> vx >> vy >> vz;
        std::string lauout;
        iss >> lauout;
        auto fiber = std::make_unique<Fiber>(id, PropertyTypeID(categoryId, subTypeId), startstate, endstate, manager);
        const auto &fiberproperties = manager->getFiberProperties(PropertyTypeID(categoryId, subTypeId));
        FiberProperties copyfiberproperties = *fiberproperties;
        int bond_id_index = fiberbonds.size();
        int sphere_id_index = fibershpereparticles.size();
        std::vector<int> spherecylinderbondID;
        for (int k = 0; k < copyfiberproperties.getNodenumber(); ++k)
        {
            // Lambda function to create a SphereParticle based on layout
            auto createSphereParticle = [&](int state, const std::string &layout)
            {
                Eigen::Vector3d position;
                if (layout == "X")
                {
                    position = Eigen::Vector3d(startx + k * copyfiberproperties.getElementlength(), starty, startz);
                }
                else if (layout == "Y")
                {
                    position = Eigen::Vector3d(startx, starty + k * copyfiberproperties.getElementlength(), startz);
                }
                else if (layout == "Z")
                {
                    position = Eigen::Vector3d(startx, starty, startz + k * copyfiberproperties.getElementlength());
                }
                return std::make_unique<SphereParticle>(
                    sphere_id_index + k, PropertyTypeID(categoryId, subTypeId), state, manager, position, Eigen::Vector3d(vx, vy, vz));
            };

            // Determine the state of the SphereParticle
            int state = (k == 0) ? startstate : (k == copyfiberproperties.getNodenumber() - 1) ? endstate
                                                                                               : 1;
            auto sp = createSphereParticle(state, lauout);
            fibershpereparticles.insert({sp->getId(), std::move(sp)});

            // Determine the SphereCylinderBond neighbors
            int neighborElement1 = (k == 0) ? -1 : bond_id_index + k - 1;
            int neighborElement2 = (k == copyfiberproperties.getNodenumber() - 2) ? -1 : bond_id_index + k + 1;

            // Create SphereCylinderBond except for the last element
            if (k != copyfiberproperties.getNodenumber() - 1)
            {
                auto bond = std::make_unique<SphereCylinderBond>(
                    bond_id_index + k,
                    PropertyTypeID(categoryId, subTypeId),
                    manager,
                    id,
                    sphere_id_index + k,
                    sphere_id_index + k + 1,
                    neighborElement1,
                    neighborElement2);
                fiberbonds.insert({bond->getId(), std::move(bond)});
                spherecylinderbondID.push_back(bond_id_index + k);
            }
        }
        fiber->setSphereCylinderBond(spherecylinderbondID);
        fibers.push_back(std::move(fiber));
    }
}
void DEMProperties::parseFiberParticle(std::istringstream &iss)
{
    std::string type;
    if (iss >> type)
    {
        if (type == "SPHERE")
        {
            int id, categoryId, subtypeId, state;
            Eigen::Vector3d position, velocity, omega;
            iss >> id >> categoryId >> subtypeId >> state;
            if (!parseVector3d(iss, position) || !parseVector3d(iss, velocity) || !parseVector3d(iss, omega))
            {
                throw std::runtime_error("Error parsing sphere vectors.");
            }
            auto sphere = std::make_unique<SphereParticle>(id, PropertyTypeID(categoryId, subtypeId), state, manager, position, velocity, omega);
            fibershpereparticles.insert({id, std::move(sphere)});
        }
        else if (type == "Bond")
        {

            int id, categoryId, subtypeId, fiberid, node1, node2, neighborelement1, neighborelement2;
            double energyDissipation;
            Eigen::Vector3d tangentialforce, twisttorque, bendtorque;
            iss >> id >> categoryId >> subtypeId >> fiberid >> node1 >> node2 >> neighborelement1 >> neighborelement2;
            iss >> energyDissipation;
            if (!parseVector3d(iss, tangentialforce) || !parseVector3d(iss, twisttorque) || !parseVector3d(iss, bendtorque))
            {
                throw std::runtime_error("Error parsing spherecylinder bond vectors.");
            }
            auto bond = std::make_unique<SphereCylinderBond>(id, PropertyTypeID(categoryId, subtypeId), manager, fiberid, node1, node2, neighborelement1,
                                                             neighborelement2, energyDissipation, tangentialforce, twisttorque, bendtorque);
            fiberbonds.insert({id, std::move(bond)});
        }
        else
        {
            int id, categoryId, subtypeId, startstate, endstate;
            id = std::stoi(type);
            iss >> categoryId >> subtypeId >> startstate >> endstate;
            PropertyTypeID fibertype = PropertyTypeID(categoryId, subtypeId);
            std::vector<int> spherecylinderbond(manager->getFiberProperties(fibertype)->getNodenumber() - 1, 0);
            for (int i = 0; i < spherecylinderbond.size(); ++i)
            {
                iss >> spherecylinderbond[i];
            }

            auto fiber = std::make_unique<Fiber>(id, fibertype, startstate, endstate, manager);
            fiber->setSphereCylinderBond(spherecylinderbond);
            fibers.push_back(std::move(fiber));
        }
    }
}
void DEMProperties::generateRemainingParticles()
{
    if (!gernerateFiberFlag || !gernerateSphereFlag)
    {

        for (const auto &sub : gernerateInfor)
        {

            std::istringstream iss(sub.second);
            parseRandomParticle(iss);
        }
    }
}
double DEMProperties::getAverageVelocity()
{
    double totalVelocity = 0;

    for (const auto &sphere : sphereparticles)
    {
        totalVelocity += sphere.second->getVelocity().norm();
    }
    for (const auto &sphere : fibershpereparticles)
    {
        totalVelocity += sphere.second->getVelocity().norm();
    }
    double averageVelocity = 0;
    if (sphereparticles.size() + fibershpereparticles.size() != 0)
    {
        averageVelocity = totalVelocity / (sphereparticles.size() + fibershpereparticles.size());
    }
    return averageVelocity;
}
void DEMProperties::initialSimulation()
{
    contactforce->addParticleProperties(manager);
    contactforce->setEnergydissipation(contactEnergydissipation);
    double maxRadius = std::numeric_limits<double>::min();
    double maxYoungmodulus = std::numeric_limits<double>::min();
    double maxBondmodulus = std::numeric_limits<double>::min();
    double minBondlength = std::numeric_limits<double>::max();
    double mindensity = std::numeric_limits<double>::max();
    double minRadius = std::numeric_limits<double>::max();
    double mindensityL = std::numeric_limits<double>::max();

    for (const auto &particleproperties : manager->getParticleProperties())
    {
        if (auto &sphereProps = std::dynamic_pointer_cast<SphereProperties>(particleproperties.second))
        {
            double radius = sphereProps->getRadius();
            if (radius > maxRadius)
            {
                maxRadius = radius;
            }
            if (radius < minRadius)
            {
                minRadius = radius;
            }
            double Youngmodulus = sphereProps->getYoungModulus();
            if (Youngmodulus > maxYoungmodulus)
            {
                maxYoungmodulus = Youngmodulus;
            }
            double density = sphereProps->getDensity();
            if (density < mindensity)
            {
                mindensity = density;
            }
        }
        else if (auto &fiberProps = std::dynamic_pointer_cast<FiberProperties>(particleproperties.second))
        {
            double radius = fiberProps->getRadius();
            if (radius > maxRadius)
            {
                maxRadius = radius;
            }
            if (radius < minRadius)
            {
                minRadius = radius;
            }
            double Youngmodulus = fiberProps->getYoungModulus();
            if (Youngmodulus > maxYoungmodulus)
            {
                maxYoungmodulus = Youngmodulus;
            }
            double Bondmodulus = fiberProps->getNormalmodulus();
            if (Bondmodulus > maxBondmodulus)
            {
                maxBondmodulus = Bondmodulus;
            }
            Bondmodulus = fiberProps->getShearmodulus();
            if (Bondmodulus > maxBondmodulus)
            {
                maxBondmodulus = Bondmodulus;
            }
            Bondmodulus = fiberProps->getTwistmodulus();
            if (Bondmodulus > maxBondmodulus)
            {
                maxBondmodulus = Bondmodulus;
            }
            Bondmodulus = fiberProps->getBendingmodulus();
            if (Bondmodulus > maxBondmodulus)
            {
                maxBondmodulus = Bondmodulus;
            }

            double Bondlength = fiberProps->getElementlength();
            if (Bondlength < minBondlength)
            {
                minBondlength = Bondlength;
            }

            double densityL = fiberProps->getNodemass() / Bondlength;
            if (densityL < mindensityL)
            {
                mindensityL = densityL;
            }

            double density = fiberProps->getDensity();
            if (density < mindensity)
            {
                mindensity = density;
            }
        }
    }

    double temp_t = 0.1 * minBondlength * sqrt(mindensity / (maxBondmodulus * PI * maxRadius * maxRadius));
    if (timestep > temp_t)
    {
        timestep = temp_t;
    }
    temp_t = (0.1 * 2 * PI * minRadius) / (0.1631 * 0.3 + 0.8766) * sqrt(1.3 * mindensity / (2 * maxYoungmodulus));
    if (timestep > temp_t)
    {
        timestep = temp_t;
    }
    gridbasedcontactdetection->initial(simulationdimensions.x(), simulationdimensions.y(), simulationdimensions.z(), 2.5 * maxRadius);
    for (const auto &planewall : planewalls)
    {
        planewall.second->setMass(manager);
    }

    for (const auto &cylinderwall : cylinderwalls)
    {
        cylinderwall.second->setMass(manager);
    }

    for (const auto &connectedgeometry : connectedgeometrys)
    {
        double mass = 0;
        for (auto &geometry : connectedgeometry.second->getGeometrys())
        {
            if (std::get<0>(geometry) == "SPHERE")
            {
                mass += manager->getSphereProperties(sphereparticles[std::get<1>(geometry)]->getType())->getMass();
            }
            else if (std::get<0>(geometry) == "FIBER")
            {
                mass += manager->getFiberProperties(fibers[std::get<1>(geometry)]->getType())->getNodemass() * manager->getFiberProperties(fibers[std::get<1>(geometry)]->getType())->getNodenumber();
            }
            else if (std::get<0>(geometry) == "PLANEWALL")
            {

                mass += planewalls[std::get<1>(geometry)]->getMass();
            }
            else if (std::get<0>(geometry) == "CYLINDERCONTAINER")
            {
                mass += cylinderwalls[std::get<1>(geometry)]->getMass();
            }
        }
        connectedgeometry.second->setMass(mass);
    }
}
void DEMProperties::handleCollisions()
{
    std::unordered_map<int, std::unordered_set<int>> sphere_sphere_contact_paris;
    std::unordered_map<int, std::unordered_set<int>> sphere_fiber_contact_paris;

    std::unordered_map<int, std::unordered_set<int>> planewall_sphere_contact_paris;
    std::unordered_map<int, std::unordered_set<int>> planewall_fiber_contact_paris;

    std::unordered_map<int, std::unordered_set<int>> fiber_fiber_contact_paris;

    std::unordered_map<int, std::unordered_set<int>> cylindercontainer_sphere_contact_paris;
    std::unordered_map<int, std::unordered_set<int>> cylindercontainer_fiber_contact_paris;

    gridbasedcontactdetection->ParticleBroadPhase(sphereparticles, fiberbonds, fibershpereparticles,
                                                  sphere_sphere_contact_paris, sphere_fiber_contact_paris, fiber_fiber_contact_paris);
    for (const auto &contact_list : sphere_sphere_contact_paris)
    {
        auto id1 = contact_list.first;
        for (auto id2 : contact_list.second)
        {
            auto &sp1 = sphereparticles[id1];
            auto &sp2 = sphereparticles[id2];

            contactforce->computeSphereSphereForce(sp1, sp2, timestep);
        }
    }

    for (const auto &contact_list : sphere_fiber_contact_paris)
    {
        auto id1 = contact_list.first;
        for (auto id2 : contact_list.second)
        {
            auto &sp1 = sphereparticles[id1];
            auto &sc2 = fiberbonds[id2];

            auto &node1 = fibershpereparticles[sc2->getNode1()];
            auto &node2 = fibershpereparticles[sc2->getNode2()];

            contactforce->computeSphereFiberForce(sp1, sc2, node1, node2, timestep);
        }
    }

    for (const auto &contact_list : fiber_fiber_contact_paris)
    {
        auto id1 = contact_list.first;
        for (auto id2 : contact_list.second)
        {
            const auto &sc1 = fiberbonds[id1];
            const auto &sc2 = fiberbonds[id2];

            auto &sconenode1 = fibershpereparticles[sc1->getNode1()];
            auto &sconenode2 = fibershpereparticles[sc1->getNode2()];
            auto &sctwonode1 = fibershpereparticles[sc2->getNode1()];
            auto &sctwonode2 = fibershpereparticles[sc2->getNode2()];

            contactforce->computeFiberFiberForce(sc1, sconenode1, sconenode2, sc2, sctwonode1, sctwonode2, timestep);
        }
    }

    gridbasedcontactdetection->planewallBroadPhase(planewalls, planewall_sphere_contact_paris, planewall_fiber_contact_paris);
    for (const auto &contact_list : planewall_sphere_contact_paris)
    {
        auto wall_id = contact_list.first;
        for (auto id : contact_list.second)
        {
            auto &wall = planewalls[wall_id];
            auto &sphere = sphereparticles[id];

            contactforce->computePlaneWallSphereForce(wall, sphere, timestep);
        }
    }
    for (const auto &contact_list : planewall_fiber_contact_paris)
    {
        auto wall_id = contact_list.first;
        for (auto id2 : contact_list.second)
        {
            auto &wall = planewalls[wall_id];
            auto &sc2 = fiberbonds[id2];

            auto &node1 = fibershpereparticles[sc2->getNode1()];
            auto &node2 = fibershpereparticles[sc2->getNode2()];

            contactforce->computePlaneWallFiberForce(wall, sc2, node1, node2, timestep);
        }
    }

    gridbasedcontactdetection->cylinderwallBroadPhase(cylinderwalls, cylindercontainer_sphere_contact_paris, cylindercontainer_fiber_contact_paris);
    for (const auto &contact_list : cylindercontainer_sphere_contact_paris)
    {
        auto wall_id = contact_list.first;
        for (auto id : contact_list.second)
        {
            auto &wall = cylinderwalls[wall_id];
            auto &sphere = sphereparticles[id];
            contactforce->computeCylinderWallSphereForce(wall, sphere, timestep);
        }
    }
    for (const auto &contact_list : cylindercontainer_fiber_contact_paris)
    {
        auto wall_id = contact_list.first;
        for (auto id2 : contact_list.second)
        {
            auto &wall = cylinderwalls[wall_id];
            auto &sc2 = fiberbonds[id2];

            auto &node1 = fibershpereparticles[sc2->getNode1()];
            auto &node2 = fibershpereparticles[sc2->getNode2()];

            contactforce->computeCylinderWallFiberForce(wall, sc2, node1, node2, timestep);
        }
    }
}
void DEMProperties::bondforce()
{
    bondEnergydissipation = 0;
    for (const auto &particle : fiberbonds)
    {
        auto &node1 = fibershpereparticles[particle.second->getNode1()];
        auto &node2 = fibershpereparticles[particle.second->getNode2()];
        particle.second->updateBond(node1, node2, timestep);
        bondEnergydissipation += particle.second->getEnergydissipation();
    }
}
void DEMProperties::applyExternalForces()
{
    for (auto &sphere : sphereparticles)
    {

        sphere.second->addForce(globalforces);
    }
    for (auto &indexforce : specificforces)
    {
        sphereparticles[indexforce.first]->addForce(indexforce.second);
    }
}
void DEMProperties::motion()
{
    // motion particles
    for (auto &sphere : sphereparticles)
    {
        if (sphere.second->getState() == 1)
        {
            double mass = manager->getSphereProperties(sphere.second->getType())->getMass();
            double moment_of_inertia = manager->getSphereProperties(sphere.second->getType())->getMomentOfInertia();

            sphere.second->updateVelocity(timestep, gravity, mass);
            sphere.second->updateOmega(timestep, moment_of_inertia);
            sphere.second->updatePosition(timestep);
        }
        sphere.second->resetForce();
        sphere.second->resetTorque();
    }

    for (auto &sphere : fibershpereparticles)
    {
        if (sphere.second->getState() == 1)
        {
            double mass = manager->getFiberProperties(sphere.second->getType())->getNodemass();
            double moment_of_inertia = manager->getFiberProperties(sphere.second->getType())->getNodemomentofinertia();

            sphere.second->updateVelocity(timestep, gravity, mass);
            sphere.second->updateOmega(timestep, moment_of_inertia);
            sphere.second->updatePosition(timestep);
        }
        sphere.second->resetForce();
        sphere.second->resetTorque();
    }
    contactEnergydissipation = contactforce->getEnergydissipation();
    contactforce->updateContactInformation();

    for (auto &planewall : planewalls)
    {
        planewall.second->resetForce();
    }
    for (auto &cylinderwall : cylinderwalls)
    {
        cylinderwall.second->resetForce();
    }
}
void DEMProperties::initial_task()
{
    if (timestep > taskTimestep)
    {
        timestep = taskTimestep;
    }
}
void DEMProperties::dotask()
{
    for (auto &subtask : dynamicBoundaryTask)
    {
        Eigen::Vector3d resultantForce = Eigen::Vector3d::Zero();
        Eigen::Vector3d v = connectedgeometrys[subtask.first]->getVelocity();
        Eigen::Vector3d center = connectedgeometrys[subtask.first]->getCenterOfMass();
        double mass = connectedgeometrys[subtask.first]->getMass();
        Eigen::Vector3d displancement;

        for (auto &geometry : connectedgeometrys[subtask.first]->getGeometrys())
        {
            if (std::get<0>(geometry) == "SPHERE")
            {
                resultantForce += sphereparticles[std::get<1>(geometry)]->getForce();
                sphereparticles[std::get<1>(geometry)]->resetForce();
                sphereparticles[std::get<1>(geometry)]->resetTorque();

                sphereparticles[std::get<1>(geometry)]->setVelocity(v);
                // sphereparticles[std::get<1>(geometry)]->setOmega(v);
            }
            else if (std::get<0>(geometry) == "FIBER")
            {
                for (int i = 0; i < fibers[std::get<1>(geometry)]->getSphereCylinderBond().size(); ++i)
                {
                    if (i == 0)
                    {
                        int node1Id = fiberbonds[fibers[std::get<1>(geometry)]->getSphereCylinderBond()[i]]->getNode1();
                        resultantForce += fibershpereparticles[node1Id]->getForce();
                        fibershpereparticles[node1Id]->resetForce();
                        fibershpereparticles[node1Id]->resetTorque();
                    }
                    int node2Id = fiberbonds[fibers[std::get<1>(geometry)]->getSphereCylinderBond()[i]]->getNode2();
                    resultantForce += fibershpereparticles[node2Id]->getForce();
                    fibershpereparticles[node2Id]->resetForce();
                    fibershpereparticles[node2Id]->resetTorque();
                }
            }
            else if (std::get<0>(geometry) == "PLANEWALL")
            {
                resultantForce += planewalls[std::get<1>(geometry)]->getForce();
                planewalls[std::get<1>(geometry)]->resetForce();
            }
            else if (std::get<0>(geometry) == "CYLINDERCONTAINER")
            {
                resultantForce += cylinderwalls[std::get<1>(geometry)]->getForce();
                cylinderwalls[std::get<1>(geometry)]->resetForce();
            }
        }

        if (auto &springoscillator = std::dynamic_pointer_cast<SpringOscillator>(subtask.second))
        {
            resultantForce += springoscillator->getForce(center.y(), v.y());
            v += (resultantForce / mass + gravity) * timestep;
            displancement = v * timestep;
        }
        else if (auto &constantmove = std::dynamic_pointer_cast<ConstantMove>(subtask.second))
        {
            v = constantmove->getVelocity();
            double period = constantmove->getPeriod();
           
            double time_since_start = fmod(currentTime, period);
            if (time_since_start >= period / 2)
            {
                v = -v;
            }
            

            displancement = v * timestep;
        }

        if (!connectedgeometrys[subtask.first]->getDegreeOfFreedom(ConnectedGeometry::DegreeOfFreedom::X))
        {
            displancement.x() = 0;
        }
        if (!connectedgeometrys[subtask.first]->getDegreeOfFreedom(ConnectedGeometry::DegreeOfFreedom::Y))
        {
            displancement.y() = 0;
        }
        if (!connectedgeometrys[subtask.first]->getDegreeOfFreedom(ConnectedGeometry::DegreeOfFreedom::Z))
        {
            displancement.z() = 0;
        }
        connectedgeometrys[subtask.first]->setCenterOfMass(center + displancement);
        connectedgeometrys[subtask.first]->setVelocity(v);

        for (auto &geometry : connectedgeometrys[subtask.first]->getGeometrys())
        {
            if (std::get<0>(geometry) == "SPHERE")
            {
                Eigen::Vector3d initialPosition = sphereparticles[std::get<1>(geometry)]->getPosition();
                sphereparticles[std::get<1>(geometry)]->setPosition(initialPosition + displancement);
            }
            else if (std::get<0>(geometry) == "FIBER")
            {
                for (int i = 0; i < fibers[std::get<1>(geometry)]->getSphereCylinderBond().size(); ++i)
                {
                    if (i == 0)
                    {
                        int node1Id = fiberbonds[fibers[std::get<1>(geometry)]->getSphereCylinderBond()[i]]->getNode1();
                        Eigen::Vector3d initialPosition = fibershpereparticles[node1Id]->getPosition();
                        fibershpereparticles[node1Id]->setPosition(initialPosition + displancement);
                    }
                    int node2Id = fiberbonds[fibers[std::get<1>(geometry)]->getSphereCylinderBond()[i]]->getNode2();
                    Eigen::Vector3d initialPosition = fibershpereparticles[node2Id]->getPosition();
                    fibershpereparticles[node2Id]->setPosition(initialPosition + displancement);
                }
            }
            else if (std::get<0>(geometry) == "PLANEWALL")
            {
                Eigen::Vector3d initialPosition1 = planewalls[std::get<1>(geometry)]->getCorner1();
                Eigen::Vector3d initialPosition2 = planewalls[std::get<1>(geometry)]->getCorner2();
                Eigen::Vector3d initialPosition3 = planewalls[std::get<1>(geometry)]->getCorner3();
                planewalls[std::get<1>(geometry)]->setCorner1(initialPosition1 + displancement);
                planewalls[std::get<1>(geometry)]->setCorner2(initialPosition2 + displancement);
                planewalls[std::get<1>(geometry)]->setCorner3(initialPosition3 + displancement);
                planewalls[std::get<1>(geometry)]->setVelocity(v);
            }
            else if (std::get<0>(geometry) == "CYLINDERCONTAINER")
            {
                Eigen::Vector3d initialPosition1 = cylinderwalls[std::get<1>(geometry)]->getStartpoint();
                Eigen::Vector3d initialPosition2 = cylinderwalls[std::get<1>(geometry)]->getEndpoint();
                cylinderwalls[std::get<1>(geometry)]->setStartpoint(initialPosition1 + displancement);
                cylinderwalls[std::get<1>(geometry)]->setEndpoint(initialPosition2 + displancement);
            }
        }
    }
}
