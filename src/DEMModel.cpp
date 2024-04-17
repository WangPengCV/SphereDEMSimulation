#include "DEMModel.h"

// Initialize DEMModel
DEMModel::DEMModel(const std::string &filename)
{
    DEMproperties = std::make_shared<DEMProperties>();

    std::string stateFolder = "DEMProperties";
    if (canResumeFromState(stateFolder))
    {
        std::string lastStateFile = findLastStateFile(stateFolder);
        std::cout << "Resuming from state: " << lastStateFile << std::endl;
        DEMproperties->loadFromFile(lastStateFile); // Make sure this function loads all needed state
        iter_num = extractNumberFromFilename(lastStateFile);
    }
    else
    {
        DEMproperties->loadFromFile(filename);
        iter_num = -1;
    }
    vis = std::make_shared<Visualization>(*DEMproperties);
    DEMproperties->initialSimulation();
}

void DEMModel::runSimulation()
{

    double currentTime = DEMproperties->getCurrentTime();
    double totalTime = DEMproperties->getTotalTime();
    double timeStep = DEMproperties->getTimestep();
    double criticalSpeed = DEMproperties->getCriticalSpeed();
    int taskShowInterval = DEMproperties->getTaskShowInterval();
    int showInterval = DEMproperties->getShowInterval();

    if (iter_num < 0)
    {
        // generate particles
        int generate_number = 0;
        while (!DEMproperties->isGenerateComplete())
        {
            // DEMproperties->applyExternalForces();
            DEMproperties->handleCollisions();
            DEMproperties->bondforce();
            DEMproperties->motion();
            if (generate_number % showInterval == 0)
            {
                DEMproperties->generateRemainingParticles();
                vis->Update();
            }
            generate_number++;
        }
        // reach a quasi-static state
       
   
        
        generate_number = 0;
        while (DEMproperties->getAverageVelocity() > criticalSpeed || generate_number < 100000)
        {
            // DEMproperties->applyExternalForces();
            DEMproperties->handleCollisions();
            DEMproperties->bondforce();
            DEMproperties->motion();
            if (generate_number % showInterval == 0)
            {
                vis->Update();
                std::cout << " average velocity is " << DEMproperties->getAverageVelocity() << std::endl;
            
            }
            generate_number++;
        }
    }

    // task such as compression, damper......
    DEMproperties->initial_task();
    timeStep = DEMproperties->getTimestep();
    std::string file_folder = "DEMProperties";
    std::filesystem::path dirPath(file_folder); // Directory path

    // Check if the directory exists
    if (!std::filesystem::exists(dirPath))
    {
        // Create the directory if it does not exist
        std::filesystem::create_directories(dirPath);
    }
    while (currentTime < totalTime)
    {
        DEMproperties->applyExternalForces();
        DEMproperties->handleCollisions();
        DEMproperties->bondforce();

        DEMproperties->dotask();

        // Update visualization at specified intervals
        if (iter_num % taskShowInterval == 0)
        {
            vis->Update();
            std::string filename = file_folder + "/DEMProperties" + std::to_string(iter_num) + ".dat";
            DEMproperties->saveToFile(filename);
        }
        currentTime += timeStep;
        DEMproperties->setCurrentTime(currentTime);
        iter_num++;
        DEMproperties->motion();
    }
}
