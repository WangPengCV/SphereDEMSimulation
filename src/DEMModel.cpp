#include "DEMModel.h"
<<<<<<< HEAD
#include <filesystem>
=======
>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200

// Initialize DEMModel
DEMModel::DEMModel(const std::string &filename)
{
    DEMproperties = std::make_shared<DEMProperties>();
    DEMproperties->loadFromFile(filename);
    vis = std::make_shared<Visualization>(*DEMproperties);


    
}

void DEMModel::runSimulation()
{
<<<<<<< HEAD
    

    double currentTime = 0.0;
    DEMproperties->initialSimulation();
    double timeStep = DEMproperties->getTimestep();
    double totalTime = DEMproperties->getTotalTime();
    int showInterval = DEMproperties->getShowInterval();
    //generate particles
    int generate_number = 0;
    while (!DEMproperties->isGenerateComplete())
    {
        //DEMproperties->applyExternalForces();
        DEMproperties->handleCollisions();
        DEMproperties->motion();
        if(generate_number % 50000 == 0)
        {
            DEMproperties->generateRemainingParticles();
            vis->Update();
        }
        generate_number++;
    }
    //reach a quasi-static state
    generate_number =0;
    while(DEMproperties->getAverageVelocity() > 0.001 || generate_number < 50000)
    {
        //DEMproperties->applyExternalForces();
        DEMproperties->handleCollisions();
        DEMproperties->motion();
        if(generate_number % 50000 == 0)
        {
            vis->Update();
            std::cout << " average velocity is " << DEMproperties->getAverageVelocity() << std::endl;
        }
        generate_number++;
    }

    // task such as compression, damper......
    DEMproperties->initial_task();
    int iter_num = 0;
    std::string file_folder = "DEMProperties";
    std::filesystem::path dirPath(file_folder); // Directory path

    // Check if the directory exists
    if (!std::filesystem::exists(dirPath)) {
        // Create the directory if it does not exist
        std::filesystem::create_directories(dirPath);
    }
=======
    double timeStep = DEMproperties->getTimestep();
    double totalTime = DEMproperties->getTotalTime();
    double currentTime = 0.0;
    DEMproperties->initialContactDetection();

    

>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200
    while (currentTime < totalTime)
    {
        DEMproperties->applyExternalForces();
        DEMproperties->handleCollisions();
<<<<<<< HEAD
        DEMproperties->dotask();

        // Update visualization at specified intervals
        if (iter_num % showInterval == 0) {
            vis->Update();
            std::string filename = file_folder + "/DEMProperties" + std::to_string(iter_num) + ".dat";
            DEMproperties->saveToFile(filename);

        }
        currentTime += timeStep;
        iter_num++;
        DEMproperties->motion();



=======
        DEMproperties->motion();
        vis->Update();
        currentTime += timeStep;
>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200
    }
}




