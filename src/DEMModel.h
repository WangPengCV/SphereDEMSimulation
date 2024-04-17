#pragma once
#include "Particle.h"
#include "SphereParticle.h"
#include "DEMProperties.h"
#include "Visualization.h"
#include <vector>
#include <memory>
#include <filesystem>
#include <regex>

class DEMModel
{
public:
    DEMModel(const std::string &filename);

    void runSimulation();

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

private:
    std::shared_ptr<DEMProperties> DEMproperties;
    std::shared_ptr<Visualization> vis;
    int iter_num;
};
