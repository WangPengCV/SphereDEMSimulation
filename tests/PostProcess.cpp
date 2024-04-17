#include "DEMProperties.h"
#include <filesystem>
#include <regex>

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

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " <stateFolder>" << std::endl;
        return 1;
    }

    std::string stateFolder = argv[1];
    if (!canResumeFromState(stateFolder))
    {
        std::cerr << "Error: Invalid or empty state folder" << std::endl;
        return 1;
    }
    std::vector<long long> files; // changed to long long to match extracted numbers
    if (canResumeFromState(stateFolder))
    {
        std::filesystem::path dirPath(stateFolder);

        for (const auto &entry : std::filesystem::directory_iterator(dirPath))
        {
            const auto &path = entry.path();
            if (entry.is_regular_file())
            {
                int number = extractNumberFromFilename(path.filename().string()); // pass filename instead of path
                files.push_back(number);
            }
        }

        std::sort(files.begin(), files.end()); // sort extracted numbers
    }
    std::string filenamepre = stateFolder + "/DEMProperties/DEMProperties";

    for (const auto &file : files)
    {
        std::string filename = filenamepre + std::to_string(file) + ".dat";
        auto DEMproperties = std::make_shared<DEMProperties>();
        DEMproperties->loadFromFile(filename); // Make sure this function loads all needed state
        auto &planewallSphereContactInformationList = DEMproperties->getContactForce()->getplanewallSphereContactInformationList();
        auto &spheres = DEMproperties->getsphereParticles();








        // const auto &fibersphere = DEMproperties->getfibersphereParticles();
        // for (const auto &bond : DEMproperties->getfiberbonds())
        // {
        //     const auto &sphere1 = fibersphere.at(bond.second->getNode1());
        //     const auto &sphere2 = fibersphere.at(bond.second->getNode2());

        //     potentialenergy += bond.second->getPotentialenergy(sphere1, sphere2);
        // }

        // const auto &planewallSphereContactInformationList = DEMproperties->getContactForce()->getplanewallSphereContactInformationList();
        // const auto &planewallFiberContactInformationList = DEMproperties->getContactForce()->getplanewallFiberContactInformationList();

        // const auto &cylinderwallFiberContactInformationList = DEMproperties->getContactForce()->getcylinderwallFiberContactInformationList();
        // const auto &cylinderwallSphereContactInformationList = DEMproperties->getContactForce()->getcylinderwallSphereContactInformationList();

        // const auto &sphereSphereContactInformationList = DEMproperties->getContactForce()->getSphereSphereContactInformationList();
        // const auto &sphereFiberContactInformationList = DEMproperties->getContactForce()->getSphereFiberContactInformationList();
        // const auto &fiberFiberContactInformationList = DEMproperties->getContactForce()->getFiberFiberContactInformationList();

        // double contact_number = 0;
        // for (const auto &contactforceLists : planewallSphereContactInformationList)
        // {
        //     contact_number += contactforceLists.second.size();
        // }
        // for (const auto &contactforceLists : planewallFiberContactInformationList)
        // {
        //     contact_number += contactforceLists.second.size();
        // }

        // for (const auto &contactforceLists : cylinderwallFiberContactInformationList)
        // {
        //     contact_number += contactforceLists.second.size();
        // }

        // for (const auto &contactforceLists : cylinderwallSphereContactInformationList)
        // {
        //     contact_number += contactforceLists.second.size();
        // }

        // for (const auto &contactforceLists : cylinderwallSphereContactInformationList)
        // {
        //     contact_number += contactforceLists.second.size();
        // }

        // std::ofstream file(stateFolder + "/postprocess.txt", std::ios::app);
        // file.precision(std::numeric_limits<double>::digits10 + 1);
        // file << DEMproperties->getCurrentTime() << " " << potentialenergy << " " << contact_number << std::endl;
        // file.close();
    }
    return 0;
}