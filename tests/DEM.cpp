#include "DEMModel.h"

int main()
{
    std::string inputfilename = "InputFile.txt";
    DEMModel dem(inputfilename);
    dem.runSimulation();

    return 0;
}