#include "ParallelDEMModel.h"


int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);
    std::string inputfilename = "InputFile.txt";
    ParallelDEMModel dem(inputfilename);

    dem.runSimulation();
    dem.finalizeMPIType();
    // Ensure to finalize the MPI environment properly before exiting
    MPI_Finalize();
    return 0;
}