#include <vtkActor.h>
#include <vtkCubeSource.h>
#include <vtkNamedColors.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>
#include <vtkProperty.h>

#include <array>
#include <cmath>
#include <DEMProperties.h>
int main()
{
    // 设置渲染器、渲染窗口和交互器
    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);
    std::string inputfilename = "InputFile.txt";

    std::shared_ptr<DEMProperties> DEMproperties = std::make_shared<DEMProperties>();
    DEMproperties->loadFromFile(inputfilename);
    DEMproperties->initialSimulation();

    int numberOfGridX = DEMproperties->getGridBasedContactDetection()->getNumberOfGridX();
    int numberOfGridY = DEMproperties->getGridBasedContactDetection()->getNumberOfGridY();
    int numberOfGridZ = DEMproperties->getGridBasedContactDetection()->getNumberOfGridZ();

    double gridSizeX = DEMproperties->getGridBasedContactDetection()->getGridSizeX();
    double gridSizeY = DEMproperties->getGridBasedContactDetection()->getGridSizeY();
    double gridSizeZ = DEMproperties->getGridBasedContactDetection()->getGridSizeZ();

    std::vector<std::vector<int>> sphereGrids;
    std::vector<std::vector<int>> fiberGrids;
    int gridNumber = numberOfGridX * numberOfGridZ * numberOfGridY;
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
    int world_size = 9;  
    std::vector<int> offset; 
    offset.resize(2*world_size);
    std::vector<int> offsetX; 
    offsetX.resize(2*world_size);
    std::vector<int> offsetY; 
    offsetY.resize(2*world_size);
    std::vector<int> offsetZ; 
    offsetZ.resize(2*world_size);
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
                count_number += sphereGrids[gridIndex].size();
                count_number += fiberGrids[gridIndex].size();

                if (0.99 * averageParticleNumber <= count_number || gridIndex == gridNumber - 1)
                {
                    if (rank_index == world_size - 1)
                    {
                        i = numberOfGridX;
                        j = numberOfGridZ;
                        k = numberOfGridY;
                        offset[rank_index * 2 + 1] = gridNumber-1;
                        offsetX[rank_index * 2 + 1] = i-1;
                        offsetY[rank_index * 2 + 1] = k-1;
                        offsetZ[rank_index * 2 + 1] = j-1;
                    }
                    else
                    {
                        offset[rank_index * 2 + 1] = gridIndex;
                        offsetX[rank_index * 2 + 1] = i;
                        offsetY[rank_index * 2 + 1] = k;
                        offsetZ[rank_index * 2 + 1] = j;

                        count_number = 0;
                        rank_index++;
                        offset[rank_index * 2] = gridIndex + 1;
                        // Correctly handle boundary conditions when setting the start of the next segment
                        if (i + 1 >= numberOfGridX)
                        {
                            offsetX[rank_index * 2] = 0; // Reset i to the beginning of the next row
                            if (k + 1 >= numberOfGridY)
                            {
                                offsetY[rank_index * 2] = 0;     // Reset k to the beginning of the next "layer"
                                offsetZ[rank_index * 2] = j + 1; // Move to the next "layer" along Z-axis
                            }
                            else
                            {
                                offsetY[rank_index * 2] = k + 1; // Move to the next row
                                offsetZ[rank_index * 2] = j;     // Stay in the same "layer" along Z-axis
                            }
                        }
                        else
                        {
                            offsetX[rank_index * 2] = i + 1; // Move to the next column
                            offsetY[rank_index * 2] = k;     // Stay in the same row
                            offsetZ[rank_index * 2] = j;     // Stay in the same "layer" along Z-axis
                        }
                    }
                }
            }
        }
    }                                  // 总进程数
    std::vector<std::unordered_map<int, std::unordered_set<int>>> ghostLayer;
    ghostLayer.resize(world_size);
    std::vector<std::unordered_set<int>> neighborRankId;
    neighborRankId.resize(world_size);
    for (int currentRank = 0; currentRank < world_size; ++currentRank)
    {

        for (int gridID = offset[currentRank * 2]; gridID <= offset[currentRank * 2 + 1]; ++gridID)
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
                                if (currentRank == otherRank)
                                    continue; // Skip comparing the same rank with itself.
                                if (neighborId >= offset[otherRank * 2] && neighborId <= offset[otherRank * 2 + 1])
                                {
                                    neighborRankId[currentRank].insert(otherRank);
                                    ghostLayer[currentRank][otherRank].insert(gridID);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    vtkSmartPointer<vtkNamedColors> colors = vtkSmartPointer<vtkNamedColors>::New();

    // 计算每个立方体的边长
    double cubeSize = 1.0;
    // 预定义颜色数组
    std::vector<std::string> colorNames = {
        "Tomato", "DodgerBlue", "LimeGreen", "Gold",
        "Violet", "SlateGray", "LightSeaGreen", "Chocolate",
        "CadetBlue", "OrangeRed"};
    // 为每个进程管理的网格创建立方体
    for (int process = 0; process < world_size; ++process)
    {
        unsigned char color[3];
        colors->GetColor(colorNames[process].c_str(), color);

        for (int id = offset[process * 2]; id <= offset[process * 2 + 1]; ++id)
        {
            // 计算立方体中心的位置
            int y = id / (numberOfGridX * numberOfGridZ);
            int temp = id % (numberOfGridX * numberOfGridZ);
            int z = temp / numberOfGridX;
            int x = temp % numberOfGridX;

            vtkSmartPointer<vtkCubeSource> cubeSource = vtkSmartPointer<vtkCubeSource>::New();
            cubeSource->SetCenter(x * cubeSize, y * cubeSize, z * cubeSize);
            cubeSource->SetXLength(cubeSize);
            cubeSource->SetYLength(cubeSize);
            cubeSource->SetZLength(cubeSize);

            vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
            mapper->SetInputConnection(cubeSource->GetOutputPort());

            vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
            actor->SetMapper(mapper);
            actor->GetProperty()->SetColor(color[0] / 255.0, color[1] / 255.0, color[2] / 255.0);
            actor->GetProperty()->EdgeVisibilityOn();
            actor->GetProperty()->SetEdgeColor(0.0, 0.0, 0.0); // 设置边线为黑色
            actor->GetProperty()->SetLineWidth(1);             // 设置边线宽度
            renderer->AddActor(actor);
        }
    }

    renderer->SetBackground(colors->GetColor3d("White").GetData());
    renderWindow->SetSize(800, 600);
    renderWindow->Render();
    renderWindowInteractor->Start();

    return 0;
}
