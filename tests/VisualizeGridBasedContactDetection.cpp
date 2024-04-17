#include "DEMProperties.h"
#include "GridBasedContactDetection.h"
#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>
#include "Visualization.h"
#include <vtkSmartPointer.h>
#include <vtkLineSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSphereSource.h>
#include <vtkPlaneSource.h>
#include <vtkProperty.h>
#include <vtkNamedColors.h>
// Function to create and add a line to the renderer
void createLine(vtkSmartPointer<vtkRenderer> &renderer, double x1, double y1, double z1, double x2, double y2, double z2)
{
        vtkSmartPointer<vtkLineSource> lineSource = vtkSmartPointer<vtkLineSource>::New();
        lineSource->SetPoint1(x1, y1, z1);
        lineSource->SetPoint2(x2, y2, z2);

        vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper->SetInputConnection(lineSource->GetOutputPort());

        vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
        actor->SetMapper(mapper);
        renderer->AddActor(actor);
}
// Function to create a 3D grid
void create3DGrid(vtkSmartPointer<vtkRenderer> &renderer, int numX, int numY, int numZ, double gridSizeX, double gridSizeY, double gridSizeZ)
{
        // Draw lines along X-axis
        for (int y = 0; y <= numY; ++y)
        {
                for (int z = 0; z <= numZ; ++z)
                {
                        createLine(renderer, 0, y * gridSizeY, z * gridSizeZ, numX * gridSizeX, y * gridSizeY, z * gridSizeZ);
                }
        }

        // Draw lines along Y-axis
        for (int x = 0; x <= numX; ++x)
        {
                for (int z = 0; z <= numZ; ++z)
                {
                        createLine(renderer, x * gridSizeX, 0, z * gridSizeZ, x * gridSizeX, numY * gridSizeY, z * gridSizeZ);
                }
        }

        // Draw lines along Z-axis
        for (int x = 0; x <= numX; ++x)
        {
                for (int y = 0; y <= numY; ++y)
                {
                        createLine(renderer, x * gridSizeX, y * gridSizeY, 0, x * gridSizeX, y * gridSizeY, numZ * gridSizeZ);
                }
        }
}
void addParticlesToRenderer(const std::unordered_map<int,std::unique_ptr<SphereParticle>> &particles,
                            vtkSmartPointer<vtkRenderer> &renderer)
{
        for (const auto &particle : particles)
        {
                // Assuming you have a method to get the radius and position of the particle
                double radius = 0.01;
                auto position = particle.second->getPosition();

                vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
                sphereSource->SetPhiResolution(10);   // Increase vertical resolution
                sphereSource->SetThetaResolution(10); // Increase horizontal resolution
                sphereSource->SetCenter(position.x(), position.y(), position.z());
                sphereSource->SetRadius(radius);
                sphereSource->Update();

                vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
                mapper->SetInputConnection(sphereSource->GetOutputPort());

                vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
                actor->SetMapper(mapper);
                renderer->AddActor(actor);
        }
}

void visualizeContactPairs(const std::unordered_map<int, std::unordered_set<int>> &contact_pairs,
                           const std::unordered_map<int,std::unique_ptr<SphereParticle>> &particles,
                           vtkSmartPointer<vtkRenderer> &renderer)
{
        vtkSmartPointer<vtkNamedColors> colors = vtkSmartPointer<vtkNamedColors>::New();

        for (const auto &pair : contact_pairs)
        {
                int id1 = pair.first;
                // Create a color for this pair
                double r = static_cast<double>(std::rand()) / RAND_MAX;
                double g = static_cast<double>(std::rand()) / RAND_MAX;
                double b = static_cast<double>(std::rand()) / RAND_MAX;
                colors->SetColor("RandomColor", r, g, b, 1.0);
                for (int id2 : pair.second)
                {
                        // Get positions of the particles
                        Eigen::Vector3d pos1 = particles.at(id1)->getPosition();
                        Eigen::Vector3d pos2 = particles.at(id2)->getPosition();

                        // Create a line between these points
                        vtkSmartPointer<vtkLineSource> lineSource = vtkSmartPointer<vtkLineSource>::New();
                        lineSource->SetPoint1(pos1.x(), pos1.y(), pos1.z());
                        lineSource->SetPoint2(pos2.x(), pos2.y(), pos2.z());

                        // Create a mapper and actor for the line
                        vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
                        mapper->SetInputConnection(lineSource->GetOutputPort());
                        vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
                        actor->SetMapper(mapper);
                        actor->GetProperty()->SetColor(colors->GetColor3d("RandomColor").GetData());
                        actor->GetProperty()->SetLineWidth(5); // Set the line width

                        // Add the actor to the renderer
                        renderer->AddActor(actor);
                }
                // break;
        }
}

void visualizeWallContactPairs(const std::unordered_map<int, std::unordered_set<int>> &contact_pairs,
                               const std::unordered_map<int,std::unique_ptr<SphereParticle>> &particles,
                               const std::unordered_map<int,std::unique_ptr<PlaneWall>> &planewalls,
                               vtkSmartPointer<vtkRenderer> &renderer)
{
        vtkSmartPointer<vtkNamedColors> colors = vtkSmartPointer<vtkNamedColors>::New();

        for (const auto &pair : contact_pairs)
        {
                int id1 = pair.first;
                // Create a color for this pair
                double r = static_cast<double>(std::rand()) / RAND_MAX;
                double g = static_cast<double>(std::rand()) / RAND_MAX;
                double b = static_cast<double>(std::rand()) / RAND_MAX;
                colors->SetColor("RandomColor", r, g, b, 1.0);
                for (int id2 : pair.second)
                {

                        // Assuming you have a method to get the radius and position of the particle
                        double radius = 0.01;
                        auto position = particles.at(id2)->getPosition();

                        vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
                        sphereSource->SetPhiResolution(10);   // Increase vertical resolution
                        sphereSource->SetThetaResolution(10); // Increase horizontal resolution
                        sphereSource->SetCenter(position.x(), position.y(), position.z());
                        sphereSource->SetRadius(radius);
                        sphereSource->Update();
                        // Create a mapper and actor for the line
                        vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
                        mapper->SetInputConnection(sphereSource->GetOutputPort());
                        vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
                        actor->SetMapper(mapper);
                        actor->GetProperty()->SetColor(colors->GetColor3d("RandomColor").GetData());

                        // Add the actor to the renderer
                        renderer->AddActor(actor);
                }
                break;
        }
}
void addPlaneWallsToRenderer(const std::unordered_map<int, std::unique_ptr<PlaneWall>> &planeWalls, vtkSmartPointer<vtkRenderer> &renderer)
{
        for (const auto &wall : planeWalls)
        {
                // Assuming you have methods to get the plane's points and orientation
                auto point1 = wall.second->getCorner1();
                auto point2 = wall.second->getCorner2();
                auto point3 = wall.second->getCorner3();

                vtkSmartPointer<vtkPlaneSource> planeSource = vtkSmartPointer<vtkPlaneSource>::New();
                planeSource->SetOrigin(point2.x(), point2.y(), point2.z());
                planeSource->SetPoint1(point1.x(), point1.y(), point1.z());
                planeSource->SetPoint2(point3.x(), point3.y(), point3.z());
                planeSource->Update();

                vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
                mapper->SetInputConnection(planeSource->GetOutputPort());

                vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
                actor->SetMapper(mapper);
                // Set properties like color, transparency, etc.
                actor->GetProperty()->SetColor(1.0, 1.0, 1.0); // Red color
                actor->GetProperty()->EdgeVisibilityOn();
                actor->GetProperty()->SetOpacity(0.1);

                renderer->AddActor(actor);
        }
}
int main(int argc, char **argv)
{
        std::string inputfilename = "InputFile.txt";
        DEMProperties demproperties;

        demproperties.loadFromFile(inputfilename);

        double gird_size = 0.02;
        Eigen::Vector3d simulationdimensions = demproperties.getSimulationDimensions();
        auto gridbasedcontactdetection = std::make_shared<GridBasedContactDetection>();
        gridbasedcontactdetection->initial(simulationdimensions.x(), simulationdimensions.y(),simulationdimensions.z(), gird_size);

        // Assuming you have these values from GridBasedContactDetection
        double gridSizeX = gridbasedcontactdetection->getGridSizeX();
        double gridSizeY = gridbasedcontactdetection->getGridSizeY();
        double gridSizeZ = gridbasedcontactdetection->getGridSizeZ();
        int numberOfGridX = gridbasedcontactdetection->getNumberOfGridX();
        int numberOfGridY = gridbasedcontactdetection->getNumberOfGridY();
        int numberOfGridZ = gridbasedcontactdetection->getNumberOfGridZ();

        vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
        create3DGrid(renderer, numberOfGridX, numberOfGridY, numberOfGridZ, gridSizeX, gridSizeY, gridSizeZ);

        auto& planewalls = demproperties.getPlaneWall();

        auto& sphereparticles = demproperties.getsphereParticles();
        std::unordered_map<int, std::unordered_set<int>> pp_contact_paris;
        std::unordered_map<int, std::unordered_set<int>> pw_contact_paris;

        // gridbasedcontactdetection->ParticleBroadPhase(sphereparticles, pp_contact_paris);
        // gridbasedcontactdetection->planewallBroadPhase(sphereparticles, planewalls, pw_contact_paris);

        // Add particles and plane walls to the renderer
        addParticlesToRenderer(sphereparticles, renderer);
        addPlaneWallsToRenderer(planewalls, renderer);

        //visualizeContactPairs(pp_contact_paris, particles, renderer);
        visualizeWallContactPairs(pw_contact_paris,sphereparticles,planewalls,renderer);
        vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
        renderWindow->AddRenderer(renderer);
        vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
        renderWindowInteractor->SetRenderWindow(renderWindow);

        // Start the interaction
        renderWindow->Render();
        renderWindowInteractor->Start();
}