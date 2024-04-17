#pragma once

#include <vtkSmartPointer.h>
#include "DEMProperties.h"
#include "LineCylinderSource.h"
class vtkPoints;
class vtkRenderer;
class vtkRenderWindow;
class vtkRenderWindowInteractor;
class vtkPolyData;
class vtkSphereSource;
class vtkCylinderSource;
class vtkActor;
class vtkPlaneSource;
class vtkCellArray;
class vtkGlyph3DMapper;

class Visualization
{
public:
    explicit Visualization(const DEMProperties &DEMproperties);
    void Update();

private:
    const DEMProperties &DEMproperties;

    // sphere
    vtkSmartPointer<vtkPoints> spherepoints;
    vtkSmartPointer<vtkPolyData> spherepolyData;

    // cylinder
    vtkSmartPointer<vtkPoints> cylinderpoints;
    vtkSmartPointer<vtkPolyData> cylinderpolydata;
    //fibersphere
    vtkSmartPointer<vtkPoints> fiberspherepoints;
    vtkSmartPointer<vtkPolyData> fiberspherepolyData;

    vtkSmartPointer<vtkRenderer> renderer;
    vtkSmartPointer<vtkRenderWindow> renderWindow;
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor;

    // Plane wall
    std::vector<vtkSmartPointer<vtkPlaneSource>> planeWallSource;
    std::vector<vtkSmartPointer<vtkActor>> planeWallActor;

    //cylinder wall
    std::vector<vtkSmartPointer<LineCylinderSource>> cylinderWallSource;
    std::vector<vtkSmartPointer<vtkActor>> cylinderWallActor;

    void UpdateCylinderWall();
    void UpdatePlaneWall();
    void UpdateSphere();
    void UpdataFibers();

    int count_numer;
};
