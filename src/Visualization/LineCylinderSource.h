#pragma once
#include <vtkPolyDataAlgorithm.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkCylinderSource.h>
#include <vtkMath.h>
#include <vtkObjectFactory.h>
#include <vtkTransform.h>
#include <vtkInformationVector.h>
#include <vtkInformation.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkMinimalStandardRandomSequence.h>

class LineCylinderSource : public vtkPolyDataAlgorithm {
public:
    vtkTypeMacro(LineCylinderSource, vtkPolyDataAlgorithm);
    static LineCylinderSource* New();

    vtkSetMacro(Radius, double);
    vtkGetMacro(Radius, double);

    vtkSetVector3Macro(StartPoint, double);
    vtkGetVectorMacro(StartPoint, double, 3);

    vtkSetVector3Macro(EndPoint, double);
    vtkGetVectorMacro(EndPoint, double, 3);

protected:
    LineCylinderSource();
    ~LineCylinderSource() override;

    int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

private:
    double Radius;
    double StartPoint[3]; // 圆柱起点
    double EndPoint[3];   // 圆柱终点

    void CreateCylinder(vtkPoints* points, vtkCellArray* polys);
};

