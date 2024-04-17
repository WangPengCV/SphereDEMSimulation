#include "LineCylinderSource.h"

vtkStandardNewMacro(LineCylinderSource);

LineCylinderSource::LineCylinderSource() : Radius(1.0)
{
    this->SetNumberOfInputPorts(0);
    this->StartPoint[0] = 0.0; // Default start point
    this->StartPoint[1] = 0.0;
    this->StartPoint[2] = 0.0;
    this->EndPoint[0] = 1.0; // Default end point
    this->EndPoint[1] = 1.0;
    this->EndPoint[2] = 1.0;
}

int LineCylinderSource::RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *outputVector)
{
    // 获取输出信息
    vtkInformation *outInfo = outputVector->GetInformationObject(0);
    vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> polys = vtkSmartPointer<vtkCellArray>::New();

    CreateCylinder(points, polys);

    output->SetPoints(points);
    output->SetPolys(polys);

    return 1;
}

void LineCylinderSource::CreateCylinder(vtkPoints *points, vtkCellArray *polys)
{
    // Compute a basis
    double normalizedX[3];
    double normalizedY[3];
    double normalizedZ[3];
    // The X axis is a vector from start to end
    vtkMath::Subtract(EndPoint, StartPoint, normalizedY);
    double height = vtkMath::Norm(normalizedY);
    vtkMath::Normalize(normalizedY);

   
    double center[3];
    for (int i = 0; i < 3; ++i)
    {
        center[i] = (this->StartPoint[i] + this->EndPoint[i]) / 2.0;
    }
    // 创建圆柱
    vtkSmartPointer<vtkCylinderSource> cylinderSrc = vtkSmartPointer<vtkCylinderSource>::New();
    cylinderSrc->SetRadius(this->Radius);
    cylinderSrc->SetHeight(height);
    cylinderSrc->SetResolution(50);
    cylinderSrc->Update();

    vtkNew<vtkMinimalStandardRandomSequence> rng;
    rng->SetSeed(8775070); // For testing.

    // The Z axis is an arbitrary vector cross Y
    double arbitrary[3];
    for (auto i = 0; i < 3; ++i)
    {
        rng->Next();
        arbitrary[i] = rng->GetRangeValue(-10, 10);
    }
    vtkMath::Cross(normalizedY, arbitrary, normalizedZ);
    vtkMath::Normalize(normalizedZ);

    // The X axis is Z cross Y
    vtkMath::Cross(normalizedZ, normalizedY, normalizedX);
    vtkNew<vtkMatrix4x4> matrix;

    // Create the direction cosine matrix
    matrix->Identity();
    for (unsigned int i = 0; i < 3; i++)
    {
        matrix->SetElement(i, 0, normalizedX[i]);
        matrix->SetElement(i, 1, normalizedY[i]);
        matrix->SetElement(i, 2, normalizedZ[i]);
    }

    // Apply the transforms
    vtkNew<vtkTransform> transform;
    transform->Translate(center); // translate to starting point
    transform->Concatenate(matrix);   // apply direction cosines

    vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    transformFilter->SetInputData(cylinderSrc->GetOutput());
    transformFilter->SetTransform(transform);
    transformFilter->Update();

    // 复制圆柱的点和多边形
    vtkPolyData *cylinderData = transformFilter->GetOutput();
    points->DeepCopy(cylinderData->GetPoints());
    polys->DeepCopy(cylinderData->GetPolys());
}

LineCylinderSource::~LineCylinderSource()
{
    // 析构函数代码
}
