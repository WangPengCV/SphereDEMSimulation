// Visualization.cpp
#include "Visualization.h"
#include <vtkCylinderSource.h>
#include <vtkSphereSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkGlyph3D.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyData.h>
#include <vtkProperty.h>
#include <vtkCubeSource.h>
#include <vtkFloatArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkPointData.h>
#include <vtkLight.h>
#include <vtkCamera.h>
#include <vtkPlaneSource.h>
#include <vtkAppendPolyData.h>
#include <vtkRegularPolygonSource.h>
#include "LineSphereSource.h"
#include <vtkWindowToImageFilter.h>
#include <vtkBMPWriter.h>
#include <vtkTubeFilter.h>
#include <vtkCellArray.h>
#include <vtkGlyph3DMapper.h>
#include <filesystem>

Visualization::Visualization(const DEMProperties &DEMproperties) : DEMproperties(DEMproperties)
{
    count_numer = 0;
    std::filesystem::path dirPath("Image"); // Directory path

    // Check if the directory exists
    if (!std::filesystem::exists(dirPath))
    {
        // Create the directory if it does not exist
        std::filesystem::create_directories(dirPath);
    }
    spherepoints = vtkSmartPointer<vtkPoints>::New();
    fiberspherepoints = vtkSmartPointer<vtkPoints>::New();
    cylinderpoints = vtkSmartPointer<vtkPoints>::New();

    UpdateSphere();

    vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
    sphereSource->SetPhiResolution(20);   // Increase vertical resolution
    sphereSource->SetThetaResolution(20); // Increase horizontal resolution
    sphereSource->SetRadius(1);

    vtkSmartPointer<vtkGlyph3D> glyph3D = vtkSmartPointer<vtkGlyph3D>::New();
    glyph3D->SetSourceConnection(sphereSource->GetOutputPort());
    glyph3D->SetInputData(spherepolyData);
    glyph3D->SetScaleModeToScaleByScalar();
    glyph3D->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "scales");
    glyph3D->SetColorModeToColorByScalar();
    glyph3D->SetInputArrayToProcess(1, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "colors");

    glyph3D->Update();

    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(glyph3D->GetOutputPort());
    mapper->SetColorModeToDefault();
    mapper->SetScalarVisibility(true);

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);

    renderer = vtkSmartPointer<vtkRenderer>::New();
    renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    renderer->AddActor(actor);

    for (const auto &wall : DEMproperties.getPlaneWall())
    {
        // Create a plane source
        vtkSmartPointer<vtkPlaneSource> planesource = vtkSmartPointer<vtkPlaneSource>::New();
        planesource->SetOrigin(wall.second->getCorner2().x(), wall.second->getCorner2().y(), wall.second->getCorner2().z());
        planesource->SetPoint1(wall.second->getCorner1().x(), wall.second->getCorner1().y(), wall.second->getCorner1().z());
        planesource->SetPoint2(wall.second->getCorner3().x(), wall.second->getCorner3().y(), wall.second->getCorner3().z());

        // Map the transformed plane to a mapper
        vtkSmartPointer<vtkPolyDataMapper> planemapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        planemapper->SetInputConnection(planesource->GetOutputPort());

        // Create an actor to represent the wall
        vtkSmartPointer<vtkActor> planeactor = vtkSmartPointer<vtkActor>::New();
        planeactor->SetMapper(planemapper);
        // Set properties like color, transparency, etc.
        planeactor->GetProperty()->SetColor(1.0, 1.0, 1.0); // Red color
        planeactor->GetProperty()->EdgeVisibilityOn();
        planeactor->GetProperty()->SetOpacity(0.1);

        // Add the actor to the renderer
        renderer->AddActor(planeactor);

        // Store the sources and actors for later updates
        planeWallSource.push_back(planesource);
        planeWallActor.push_back(planeactor);
    }
    

    for (const auto &wall : DEMproperties.getCylinderWall())
    {

        double radius = wall.second->getRadius();

        // Calculate cylinder direction vector
        Eigen::Vector3d startpoint = wall.second->getStartpoint();
        Eigen::Vector3d endpoint = wall.second->getEndpoint();

        vtkSmartPointer<LineCylinderSource> cylinderSource = vtkSmartPointer<LineCylinderSource>::New();
        // 设置圆柱的高度和半径
        cylinderSource->SetRadius(radius);
        cylinderSource->SetStartPoint(startpoint.x(), startpoint.y(), startpoint.z());
        cylinderSource->SetEndPoint(endpoint.x(), endpoint.y(), endpoint.z());
        cylinderSource->Update();

        // 创建圆柱的 Mapper
        vtkSmartPointer<vtkPolyDataMapper> cylinderMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        cylinderMapper->SetInputConnection(cylinderSource->GetOutputPort());

        // 创建圆柱的 Actor
        vtkSmartPointer<vtkActor> cylinderActor = vtkSmartPointer<vtkActor>::New();
        cylinderActor->SetMapper(cylinderMapper);
        vtkSmartPointer<vtkProperty> cylinderProperty = vtkSmartPointer<vtkProperty>::New();
        cylinderProperty->SetOpacity(0.5); // 设置不透明度为0.5，即半透明
        cylinderActor->SetProperty(cylinderProperty);
        // 将圆柱的 Actor 添加到渲染器
        renderer->AddActor(cylinderActor);

        // 存储圆柱的 Source 和 Actor
        cylinderWallSource.push_back(cylinderSource);
        cylinderWallActor.push_back(cylinderActor);
    }

    UpdataFibers();
    vtkSmartPointer<vtkCylinderSource> cylinder = vtkSmartPointer<vtkCylinderSource>::New();
    cylinder->SetResolution(20);
    cylinder->SetHeight(1.0); // Default height, will be scaled
    cylinder->SetRadius(1.0); // Default radius, will be scaled
    cylinder->Update();
    vtkSmartPointer<vtkGlyph3DMapper> cylinderglyph3D = vtkSmartPointer<vtkGlyph3DMapper>::New();

    cylinderglyph3D->SetSourceConnection(cylinder->GetOutputPort());
    cylinderglyph3D->SetInputData(cylinderpolydata);
    cylinderglyph3D->SetScalarModeToUsePointFieldData();
    cylinderglyph3D->SetScaleArray("Scales");
    cylinderglyph3D->SetScaleModeToScaleByVectorComponents();
    cylinderglyph3D->SelectColorArray("Colors");
    cylinderglyph3D->SetOrientationArray("Orientation");
    cylinderglyph3D->SetOrientationModeToQuaternion();
    cylinderglyph3D->Update();

    vtkSmartPointer<vtkActor> cylinderactor = vtkSmartPointer<vtkActor>::New();
    cylinderactor->SetMapper(cylinderglyph3D);
    renderer->AddActor(cylinderactor);

    vtkSmartPointer<vtkGlyph3D> fibersphereglyph3D = vtkSmartPointer<vtkGlyph3D>::New();
    fibersphereglyph3D->SetSourceConnection(sphereSource->GetOutputPort());
    fibersphereglyph3D->SetInputData(fiberspherepolyData);
    fibersphereglyph3D->SetScaleModeToScaleByScalar();
    fibersphereglyph3D->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "scales");
    fibersphereglyph3D->SetColorModeToColorByScalar();
    fibersphereglyph3D->SetInputArrayToProcess(1, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "colors");
    fibersphereglyph3D->Update();

    vtkSmartPointer<vtkPolyDataMapper> fiberspheremapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    fiberspheremapper->SetInputConnection(fibersphereglyph3D->GetOutputPort());
    fiberspheremapper->SetColorModeToDefault();
    fiberspheremapper->SetScalarVisibility(true);

    vtkSmartPointer<vtkActor> fibersphereactor = vtkSmartPointer<vtkActor>::New();
    fibersphereactor->SetMapper(fiberspheremapper);
    renderer->AddActor(fibersphereactor);

    renderer->SetBackground(1, 1, 1);
    renderWindow->SetSize(540, 540);
    renderWindow->SetPosition(0, 0);
    renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);
}

void Visualization::Update()
{
    UpdateSphere();
    UpdatePlaneWall();
    UpdateCylinderWall();
    UpdataFibers();
    renderWindowInteractor->GetRenderWindow()->Render();
    if (renderWindowInteractor)
    {
        renderWindowInteractor->ProcessEvents();
    }
    vtkSmartPointer<vtkWindowToImageFilter> windowto_image_filter = vtkSmartPointer<vtkWindowToImageFilter>::New();
    windowto_image_filter->SetInput(renderWindow);
    windowto_image_filter->SetScale(1);
    windowto_image_filter->SetInputBufferTypeToRGB();
    windowto_image_filter->ReadFrontBufferOff();
    windowto_image_filter->Update();
    vtkSmartPointer<vtkBMPWriter> writer = vtkSmartPointer<vtkBMPWriter>::New();

    std::string bmpname = "Image/out" + std::to_string(count_numer) + ".jpg";

    const char *cfirst = bmpname.c_str();
    writer->SetFileName(cfirst);
    writer->SetInputConnection(windowto_image_filter->GetOutputPort());
    writer->Write();
    count_numer++;
}

void Visualization::UpdateSphere()
{
    spherepoints->Reset();

    // Create scale and color arrays
    vtkSmartPointer<vtkFloatArray> scales = vtkSmartPointer<vtkFloatArray>::New();
    scales->SetName("scales");
    vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
    colors->SetName("colors");
    colors->SetNumberOfComponents(3); // 3 components for R, G, B
    for (const auto &particle : DEMproperties.getsphereParticles())
    {
        auto manager = particle.second->getParticlePropertyManager();

        spherepoints->InsertNextPoint(particle.second->getPosition().x(), particle.second->getPosition().y(), particle.second->getPosition().z());
        float scale = manager->getSphereProperties(particle.second->getType())->getRadius();
        scales->InsertNextValue(scale);
        unsigned char color[3] = {0, 0, 255};
        colors->InsertNextTypedTuple(color);
    }
    // spherepoints->Modified();

    if (!spherepolyData)
    {
        spherepolyData = vtkSmartPointer<vtkPolyData>::New();
    }
    spherepolyData->SetPoints(spherepoints);
    spherepolyData->GetPointData()->SetScalars(colors);
    spherepolyData->GetPointData()->AddArray(scales);
}

void Visualization::UpdataFibers()
{
    cylinderpoints->Reset();
    fiberspherepoints->Reset();

    vtkSmartPointer<vtkFloatArray> cylinderscales = vtkSmartPointer<vtkFloatArray>::New();
    cylinderscales->SetNumberOfComponents(3);
    cylinderscales->SetName("Scales");

    vtkSmartPointer<vtkFloatArray> cylinderorientation = vtkSmartPointer<vtkFloatArray>::New();
    cylinderorientation->SetNumberOfComponents(4);
    cylinderorientation->SetName("Orientation");

    vtkSmartPointer<vtkUnsignedCharArray> cylindercolors = vtkSmartPointer<vtkUnsignedCharArray>::New();
    cylindercolors->SetNumberOfComponents(3);
    cylindercolors->SetName("Colors");

    // Create scale and color arrays
    vtkSmartPointer<vtkFloatArray> spherescales = vtkSmartPointer<vtkFloatArray>::New();
    spherescales->SetName("scales");
    vtkSmartPointer<vtkUnsignedCharArray> spherecolors = vtkSmartPointer<vtkUnsignedCharArray>::New();
    spherecolors->SetName("colors");
    spherecolors->SetNumberOfComponents(3); // 3 components for R, G, B

    const auto &nodes = DEMproperties.getfibersphereParticles();
    const auto &manger = DEMproperties.getParticleManager();
    Eigen::Vector3d referenceVector(0, 1, 0); // cylinders are aligned along y-axis initially

    for (const auto &fiberbond : DEMproperties.getfiberbonds())
    {

        const auto &type = fiberbond.second->getType();
        const auto &node1 = nodes.at(fiberbond.second->getNode1())->getPosition();
        const auto &node2 = nodes.at(fiberbond.second->getNode2())->getPosition();
        cylinderpoints->InsertNextPoint(0.5 * (node1.x() + node2.x()),
                                        0.5 * (node1.y() + node2.y()),
                                        0.5 * (node1.z() + node2.z()));

        cylinderscales->InsertNextTuple3(manger->getFiberProperties(type)->getRadius(), manger->getFiberProperties(type)->getElementlength(),
                                         manger->getFiberProperties(type)->getRadius());
        // Calculate the direction vector and normalize
        Eigen::Vector3d dir = (node2 - node1).normalized();

        // Compute rotation axis (cross product) and angle
        Eigen::Vector3d axis = referenceVector.cross(dir);
        double angle = acos(referenceVector.dot(dir)); // In radians

        // Convert angle and axis to a quaternion
        Eigen::Quaterniond quat;
        quat = Eigen::AngleAxisd(angle, axis.normalized());

        // Add quaternion components to the orientation array
        cylinderorientation->InsertNextTuple4(quat.w(), quat.x(), quat.y(), quat.z());

        unsigned char color[3] = {255, 0, 0};
        cylindercolors->InsertNextTypedTuple(color);
        fiberspherepoints->InsertNextPoint(node1.x(), node1.y(), node1.z());
        spherescales->InsertNextValue(manger->getFiberProperties(type)->getRadius());
        spherecolors->InsertNextTypedTuple(color);
        if (fiberbond.second->getNeighborelement2() == -1)
        {
            fiberspherepoints->InsertNextPoint(node2.x(), node2.y(), node2.z());
            spherescales->InsertNextValue(manger->getFiberProperties(type)->getRadius());
            spherecolors->InsertNextTypedTuple(color);
        }
    }
    if (!cylinderpolydata)
    {
        cylinderpolydata = vtkSmartPointer<vtkPolyData>::New();
    }
    cylinderpolydata->SetPoints(cylinderpoints);
    cylinderpolydata->GetPointData()->SetScalars(cylinderscales);
    cylinderpolydata->GetPointData()->AddArray(cylinderorientation);
    cylinderpolydata->GetPointData()->AddArray(cylindercolors);

    if (!fiberspherepolyData)
    {
        fiberspherepolyData = vtkSmartPointer<vtkPolyData>::New();
    }
    fiberspherepolyData->SetPoints(fiberspherepoints);
    fiberspherepolyData->GetPointData()->SetScalars(spherecolors);
    fiberspherepolyData->GetPointData()->AddArray(spherescales);
}

void Visualization::UpdatePlaneWall()
{
    for (const auto &wall : DEMproperties.getPlaneWall())
    {

        planeWallSource[wall.first]->SetOrigin(wall.second->getCorner2().x(), wall.second->getCorner2().y(), wall.second->getCorner2().z());
        planeWallSource[wall.first]->SetPoint1(wall.second->getCorner1().x(), wall.second->getCorner1().y(), wall.second->getCorner1().z());
        planeWallSource[wall.first]->SetPoint2(wall.second->getCorner3().x(), wall.second->getCorner3().y(), wall.second->getCorner3().z());
    }
    
}
void Visualization::UpdateCylinderWall()
{
    for (const auto &wall : DEMproperties.getCylinderWall())
    {

        double radius = wall.second->getRadius();

        // Calculate cylinder direction vector
        Eigen::Vector3d startpoint = wall.second->getStartpoint();
        Eigen::Vector3d endpoint = wall.second->getEndpoint();

        cylinderWallSource[wall.first]->SetStartPoint(startpoint.x(), startpoint.y(), startpoint.z());
        cylinderWallSource[wall.first]->SetEndPoint(endpoint.x(), endpoint.y(), endpoint.z());
    }
}
