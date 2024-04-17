#include "ContactForce.h"

ContactForce::ContactForce()
{
    energyDissipation = 0;
}

void ContactForce::setEnergydissipation(double energyDissipation)
{
    this->energyDissipation = energyDissipation;
}

void ContactForce::computerSphereSphereEffective(std::shared_ptr<SphereProperties> &sphereProps1, PropertyTypeID &type1,
                                                 std::shared_ptr<SphereProperties> &sphereProps2, PropertyTypeID &type2)
{
    // Handle sphere property
    const double radius1 = sphereProps1->getRadius();
    const double radius2 = sphereProps2->getRadius();

    const double youngs_modulus1 = sphereProps1->getYoungModulus();
    const double youngs_modulus2 = sphereProps2->getYoungModulus();

    const double poisson_ratio1 = sphereProps1->getPoissonRatio();
    const double poisson_ratio2 = sphereProps2->getPoissonRatio();

    const double restitution1 = sphereProps1->getRestitution();
    const double restitution2 = sphereProps2->getRestitution();

    const double contact_damping_coefficient1 = -log(restitution1) / sqrt(PI * PI + log(restitution1) * log(restitution1));
    const double contact_damping_coefficient2 = -log(restitution2) / sqrt(PI * PI + log(restitution2) * log(restitution2));

    const double sliding_friction1 = sphereProps1->getSlidingFriction();
    const double sliding_friction2 = sphereProps2->getSlidingFriction();

    const double rolling_friction1 = sphereProps1->getRollingFriction();
    const double rolling_friction2 = sphereProps2->getRollingFriction();

    const double mass1 = sphereProps1->getMass();
    const double mass2 = sphereProps2->getMass();

    effectiveradius[type1][type2] = 2 * (radius1 * radius2) / (radius1 + radius2);

    effectivemass[type1][type2] = 2 * (mass1 * mass2) / (mass1 + mass2);

    effectiveyoungsmodulus[type1][type2] = (youngs_modulus1 * youngs_modulus2) /
                                           ((youngs_modulus2 * (1 - poisson_ratio1 * poisson_ratio1)) +
                                            (youngs_modulus1 * (1 - poisson_ratio2 * poisson_ratio2)));

    effectiveshearmodulus[type1][type2] = (youngs_modulus1 * youngs_modulus2) /
                                          (2 * ((youngs_modulus2 * (2 - poisson_ratio1) *
                                                 (1 + poisson_ratio1)) +
                                                (youngs_modulus1 * (2 - poisson_ratio2) *
                                                 (1 + poisson_ratio2))));

    modelparameterbeta[type1][type2] = 2 * (contact_damping_coefficient1 * contact_damping_coefficient2) /
                                       (contact_damping_coefficient1 + contact_damping_coefficient2 + DBL_MIN);

    effectiveslidingfriction[type1][type2] = 2 * sliding_friction1 * sliding_friction2 /
                                             (sliding_friction1 + sliding_friction2 + DBL_MIN);

    effectiverollingfriction[type1][type2] = 2 * rolling_friction1 * rolling_friction2 /
                                             (rolling_friction1 + rolling_friction2 + DBL_MIN);
}
void ContactForce::computerSpherePlanewallEffective(std::shared_ptr<PlanewallProperties> &planeWallProps1, PropertyTypeID &type1,
                                                    std::shared_ptr<SphereProperties> &sphereProps2, PropertyTypeID &type2)
{
    const double youngs_modulus2 = sphereProps2->getYoungModulus();
    const double youngs_modulus1 = planeWallProps1->getYoungModulus();

    const double poisson_ratio2 = sphereProps2->getPoissonRatio();
    const double poisson_ratio1 = planeWallProps1->getPoissonRatio();

    const double restitution2 = sphereProps2->getRestitution();
    const double restitution1 = planeWallProps1->getRestitution();

    const double contact_damping_coefficient2 = -log(restitution2) / sqrt(PI * PI + log(restitution2) * log(restitution2));
    const double contact_damping_coefficient1 = -log(restitution1) / sqrt(PI * PI + log(restitution1) * log(restitution1));

    const double sliding_friction2 = sphereProps2->getSlidingFriction();
    const double sliding_friction1 = planeWallProps1->getSlidingFriction();

    const double rolling_friction2 = sphereProps2->getRollingFriction();
    const double rolling_friction1 = planeWallProps1->getRollingFriction();

    effectiveyoungsmodulus[type1][type2] = (youngs_modulus1 * youngs_modulus2) /
                                           ((youngs_modulus2 * (1 - poisson_ratio1 * poisson_ratio1)) +
                                            (youngs_modulus1 * (1 - poisson_ratio2 * poisson_ratio2)));

    effectiveshearmodulus[type1][type2] = (youngs_modulus1 * youngs_modulus2) /
                                          (2 * ((youngs_modulus2 * (2 - poisson_ratio1) *
                                                 (1 + poisson_ratio1)) +
                                                (youngs_modulus1 * (2 - poisson_ratio2) *
                                                 (1 + poisson_ratio2))));

    modelparameterbeta[type1][type2] = 2 * (contact_damping_coefficient1 * contact_damping_coefficient2) /
                                       (contact_damping_coefficient1 + contact_damping_coefficient2 + DBL_MIN);

    effectiveslidingfriction[type1][type2] = 2 * sliding_friction1 * sliding_friction2 /
                                             (sliding_friction1 + sliding_friction2 + DBL_MIN);

    effectiverollingfriction[type1][type2] = 2 * rolling_friction1 * rolling_friction2 /
                                             (rolling_friction1 + rolling_friction2 + DBL_MIN);
}
void ContactForce::computerSphereCylinderwallEffective(std::shared_ptr<CylinderwallProperties> &cylinderwallProps1, PropertyTypeID &type1,
                                                       std::shared_ptr<SphereProperties> &sphereProps2, PropertyTypeID &type2)
{
    const double youngs_modulus2 = sphereProps2->getYoungModulus();
    const double youngs_modulus1 = cylinderwallProps1->getYoungModulus();

    const double poisson_ratio2 = sphereProps2->getPoissonRatio();
    const double poisson_ratio1 = cylinderwallProps1->getPoissonRatio();

    const double restitution2 = sphereProps2->getRestitution();
    const double restitution1 = cylinderwallProps1->getRestitution();

    const double contact_damping_coefficient2 = -log(restitution2) / sqrt(PI * PI + log(restitution2) * log(restitution2));
    const double contact_damping_coefficient1 = -log(restitution1) / sqrt(PI * PI + log(restitution1) * log(restitution1));

    const double sliding_friction2 = sphereProps2->getSlidingFriction();
    const double sliding_friction1 = cylinderwallProps1->getSlidingFriction();

    const double rolling_friction2 = sphereProps2->getRollingFriction();
    const double rolling_friction1 = cylinderwallProps1->getRollingFriction();

    effectiveyoungsmodulus[type1][type2] = (youngs_modulus1 * youngs_modulus2) /
                                           ((youngs_modulus2 * (1 - poisson_ratio1 * poisson_ratio1)) +
                                            (youngs_modulus1 * (1 - poisson_ratio2 * poisson_ratio2)));

    effectiveshearmodulus[type1][type2] = (youngs_modulus1 * youngs_modulus2) /
                                          (2 * ((youngs_modulus2 * (2 - poisson_ratio1) *
                                                 (1 + poisson_ratio1)) +
                                                (youngs_modulus1 * (2 - poisson_ratio2) *
                                                 (1 + poisson_ratio2))));

    modelparameterbeta[type1][type2] = 2 * (contact_damping_coefficient1 * contact_damping_coefficient2) /
                                       (contact_damping_coefficient1 + contact_damping_coefficient2 + DBL_MIN);

    effectiveslidingfriction[type1][type2] = 2 * sliding_friction1 * sliding_friction2 /
                                             (sliding_friction1 + sliding_friction2 + DBL_MIN);

    effectiverollingfriction[type1][type2] = 2 * rolling_friction1 * rolling_friction2 /
                                             (rolling_friction1 + rolling_friction2 + DBL_MIN);
}
void ContactForce::computerFiberFiberEffective(std::shared_ptr<FiberProperties> &fibereproperties1, PropertyTypeID &type1,
                                               std::shared_ptr<FiberProperties> &fiberproperties2, PropertyTypeID &type2)
{
    // Handle sphere property
    const double radius1 = fibereproperties1->getRadius();
    const double radius2 = fiberproperties2->getRadius();

    const double youngs_modulus1 = fibereproperties1->getYoungModulus();
    const double youngs_modulus2 = fiberproperties2->getYoungModulus();

    const double poisson_ratio1 = fibereproperties1->getPoissonRatio();
    const double poisson_ratio2 = fiberproperties2->getPoissonRatio();

    const double restitution1 = fibereproperties1->getRestitution();
    const double restitution2 = fiberproperties2->getRestitution();

    const double contact_damping_coefficient1 = -log(restitution1) / sqrt(PI * PI + log(restitution1) * log(restitution1));
    const double contact_damping_coefficient2 = -log(restitution2) / sqrt(PI * PI + log(restitution2) * log(restitution2));

    const double sliding_friction1 = fibereproperties1->getSlidingFriction();
    const double sliding_friction2 = fiberproperties2->getSlidingFriction();

    const double rolling_friction1 = fibereproperties1->getRollingFriction();
    const double rolling_friction2 = fiberproperties2->getRollingFriction();

    const double mass1 = fibereproperties1->getNodemass();
    const double mass2 = fiberproperties2->getNodemass();

    effectiveradius[type1][type2] = 2 * (radius1 * radius2) / (radius1 + radius2);

    effectivemass[type1][type2] = 2 * (mass1 * mass2) / (mass1 + mass2);

    effectiveyoungsmodulus[type1][type2] = (youngs_modulus1 * youngs_modulus2) /
                                           ((youngs_modulus2 * (1 - poisson_ratio1 * poisson_ratio1)) +
                                            (youngs_modulus1 * (1 - poisson_ratio2 * poisson_ratio2)));

    effectiveshearmodulus[type1][type2] = (youngs_modulus1 * youngs_modulus2) /
                                          (2 * ((youngs_modulus2 * (2 - poisson_ratio1) *
                                                 (1 + poisson_ratio1)) +
                                                (youngs_modulus1 * (2 - poisson_ratio2) *
                                                 (1 + poisson_ratio2))));

    modelparameterbeta[type1][type2] = 2 * (contact_damping_coefficient1 * contact_damping_coefficient2) /
                                       (contact_damping_coefficient1 + contact_damping_coefficient2 + DBL_MIN);

    effectiveslidingfriction[type1][type2] = 2 * sliding_friction1 * sliding_friction2 /
                                             (sliding_friction1 + sliding_friction2 + DBL_MIN);

    effectiverollingfriction[type1][type2] = 2 * rolling_friction1 * rolling_friction2 /
                                             (rolling_friction1 + rolling_friction2 + DBL_MIN);
}
void ContactForce::computerSphereFiberEffective(std::shared_ptr<SphereProperties> &sphereProps1, PropertyTypeID &type1,
                                                std::shared_ptr<FiberProperties> &fiberproperties2, PropertyTypeID &type2)
{
    const double radius1 = sphereProps1->getRadius();
    const double radius2 = fiberproperties2->getRadius();

    const double youngs_modulus1 = sphereProps1->getYoungModulus();
    const double youngs_modulus2 = fiberproperties2->getYoungModulus();

    const double poisson_ratio1 = sphereProps1->getPoissonRatio();
    const double poisson_ratio2 = fiberproperties2->getPoissonRatio();

    const double restitution1 = sphereProps1->getRestitution();
    const double restitution2 = fiberproperties2->getRestitution();

    const double contact_damping_coefficient1 = -log(restitution1) / sqrt(PI * PI + log(restitution1) * log(restitution1));
    const double contact_damping_coefficient2 = -log(restitution2) / sqrt(PI * PI + log(restitution2) * log(restitution2));

    const double sliding_friction1 = sphereProps1->getSlidingFriction();
    const double sliding_friction2 = fiberproperties2->getSlidingFriction();

    const double rolling_friction1 = sphereProps1->getRollingFriction();
    const double rolling_friction2 = fiberproperties2->getRollingFriction();

    const double mass1 = sphereProps1->getMass();
    const double mass2 = fiberproperties2->getNodemass();

    effectiveradius[type1][type2] = 2 * (radius1 * radius2) / (radius1 + radius2);

    effectivemass[type1][type2] = 2 * (mass1 * mass2) / (mass1 + mass2);

    effectiveyoungsmodulus[type1][type2] = (youngs_modulus1 * youngs_modulus2) /
                                           ((youngs_modulus2 * (1 - poisson_ratio1 * poisson_ratio1)) +
                                            (youngs_modulus1 * (1 - poisson_ratio2 * poisson_ratio2)));

    effectiveshearmodulus[type1][type2] = (youngs_modulus1 * youngs_modulus2) /
                                          (2 * ((youngs_modulus2 * (2 - poisson_ratio1) *
                                                 (1 + poisson_ratio1)) +
                                                (youngs_modulus1 * (2 - poisson_ratio2) *
                                                 (1 + poisson_ratio2))));

    modelparameterbeta[type1][type2] = 2 * (contact_damping_coefficient1 * contact_damping_coefficient2) /
                                       (contact_damping_coefficient1 + contact_damping_coefficient2 + DBL_MIN);

    effectiveslidingfriction[type1][type2] = 2 * sliding_friction1 * sliding_friction2 /
                                             (sliding_friction1 + sliding_friction2 + DBL_MIN);

    effectiverollingfriction[type1][type2] = 2 * rolling_friction1 * rolling_friction2 /
                                             (rolling_friction1 + rolling_friction2 + DBL_MIN);
}
void ContactForce::computerFiberPlanewallEffective(std::shared_ptr<PlanewallProperties> &planeWallProps1, PropertyTypeID &type1,
                                                   std::shared_ptr<FiberProperties> &fiberproperties2, PropertyTypeID &type2)
{
    const double youngs_modulus2 = fiberproperties2->getYoungModulus();
    const double youngs_modulus1 = planeWallProps1->getYoungModulus();

    const double poisson_ratio2 = fiberproperties2->getPoissonRatio();
    const double poisson_ratio1 = planeWallProps1->getPoissonRatio();

    const double restitution2 = fiberproperties2->getRestitution();
    const double restitution1 = planeWallProps1->getRestitution();

    const double contact_damping_coefficient2 = -log(restitution2) / sqrt(PI * PI + log(restitution2) * log(restitution2));
    const double contact_damping_coefficient1 = -log(restitution1) / sqrt(PI * PI + log(restitution1) * log(restitution1));

    const double sliding_friction2 = fiberproperties2->getSlidingFriction();
    const double sliding_friction1 = planeWallProps1->getSlidingFriction();

    const double rolling_friction2 = fiberproperties2->getRollingFriction();
    const double rolling_friction1 = planeWallProps1->getRollingFriction();

    effectiveyoungsmodulus[type1][type2] = (youngs_modulus1 * youngs_modulus2) /
                                           ((youngs_modulus2 * (1 - poisson_ratio1 * poisson_ratio1)) +
                                            (youngs_modulus1 * (1 - poisson_ratio2 * poisson_ratio2)));

    effectiveshearmodulus[type1][type2] = (youngs_modulus1 * youngs_modulus2) /
                                          (2 * ((youngs_modulus2 * (2 - poisson_ratio1) *
                                                 (1 + poisson_ratio1)) +
                                                (youngs_modulus1 * (2 - poisson_ratio2) *
                                                 (1 + poisson_ratio2))));

    modelparameterbeta[type1][type2] = 2 * (contact_damping_coefficient1 * contact_damping_coefficient2) /
                                       (contact_damping_coefficient1 + contact_damping_coefficient2 + DBL_MIN);

    effectiveslidingfriction[type1][type2] = 2 * sliding_friction1 * sliding_friction2 /
                                             (sliding_friction1 + sliding_friction2 + DBL_MIN);

    effectiverollingfriction[type1][type2] = 2 * rolling_friction1 * rolling_friction2 /
                                             (rolling_friction1 + rolling_friction2 + DBL_MIN);
}

void ContactForce::computerFiberCylinderwallEffective(std::shared_ptr<CylinderwallProperties> &cylinderwallproperties1, PropertyTypeID &type1,
                                                      std::shared_ptr<FiberProperties> &fiberproperties2, PropertyTypeID &type2)
{
    const double youngs_modulus2 = fiberproperties2->getYoungModulus();
    const double youngs_modulus1 = cylinderwallproperties1->getYoungModulus();

    const double poisson_ratio2 = fiberproperties2->getPoissonRatio();
    const double poisson_ratio1 = cylinderwallproperties1->getPoissonRatio();

    const double restitution2 = fiberproperties2->getRestitution();
    const double restitution1 = cylinderwallproperties1->getRestitution();

    const double contact_damping_coefficient2 = -log(restitution2) / sqrt(PI * PI + log(restitution2) * log(restitution2));
    const double contact_damping_coefficient1 = -log(restitution1) / sqrt(PI * PI + log(restitution1) * log(restitution1));

    const double sliding_friction2 = fiberproperties2->getSlidingFriction();
    const double sliding_friction1 = cylinderwallproperties1->getSlidingFriction();

    const double rolling_friction2 = fiberproperties2->getRollingFriction();
    const double rolling_friction1 = cylinderwallproperties1->getRollingFriction();

    effectiveyoungsmodulus[type1][type2] = (youngs_modulus1 * youngs_modulus2) /
                                           ((youngs_modulus2 * (1 - poisson_ratio1 * poisson_ratio1)) +
                                            (youngs_modulus1 * (1 - poisson_ratio2 * poisson_ratio2)));

    effectiveshearmodulus[type1][type2] = (youngs_modulus1 * youngs_modulus2) /
                                          (2 * ((youngs_modulus2 * (2 - poisson_ratio1) *
                                                 (1 + poisson_ratio1)) +
                                                (youngs_modulus1 * (2 - poisson_ratio2) *
                                                 (1 + poisson_ratio2))));

    modelparameterbeta[type1][type2] = 2 * (contact_damping_coefficient1 * contact_damping_coefficient2) /
                                       (contact_damping_coefficient1 + contact_damping_coefficient2 + DBL_MIN);

    effectiveslidingfriction[type1][type2] = 2 * sliding_friction1 * sliding_friction2 /
                                             (sliding_friction1 + sliding_friction2 + DBL_MIN);

    effectiverollingfriction[type1][type2] = 2 * rolling_friction1 * rolling_friction2 /
                                             (rolling_friction1 + rolling_friction2 + DBL_MIN);
}

void ContactForce::addParticleProperties(const std::shared_ptr<ParticlePropertyManager> &manager)
{

    for (auto &property1 : manager->getParticleProperties())
    {
        PropertyTypeID type1 = property1.first;

        for (auto &property2 : manager->getParticleProperties())
        {
            PropertyTypeID type2 = property2.first;
            if (auto &sphereProps1 = std::dynamic_pointer_cast<SphereProperties>(property1.second))
            {
                if (auto &sphereProps2 = std::dynamic_pointer_cast<SphereProperties>(property2.second))
                {
                    computerSphereSphereEffective(sphereProps1, type1, sphereProps2, type2);
                }
                else if (auto &fiberProps2 = std::dynamic_pointer_cast<FiberProperties>(property2.second))
                {
                    computerSphereFiberEffective(sphereProps1, type1, fiberProps2, type2);
                }
                else if (auto &planeWallProps2 = std::dynamic_pointer_cast<PlanewallProperties>(property2.second))
                {
                    computerSpherePlanewallEffective(planeWallProps2, type2, sphereProps1, type1);
                }
                else if (auto &cylinderWallProps2 = std::dynamic_pointer_cast<CylinderwallProperties>(property2.second))
                {
                    computerSphereCylinderwallEffective(cylinderWallProps2, type2, sphereProps1, type1);
                }
            }
            else if (auto &fiberProps1 = std::dynamic_pointer_cast<FiberProperties>(property1.second))
            {
                if (auto &sphereProps2 = std::dynamic_pointer_cast<SphereProperties>(property2.second))
                {
                    computerSphereFiberEffective(sphereProps2, type2, fiberProps1, type1);
                }
                else if (auto &fiberProps2 = std::dynamic_pointer_cast<FiberProperties>(property2.second))
                {
                    computerFiberFiberEffective(fiberProps1, type1, fiberProps2, type2);
                }
                else if (auto &planeWallProps2 = std::dynamic_pointer_cast<PlanewallProperties>(property2.second))
                {
                    computerFiberPlanewallEffective(planeWallProps2, type2, fiberProps1, type1);
                }
                else if (auto &cylinderWallProps2 = std::dynamic_pointer_cast<CylinderwallProperties>(property2.second))
                {
                    computerFiberCylinderwallEffective(cylinderWallProps2, type2, fiberProps1, type1);
                }
            }
        }
    }
}

void ContactForce::computeSphereSphereForce(const std::unique_ptr<SphereParticle> &sphere1, const std::unique_ptr<SphereParticle> &sphere2, double timeStep)
{

    auto &manager = sphere1->getParticlePropertyManager();

    PropertyTypeID particletype_one = sphere1->getType();
    PropertyTypeID particletype_two = sphere2->getType();

    auto radius1 = manager->getSphereProperties(particletype_one)->getRadius();
    auto radius2 = manager->getSphereProperties(particletype_two)->getRadius();

    // Compute the distance between the centers of the two spheres
    Eigen::Vector3d position1 = sphere1->getPosition();
    Eigen::Vector3d position2 = sphere2->getPosition();

    Eigen::Vector3d normal_unit_vector = position2 - position1;
    double distance = normal_unit_vector.norm();

    // Compute the overlap (positive if spheres intersect, zero or negative otherwise)
    double normal_overlap = (radius1 + radius2) - distance;

    // The Hertz-Mindlin  contact model

    if (normal_overlap > 0)
    {
        normal_unit_vector.normalize();

        Eigen::Vector3d v1 = sphere1->getVelocity();

        Eigen::Vector3d v2 = sphere2->getVelocity();

        Eigen::Vector3d w1 = sphere1->getOmega();

        Eigen::Vector3d w2 = sphere2->getOmega();

        Eigen::Vector3d contact_relative_velocity = (v1 - v2) + (radius1 * w1 + radius2 * w2).cross(normal_unit_vector);

        double normal_relative_velocity_value = contact_relative_velocity.dot(normal_unit_vector);

        // Calculation of normal relative velocity
        Eigen::Vector3d normal_relative_velocity = normal_relative_velocity_value * normal_unit_vector;
        // Calculation of tangential relative velocity
        Eigen::Vector3d tangential_relative_velocity = contact_relative_velocity - normal_relative_velocity;
        // modify tangential_overlap
        int id1 = sphere1->getId();
        int id2 = sphere2->getId();

        auto &innerMap = spherespherecontactinformationlist[id1]; // Access or create the inner map
        auto innerIt = innerMap.find(id2);

        if (innerIt == innerMap.end())
        {
            // Contact information for this pair does not exist, create a new one
            ContactInformation newContactInfo(id1, id2);
            innerMap[id2] = newContactInfo; // Insert the new contact information
        }

        ContactInformation &contact_info = spherespherecontactinformationlist[id1][id2];

        Eigen::Vector3d modified_tangential_displacement = contact_info.tangentialDisplacement + contact_info.previousTangentialVelocity * timeStep;
        // Updating the contact_info container based on the new calculated values
        contact_info.tangentialDisplacement = modified_tangential_displacement;
        contact_info.previousTangentialVelocity = tangential_relative_velocity;

        double radius_times_overlap_sqrt = sqrt(effectiveradius[particletype_one][particletype_two] * normal_overlap);

        double sn = 2 * effectiveyoungsmodulus[particletype_one][particletype_two] * radius_times_overlap_sqrt;

        double normal_contact_stiffness = 0.666666 * sn;

        double normal_damping_constant = 1.8257 * modelparameterbeta[particletype_one][particletype_two] * sqrt(sn * effectivemass[particletype_one][particletype_two]);

        // tangential stiffness
        double tangential_contact_stiffness = 8 * effectiveshearmodulus[particletype_one][particletype_two] * radius_times_overlap_sqrt;

        double tangential_damping_constant = 1.8257 * modelparameterbeta[particletype_one][particletype_two] * sqrt(tangential_contact_stiffness * effectivemass[particletype_one][particletype_two]);
        Eigen::Vector3d normal_damping_force = (normal_damping_constant * normal_relative_velocity_value) * normal_unit_vector;
        // Calculation of normal force using spring and dashpot normal forces
        Eigen::Vector3d normal_force =
            ((normal_contact_stiffness * normal_overlap) * normal_unit_vector) + normal_damping_force;

        energyDissipation += std::fabs(normal_damping_force.dot(normal_relative_velocity * timeStep));
        // Calculation of tangential force using spring and dashpot tangential
        // forces. Since we need dashpot tangential force in the gross sliding again,
        // we define it as a separate variable
        Eigen::Vector3d dashpot_tangential_force = tangential_damping_constant * tangential_relative_velocity;
        Eigen::Vector3d tangential_force = (tangential_contact_stiffness * modified_tangential_displacement) + dashpot_tangential_force;

        double coulomb_threshold = effectiveslidingfriction[particletype_one][particletype_one] * normal_force.norm();

        double tangential_force_value = tangential_force.norm();
        // Check for gross sliding
        if (tangential_force_value > coulomb_threshold)
        {
            // Gross sliding occurs and the tangential overlap and tangnetial
            // force are limited to Coulumb's criterion
            Eigen::Vector3d tangential_overlap = (coulomb_threshold * (tangential_force / (tangential_force_value + DBL_MIN)) - dashpot_tangential_force) /
                                                 (tangential_contact_stiffness + DBL_MIN);

            contact_info.tangentialDisplacement = tangential_overlap;
            tangential_force = (tangential_contact_stiffness * tangential_overlap) + dashpot_tangential_force;
        }
        energyDissipation += std::fabs(tangential_force.dot(tangential_relative_velocity * timeStep));

        // Calculation of torque
        // Torque caused by tangential force (tangential_torque)
        Eigen::Vector3d tangential_torque1 = (radius1 * normal_unit_vector).cross(tangential_force);
        Eigen::Vector3d tangential_torque2 = (radius2 * normal_unit_vector).cross(tangential_force);

        Eigen::Vector3d results_force = normal_force + tangential_force;

        sphere1->addForce(-results_force);
        sphere2->addForce(results_force);

        sphere1->addTorque(-tangential_torque1);
        sphere2->addTorque(-tangential_torque2);
        contact_info.normalForce = -normal_force;
        contact_info.tangentialForce = -tangential_force;
        nextspherespherecontactinformationlist[id1][id2] = contact_info;
    }
}

void ContactForce::computeSphereFiberForce(const std::unique_ptr<SphereParticle> &sphere, const std::unique_ptr<SphereCylinderBond> &fiberbond,
                                           const std::unique_ptr<SphereParticle> &node1, const std::unique_ptr<SphereParticle> &node2, double timeStep)
{
    auto &manager = fiberbond->getParticlePropertyManager();

    const PropertyTypeID spheretype = sphere->getType();
    const PropertyTypeID fibertype = fiberbond->getType();

    double sphereRadius = manager->getSphereProperties(spheretype)->getRadius();
    double fiberRadius = manager->getFiberProperties(fibertype)->getRadius();

    const Eigen::Vector3d &spherePos = sphere->getPosition();
    const Eigen::Vector3d &node1Pos = node1->getPosition();
    const Eigen::Vector3d &node2Pos = node2->getPosition();

    double t;
    Eigen::Vector3d projection;
    SphereCylinderBond::computeOverlap(node1Pos, node2Pos, spherePos, t, projection);
    Eigen::Vector3d normal_unit_vector = projection - spherePos;
    double distance = normal_unit_vector.norm();
    double normal_overlap = (sphereRadius + fiberRadius) - distance;

    if (normal_overlap > 0)
    {
        normal_unit_vector.normalize();

        const Eigen::Vector3d &spherevelocity = sphere->getVelocity();
        const Eigen::Vector3d &node1Velocity = node1->getVelocity();
        const Eigen::Vector3d &node2Velocity = node2->getVelocity();

        const Eigen::Vector3d &sphereomega = sphere->getOmega();
        const Eigen::Vector3d &node1Omega = node1->getOmega();
        const Eigen::Vector3d &node2Omega = node2->getOmega();

        Eigen::Vector3d fibervelocity = t * (node2Velocity - node1Velocity) + node1Velocity;
        Eigen::Vector3d fiberomega = t * (node2Omega - node1Omega) + node1Omega;

        Eigen::Vector3d contact_relative_velocity = (spherevelocity - fibervelocity) + (sphereRadius * sphereomega + fiberRadius * fiberomega).cross(normal_unit_vector);
        double normal_relative_velocity_value = contact_relative_velocity.dot(normal_unit_vector);

        // Calculation of normal relative velocity
        Eigen::Vector3d normal_relative_velocity = normal_relative_velocity_value * normal_unit_vector;
        // Calculation of tangential relative velocity
        Eigen::Vector3d tangential_relative_velocity = contact_relative_velocity - normal_relative_velocity;
        // modify tangential_overlap
        int id1 = sphere->getId();
        int id2 = fiberbond->getId();

        auto &innerMap = spherefibercontactinformationlist[id1]; // Access or create the inner map
        auto innerIt = innerMap.find(id2);

        if (innerIt == innerMap.end())
        {
            // Contact information for this pair does not exist, create a new one
            ContactInformation newContactInfo(id1, id2);
            innerMap[id2] = newContactInfo; // Insert the new contact information
        }

        ContactInformation &contact_info = spherefibercontactinformationlist[id1][id2];

        Eigen::Vector3d modified_tangential_displacement = contact_info.tangentialDisplacement + contact_info.previousTangentialVelocity * timeStep;
        // Updating the contact_info container based on the new calculated values
        contact_info.tangentialDisplacement = modified_tangential_displacement;
        contact_info.previousTangentialVelocity = tangential_relative_velocity;

        double radius_times_overlap_sqrt = sqrt(effectiveradius[spheretype][fibertype] * normal_overlap);

        double sn = 2 * effectiveyoungsmodulus[spheretype][fibertype] * radius_times_overlap_sqrt;

        double normal_contact_stiffness = 0.666666 * sn;

        double normal_damping_constant = 1.8257 * modelparameterbeta[spheretype][fibertype] * sqrt(sn * effectivemass[spheretype][fibertype]);

        // tangential stiffness
        double tangential_contact_stiffness = 8 * effectiveshearmodulus[spheretype][fibertype] * radius_times_overlap_sqrt;

        double tangential_damping_constant = 1.8257 * modelparameterbeta[spheretype][fibertype] * sqrt(tangential_contact_stiffness * effectivemass[spheretype][fibertype]);

        Eigen::Vector3d normal_damping_force = (normal_damping_constant * normal_relative_velocity_value) * normal_unit_vector;
        // Calculation of normal force using spring and dashpot normal forces
        Eigen::Vector3d normal_force =
            ((normal_contact_stiffness * normal_overlap) * normal_unit_vector) + normal_damping_force;

        energyDissipation += std::fabs(normal_damping_force.dot(normal_relative_velocity * timeStep));

        // Calculation of tangential force using spring and dashpot tangential
        // forces. Since we need dashpot tangential force in the gross sliding again,
        // we define it as a separate variable
        Eigen::Vector3d dashpot_tangential_force = tangential_damping_constant * tangential_relative_velocity;
        Eigen::Vector3d tangential_force = (tangential_contact_stiffness * modified_tangential_displacement) + dashpot_tangential_force;

        double coulomb_threshold = effectiveslidingfriction[spheretype][fibertype] * normal_force.norm();

        double tangential_force_value = tangential_force.norm();
        // Check for gross sliding
        if (tangential_force_value > coulomb_threshold)
        {
            // Gross sliding occurs and the tangential overlap and tangnetial
            // force are limited to Coulumb's criterion
            Eigen::Vector3d tangential_overlap = (coulomb_threshold * (tangential_force / (tangential_force_value + DBL_MIN)) - dashpot_tangential_force) /
                                                 (tangential_contact_stiffness + DBL_MIN);

            contact_info.tangentialDisplacement = tangential_overlap;
            tangential_force = (tangential_contact_stiffness * tangential_overlap) + dashpot_tangential_force;
        }
        energyDissipation += std::fabs(tangential_force.dot(tangential_relative_velocity * timeStep));

        // Calculation of torque
        // Torque caused by tangential force (tangential_torque)
        Eigen::Vector3d tangential_torque1 = (sphereRadius * normal_unit_vector).cross(tangential_force);
        Eigen::Vector3d tangential_torque2 = (fiberRadius * normal_unit_vector).cross(tangential_force);

        Eigen::Vector3d results_force = normal_force + tangential_force;

        sphere->addForce(-results_force);
        node1->addForce((1 - t) * results_force);
        node2->addForce(t * results_force);

        sphere->addTorque(-tangential_torque1);
        node1->addTorque((t - 1) * tangential_torque2);
        node2->addTorque(-t * tangential_torque2);
        contact_info.normalForce = -normal_force;
        contact_info.tangentialForce = -tangential_force;
        nextspherefibercontactinformationlist[id1][id2] = contact_info;
    }
}

void ContactForce::computePlaneWallSphereForce(const std::unique_ptr<PlaneWall> &planewall, const std::unique_ptr<SphereParticle> &sphere, double timeStep)
{
    // Retrieve the sphere's center and radius
    Eigen::Vector3d sphereCenter = sphere->getPosition();

    auto &manager = sphere->getParticlePropertyManager();

    PropertyTypeID spheretype = sphere->getType();
    PropertyTypeID walltype = planewall->getType();

    double sphereRadius = manager->getSphereProperties(spheretype)->getRadius();

    // Retrieve the plane wall's point and normal
    Eigen::Vector3d planePoint = planewall->getCorner1();

    Eigen::Vector3d planeNormal = planewall->getNormal();

    // Compute the vector from a point on the plane to the sphere's center
    Eigen::Vector3d vecToSphereCenter = sphereCenter - planePoint;
    // Compute the distance from the sphere's center to the plane
    double distanceToPlane = vecToSphereCenter.dot(planeNormal);

    if (distanceToPlane < 0)
    {
        planeNormal = -planeNormal;
        distanceToPlane = -distanceToPlane;
    }

    double normal_overlap = sphereRadius - distanceToPlane;

    if (normal_overlap > 0)
    {
        // Using contact_vector, the contact normal vector is obtained
        Eigen::Vector3d normal_unit_vector = -planeNormal;
        // Defining velocities and angular velocities of sphere
        Eigen::Vector3d sphere_velocity = sphere->getVelocity();

        Eigen::Vector3d sphere_omega = sphere->getOmega();
        // Calculation of contact relative velocity
        Eigen::Vector3d contact_relative_velocity = sphere_velocity + (sphereRadius * sphere_omega).cross(normal_unit_vector);

        double normal_relative_velocity_value = contact_relative_velocity.dot(normal_unit_vector);
        // Calculation of normal relative velocity
        Eigen::Vector3d normal_relative_velocity = normal_relative_velocity_value * normal_unit_vector;
        // Calculation of tangential relative velocity
        Eigen::Vector3d tangential_relative_velocity = contact_relative_velocity - normal_relative_velocity;
        // modify tangential_overlap
        int wallId = planewall->getId();
        int sphereId = sphere->getId();

        auto &innerMap = planewallspherecontactinformationlist[wallId]; // Access or create the inner map
        auto innerIt = innerMap.find(sphereId);

        if (innerIt == innerMap.end())
        {
            // Contact information for this pair does not exist, create a new one
            ContactInformation newContactInfo(wallId, sphereId);
            innerMap[sphereId] = newContactInfo; // Insert the new contact information
        }

        ContactInformation &contact_info = planewallspherecontactinformationlist[wallId][sphereId];

        Eigen::Vector3d modified_tangential_displacement = contact_info.tangentialDisplacement + contact_info.previousTangentialVelocity * timeStep;
        // Updating the contact_info container based on the new calculated values

        contact_info.tangentialDisplacement = modified_tangential_displacement;
        contact_info.previousTangentialVelocity = tangential_relative_velocity;

        double radius_times_overlap_sqrt = sqrt(manager->getSphereProperties(spheretype)->getRadius() * normal_overlap);

        double sn = 2 * effectiveyoungsmodulus[walltype][spheretype] * radius_times_overlap_sqrt;

        double normal_contact_stiffness = 0.666666 * sn;

        double normal_damping_constant = 1.8257 * modelparameterbeta[walltype][spheretype] * sqrt(sn * manager->getSphereProperties(spheretype)->getMass());

        // tangential stiffness
        double tangential_contact_stiffness = 8 * effectiveshearmodulus[walltype][spheretype] * radius_times_overlap_sqrt;

        double tangential_damping_constant = 1.8257 * modelparameterbeta[walltype][spheretype] * sqrt(tangential_contact_stiffness * manager->getSphereProperties(spheretype)->getMass());

        Eigen::Vector3d normal_damping_force = (normal_damping_constant * normal_relative_velocity_value) * normal_unit_vector;
        // Calculation of normal force using spring and dashpot normal forces
        Eigen::Vector3d normal_force =
            ((normal_contact_stiffness * normal_overlap) * normal_unit_vector) + normal_damping_force;

        energyDissipation += std::fabs(normal_damping_force.dot(normal_relative_velocity * timeStep));

        // Calculation of tangential force using spring and dashpot tangential
        // forces. Since we need dashpot tangential force in the gross sliding again,
        // we define it as a separate variable
        Eigen::Vector3d dashpot_tangential_force = tangential_damping_constant * tangential_relative_velocity;
        Eigen::Vector3d tangential_force = (tangential_contact_stiffness * modified_tangential_displacement) + dashpot_tangential_force;

        double coulomb_threshold = effectiveslidingfriction[walltype][spheretype] * normal_force.norm();

        double tangential_force_value = tangential_force.norm();
        // Check for gross sliding
        if (tangential_force_value > coulomb_threshold)
        {
            // Gross sliding occurs and the tangential overlap and tangnetial
            // force are limited to Coulumb's criterion
            Eigen::Vector3d tangential_overlap = (coulomb_threshold * (tangential_force / (tangential_force_value + DBL_MIN)) - dashpot_tangential_force) /
                                                 (tangential_contact_stiffness + DBL_MIN);

            contact_info.tangentialDisplacement = tangential_overlap;
            tangential_force = (tangential_contact_stiffness * tangential_overlap) + dashpot_tangential_force;
        }
        energyDissipation += std::fabs(tangential_force.dot(tangential_relative_velocity * timeStep));

        // Calculation of torque
        // Torque caused by tangential force (tangential_torque)
        Eigen::Vector3d tangential_torque1 = (sphereRadius * normal_unit_vector).cross(tangential_force);

        Eigen::Vector3d results_force = normal_force + tangential_force;

        sphere->addForce(-results_force);
        sphere->addTorque(-tangential_torque1);
        planewall->addForce(results_force);

        contact_info.normalForce = normal_force;
        contact_info.tangentialForce = tangential_force;
        nextplanewallspherecontactinformationlist[wallId][sphereId] = contact_info;
    }
}

void ContactForce::computeCylinderWallSphereForce(const std::unique_ptr<CylinderContainer> &cylinderWall, const std::unique_ptr<SphereParticle> &sphere, double timeStep)
{
    // Retrieve the sphere's center and radius
    Eigen::Vector3d sphereCenter = sphere->getPosition();

    auto &manager = sphere->getParticlePropertyManager();

    PropertyTypeID spheretype = sphere->getType();
    PropertyTypeID walltype = cylinderWall->getType();

    double sphereRadius = manager->getSphereProperties(spheretype)->getRadius();
    double cylinderRadius = cylinderWall->getRadius();
    // Retrieve the cylinder wall's point
    Eigen::Vector3d cylinderEnd = cylinderWall->getEndpoint();

    Eigen::Vector3d cylinderStart = cylinderWall->getStartpoint();

    // Calculate the axis direction of the cylinder
    Eigen::Vector3d cylinderAxis = (cylinderEnd - cylinderStart).normalized();

    // Calculate the vector from the start point of the cylinder to the sphere center
    Eigen::Vector3d startToCenter = sphereCenter - cylinderStart;

    // Project the vector onto the cylinder axis to find the closest point on the axis to the sphere center
    double projection = startToCenter.dot(cylinderAxis);
    projection = std::clamp(projection, 0.0, 1.0);
    Eigen::Vector3d closestPoint = cylinderStart + projection * cylinderAxis;
    Eigen::Vector3d normal_unit_vector = sphereCenter - closestPoint;

    // Calculate the distance from the sphere center to the closest point on the axis
    double distanceToAxis = normal_unit_vector.norm();
    double normal_overlap = 0;
    if (distanceToAxis > cylinderRadius)
    {
        normal_unit_vector = -normal_unit_vector / distanceToAxis;
        normal_overlap = sphereRadius + cylinderRadius - distanceToAxis;
    }
    else
    {
        normal_unit_vector = normal_unit_vector / distanceToAxis;
        normal_overlap = distanceToAxis - (cylinderRadius - sphereRadius);
    }
    if (normal_overlap > 0)
    {
        // Defining velocities and angular velocities of sphere
        Eigen::Vector3d sphere_velocity = sphere->getVelocity();

        Eigen::Vector3d sphere_omega = sphere->getOmega();
        // Calculation of contact relative velocity
        Eigen::Vector3d contact_relative_velocity = sphere_velocity + (sphereRadius * sphere_omega).cross(normal_unit_vector);

        double normal_relative_velocity_value = contact_relative_velocity.dot(normal_unit_vector);
        // Calculation of normal relative velocity
        Eigen::Vector3d normal_relative_velocity = normal_relative_velocity_value * normal_unit_vector;
        // Calculation of tangential relative velocity
        Eigen::Vector3d tangential_relative_velocity = contact_relative_velocity - normal_relative_velocity;
        // modify tangential_overlap
        int wallId = cylinderWall->getId();
        int sphereId = sphere->getId();

        auto &innerMap = cylinderwallspherecontactinformationlist[wallId]; // Access or create the inner map
        auto innerIt = innerMap.find(sphereId);

        if (innerIt == innerMap.end())
        {
            // Contact information for this pair does not exist, create a new one
            ContactInformation newContactInfo(wallId, sphereId);
            innerMap[sphereId] = newContactInfo; // Insert the new contact information
        }

        ContactInformation &contact_info = cylinderwallspherecontactinformationlist[wallId][sphereId];

        Eigen::Vector3d modified_tangential_displacement = contact_info.tangentialDisplacement + contact_info.previousTangentialVelocity * timeStep;
        // Updating the contact_info container based on the new calculated values

        contact_info.tangentialDisplacement = modified_tangential_displacement;
        contact_info.previousTangentialVelocity = tangential_relative_velocity;

        double radius_times_overlap_sqrt = sqrt(manager->getSphereProperties(spheretype)->getRadius() * normal_overlap);

        double sn = 2 * effectiveyoungsmodulus[walltype][spheretype] * radius_times_overlap_sqrt;

        double normal_contact_stiffness = 0.666666 * sn;

        double normal_damping_constant = 1.8257 * modelparameterbeta[walltype][spheretype] * sqrt(sn * manager->getSphereProperties(spheretype)->getMass());

        // tangential stiffness
        double tangential_contact_stiffness = 8 * effectiveshearmodulus[walltype][spheretype] * radius_times_overlap_sqrt;

        double tangential_damping_constant = 1.8257 * modelparameterbeta[walltype][spheretype] * sqrt(tangential_contact_stiffness * manager->getSphereProperties(spheretype)->getMass());

        Eigen::Vector3d normal_damping_force = (normal_damping_constant * normal_relative_velocity_value) * normal_unit_vector;
        // Calculation of normal force using spring and dashpot normal forces
        Eigen::Vector3d normal_force =
            ((normal_contact_stiffness * normal_overlap) * normal_unit_vector) + normal_damping_force;

        energyDissipation += std::fabs(normal_damping_force.dot(normal_relative_velocity * timeStep));

        // Calculation of tangential force using spring and dashpot tangential
        // forces. Since we need dashpot tangential force in the gross sliding again,
        // we define it as a separate variable
        Eigen::Vector3d dashpot_tangential_force = tangential_damping_constant * tangential_relative_velocity;
        Eigen::Vector3d tangential_force = (tangential_contact_stiffness * modified_tangential_displacement) + dashpot_tangential_force;

        double coulomb_threshold = effectiveslidingfriction[walltype][spheretype] * normal_force.norm();

        double tangential_force_value = tangential_force.norm();
        // Check for gross sliding
        if (tangential_force_value > coulomb_threshold)
        {
            // Gross sliding occurs and the tangential overlap and tangnetial
            // force are limited to Coulumb's criterion
            Eigen::Vector3d tangential_overlap = (coulomb_threshold * (tangential_force / (tangential_force_value + DBL_MIN)) - dashpot_tangential_force) /
                                                 (tangential_contact_stiffness + DBL_MIN);

            contact_info.tangentialDisplacement = tangential_overlap;
            tangential_force = (tangential_contact_stiffness * tangential_overlap) + dashpot_tangential_force;
        }
        energyDissipation += std::fabs(tangential_force.dot(tangential_relative_velocity * timeStep));

        // Calculation of torque
        // Torque caused by tangential force (tangential_torque)
        Eigen::Vector3d tangential_torque1 = (sphereRadius * normal_unit_vector).cross(tangential_force);

        Eigen::Vector3d results_force = normal_force + tangential_force;

        sphere->addForce(-results_force);
        sphere->addTorque(-tangential_torque1);
        cylinderWall->addForce(results_force);

        contact_info.normalForce = normal_force;
        contact_info.tangentialForce = tangential_force;
        nextcylinderwallspherecontactinformationlist[wallId][sphereId] = contact_info;
    }
}

void ContactForce::computeFiberFiberForce(const std::unique_ptr<SphereCylinderBond> &fiberbond1, const std::unique_ptr<SphereParticle> &bondonenode1, const std::unique_ptr<SphereParticle> &bondonenode2,
                                          const std::unique_ptr<SphereCylinderBond> &fiberbond2, const std::unique_ptr<SphereParticle> &bondtwonode1, const std::unique_ptr<SphereParticle> &bondtwonode2,
                                          double timeStep)
{
    auto &manager = fiberbond1->getParticlePropertyManager();

    const PropertyTypeID particletype_one = fiberbond1->getType();
    const PropertyTypeID particletype_two = fiberbond2->getType();

    double fiber1Radius = manager->getFiberProperties(particletype_one)->getRadius();
    double fiber2Radius = manager->getFiberProperties(particletype_two)->getRadius();

    const Eigen::Vector3d &bondonenode1Pos = bondonenode1->getPosition();
    const Eigen::Vector3d &bondonenode2Pos = bondonenode2->getPosition();

    const Eigen::Vector3d &bondtwonode1Pos = bondtwonode1->getPosition();
    const Eigen::Vector3d &bondtwonode2Pos = bondtwonode2->getPosition();
    double s = 0, t = 0;

    SphereCylinderBond::computeOverlap(bondonenode1Pos, bondonenode2Pos, s, bondtwonode1Pos, bondtwonode2Pos, t);
    Eigen::Vector3d projection1 = (1 - s) * bondonenode1Pos + s * bondonenode2Pos;
    Eigen::Vector3d projection2 = (1 - t) * bondtwonode1Pos + t * bondtwonode2Pos;
    double distance = (projection1 - projection2).norm();
    double fixeddis = fiber1Radius + fiber2Radius;
    double normal_overlap = fixeddis - distance;

    if (normal_overlap > 0)
    {
        double weights = normal_overlap * s;
        double weightt = normal_overlap * t;
        double weight = normal_overlap;
        Eigen::Vector3d directionone = bondonenode2Pos - bondonenode1Pos;
        Eigen::Vector3d directiontwo = bondtwonode2Pos - bondtwonode1Pos;
        double L1 = directionone.squaredNorm();
        double L2 = directiontwo.squaredNorm();

        double s1 = (bondtwonode1Pos - bondonenode1Pos).dot(directionone) / L1;
        if (s1 >= 0 && s1 <= 1 && s1 != s)
        {
            double overlap1 = fixeddis - (bondtwonode1Pos - (bondonenode1Pos + s1 * directionone)).norm();
            if (overlap1 > 0)
            {
                weights += overlap1 * s1;
                weight += overlap1;
            }
        }

        double s2 = (bondtwonode2Pos - bondonenode1Pos).dot(directionone) / L1;
        if (s2 >= 0 && s2 <= 1 && s2 != s)
        {
            double overlap2 = fixeddis - (bondtwonode2Pos - (bondonenode1Pos + s2 * directionone)).norm();
            if (overlap2 > 0)
            {
                weights += overlap2 * s2;
                weightt += overlap2;
                weight += overlap2;
            }
        }

        double t1 = (bondonenode1Pos - bondtwonode1Pos).dot(directiontwo) / L2;
        if (t1 > 0 && t1 < 1 && t1 != t)
        {
            double overlap3 = fixeddis - (bondonenode1Pos - (bondtwonode1Pos + t1 * directiontwo)).norm();
            if (overlap3 > 0)
            {
                weightt += overlap3 * t1;
                weight += overlap3;
            }
        }

        double t2 = (bondonenode2Pos - bondtwonode1Pos).dot(directiontwo) / L2;
        if (t2 > 0 && t2 < 1 && t2 != t)
        {
            double overlap4 = fixeddis - (bondonenode2Pos - (bondtwonode1Pos + t2 * directiontwo)).norm();
            if (overlap4 > 0)
            {
                weights += overlap4;
                weightt += overlap4 * t2;
                weight += overlap4;
            }
        }

        weights /= weight;
        weightt /= weight;

        // contact point
        projection1 = (1 - weights) * bondonenode1Pos + weights * bondonenode2Pos;
        projection2 = (1 - weightt) * bondtwonode1Pos + weightt * bondtwonode2Pos;

        const Eigen::Vector3d &bondonenode1Velocity = bondonenode1->getVelocity();
        const Eigen::Vector3d &bondonenode2Velocity = bondonenode2->getVelocity();

        const Eigen::Vector3d &bondtwonode1Velocity = bondtwonode1->getVelocity();
        const Eigen::Vector3d &bondtwonode2Velocity = bondtwonode2->getVelocity();

        const Eigen::Vector3d &bondonenode1Omega = bondonenode1->getOmega();
        const Eigen::Vector3d &bondonenode2Omega = bondonenode2->getOmega();

        const Eigen::Vector3d &bondtwonode1Omega = bondtwonode1->getOmega();
        const Eigen::Vector3d &bondtwonode2Omega = bondtwonode2->getOmega();

        Eigen::Vector3d v1 = weights * (bondonenode2Velocity - bondonenode1Velocity) + bondonenode1Velocity;

        Eigen::Vector3d w1 = weights * (bondonenode2Omega - bondonenode1Omega) + bondonenode1Omega;

        Eigen::Vector3d v2 = weightt * (bondtwonode2Velocity - bondtwonode1Velocity) + bondtwonode1Velocity;

        Eigen::Vector3d w2 = weightt * (bondtwonode2Omega - bondtwonode1Omega) + bondtwonode1Omega;

        Eigen::Vector3d normal_unit_vector = projection2 - projection1;
        normal_unit_vector.normalize();
        // Calculation of contact relative velocity
        Eigen::Vector3d contact_relative_velocity = (v1 - v2) + (fiber1Radius * w1 + fiber2Radius * w2).cross(normal_unit_vector);

        double normal_relative_velocity_value = contact_relative_velocity.dot(normal_unit_vector);

        // Calculation of normal relative velocity
        Eigen::Vector3d normal_relative_velocity = normal_relative_velocity_value * normal_unit_vector;
        // Calculation of tangential relative velocity
        Eigen::Vector3d tangential_relative_velocity = contact_relative_velocity - normal_relative_velocity;
        // modify tangential_overlap
        int id1 = fiberbond1->getId();
        int id2 = fiberbond2->getId();

        auto &innerMap = fiberfibercontactinformationlist[id1]; // Access or create the inner map
        auto innerIt = innerMap.find(id2);

        if (innerIt == innerMap.end())
        {
            // Contact information for this pair does not exist, create a new one
            ContactInformation newContactInfo(id1, id2);
            innerMap[id2] = newContactInfo; // Insert the new contact information
        }

        ContactInformation &contact_info = fiberfibercontactinformationlist[id1][id2];

        Eigen::Vector3d modified_tangential_displacement = contact_info.tangentialDisplacement + contact_info.previousTangentialVelocity * timeStep;
        // Updating the contact_info container based on the new calculated values
        contact_info.tangentialDisplacement = modified_tangential_displacement;
        contact_info.previousTangentialVelocity = tangential_relative_velocity;

        double radius_times_overlap_sqrt = sqrt(effectiveradius[particletype_one][particletype_two] * normal_overlap);

        double sn = 2 * effectiveyoungsmodulus[particletype_one][particletype_two] * radius_times_overlap_sqrt;

        double normal_contact_stiffness = 0.666666 * sn;

        double normal_damping_constant = 1.8257 * modelparameterbeta[particletype_one][particletype_two] * sqrt(sn * effectivemass[particletype_one][particletype_two]);

        // tangential stiffness
        double tangential_contact_stiffness = 8 * effectiveshearmodulus[particletype_one][particletype_two] * radius_times_overlap_sqrt;

        double tangential_damping_constant = 1.8257 * modelparameterbeta[particletype_one][particletype_two] * sqrt(tangential_contact_stiffness * effectivemass[particletype_one][particletype_two]);

        Eigen::Vector3d normal_damping_force = (normal_damping_constant * normal_relative_velocity_value) * normal_unit_vector;
        // Calculation of normal force using spring and dashpot normal forces
        Eigen::Vector3d normal_force =
            ((normal_contact_stiffness * normal_overlap) * normal_unit_vector) + normal_damping_force;

        energyDissipation += std::fabs(normal_damping_force.dot(normal_relative_velocity * timeStep));

        // Calculation of tangential force using spring and dashpot tangential
        // forces. Since we need dashpot tangential force in the gross sliding again,
        // we define it as a separate variable
        Eigen::Vector3d dashpot_tangential_force = tangential_damping_constant * tangential_relative_velocity;
        Eigen::Vector3d tangential_force = (tangential_contact_stiffness * modified_tangential_displacement) + dashpot_tangential_force;

        double coulomb_threshold = effectiveslidingfriction[particletype_one][particletype_two] * normal_force.norm();

        double tangential_force_value = tangential_force.norm();
        // Check for gross sliding
        if (tangential_force_value > coulomb_threshold)
        {
            // Gross sliding occurs and the tangential overlap and tangnetial
            // force are limited to Coulumb's criterion
            Eigen::Vector3d tangential_overlap = (coulomb_threshold * (tangential_force / (tangential_force_value + DBL_MIN)) - dashpot_tangential_force) /
                                                 (tangential_contact_stiffness + DBL_MIN);

            contact_info.tangentialDisplacement = tangential_overlap;
            tangential_force = (tangential_contact_stiffness * tangential_overlap) + dashpot_tangential_force;
        }
        energyDissipation += std::fabs(tangential_force.dot(tangential_relative_velocity * timeStep));

        // Calculation of torque
        // Torque caused by tangential force (tangential_torque)
        Eigen::Vector3d tangential_torque1 = (fiber1Radius * normal_unit_vector).cross(tangential_force);
        Eigen::Vector3d tangential_torque2 = (fiber2Radius * normal_unit_vector).cross(tangential_force);

        Eigen::Vector3d results_force = normal_force + tangential_force;

        bondonenode1->addForce((weights - 1) * results_force);
        bondonenode2->addForce(-weights * results_force);

        bondtwonode1->addForce((1 - weightt) * results_force);
        bondtwonode2->addForce(weightt * results_force);

        bondonenode1->addTorque((weights - 1) * tangential_torque1);
        bondonenode2->addTorque(-weights * tangential_torque1);

        bondtwonode1->addTorque((weightt - 1) * tangential_torque2);
        bondtwonode2->addTorque(-weightt * tangential_torque2);

        contact_info.normalForce = -normal_force;
        contact_info.tangentialForce = -tangential_force;
        nextfiberfibercontactinformationlist[id1][id2] = contact_info;
    }
}
void ContactForce::computePlaneWallFiberForce(const std::unique_ptr<PlaneWall> &planewall, const std::unique_ptr<SphereCylinderBond> &fiberbond,
                                              const std::unique_ptr<SphereParticle> &node1, const std::unique_ptr<SphereParticle> &node2, double timeStep)
{
    Eigen::Vector3d sphereCenter1 = node1->getPosition();
    Eigen::Vector3d sphereCenter2 = node2->getPosition();

    auto &manager = fiberbond->getParticlePropertyManager();

    PropertyTypeID fibertype = fiberbond->getType();
    PropertyTypeID walltype = planewall->getType();

    auto &fiberproperty = manager->getFiberProperties(fibertype);
    double fiberRadius = fiberproperty->getRadius();

    // Retrieve the plane wall's point and normal
    Eigen::Vector3d planePoint = planewall->getCorner1();

    Eigen::Vector3d planeNormal = planewall->getNormal();

    // Compute the vector from a point on the plane to the sphere's center
    Eigen::Vector3d vecToSphereCenter1 = sphereCenter1 - planePoint;
    // Compute the distance from the sphere's center to the plane
    double distanceToPlane1 = vecToSphereCenter1.dot(planeNormal);

    if (distanceToPlane1 < 0)
    {
        planeNormal = -planeNormal;
        distanceToPlane1 = -distanceToPlane1;
    }
    Eigen::Vector3d vecToSphereCenter2 = sphereCenter2 - planePoint;
    double distanceToPlane2 = vecToSphereCenter2.dot(planeNormal);

    double weight_s = 0;
    double weight = 0;

    double normal_overlap1 = fiberRadius - distanceToPlane1;
    if (normal_overlap1 > 0)
    {
        weight = normal_overlap1;
    }
    double normal_overlap2 = fiberRadius - distanceToPlane2;
    if (normal_overlap2 > 0)
    {
        weight_s += normal_overlap2;
        weight += normal_overlap2;
    }
    double normal_overlap = normal_overlap1 > normal_overlap2 ? normal_overlap1 : normal_overlap2;

    if (normal_overlap > 0)
    {
        weight_s = weight_s / weight;
        // Using contact_vector, the contact normal vector is obtained
        Eigen::Vector3d normal_unit_vector = -planeNormal;

        const Eigen::Vector3d &node1Velocity = node1->getVelocity();
        const Eigen::Vector3d &node2Velocity = node2->getVelocity();

        const Eigen::Vector3d &node1Omega = node1->getOmega();
        const Eigen::Vector3d &node2Omega = node2->getOmega();

        Eigen::Vector3d fiber_velocity = weight_s * (node2Velocity - node1Velocity) + node1Velocity;
        Eigen::Vector3d fiber_omega = weight_s * (node2Omega - node1Omega) + node1Omega;
        // Calculation of contact relative velocity
        Eigen::Vector3d contact_relative_velocity = fiber_velocity + (fiberRadius * fiber_omega).cross(normal_unit_vector);

        double normal_relative_velocity_value = contact_relative_velocity.dot(normal_unit_vector);
        // Calculation of normal relative velocity
        Eigen::Vector3d normal_relative_velocity = normal_relative_velocity_value * normal_unit_vector;
        // Calculation of tangential relative velocity
        Eigen::Vector3d tangential_relative_velocity = contact_relative_velocity - normal_relative_velocity;
        // modify tangential_overlap
        int wallId = planewall->getId();
        int fiberId = fiberbond->getId();

        auto &innerMap = planewallfibercontactinformationlist[wallId]; // Access or create the inner map
        auto innerIt = innerMap.find(fiberId);

        if (innerIt == innerMap.end())
        {
            // Contact information for this pair does not exist, create a new one
            ContactInformation newContactInfo(wallId, fiberId);
            innerMap[fiberId] = newContactInfo; // Insert the new contact information
        }

        ContactInformation &contact_info = planewallfibercontactinformationlist[wallId][fiberId];

        Eigen::Vector3d modified_tangential_displacement = contact_info.tangentialDisplacement + contact_info.previousTangentialVelocity * timeStep;
        // Updating the contact_info container based on the new calculated values

        contact_info.tangentialDisplacement = modified_tangential_displacement;
        contact_info.previousTangentialVelocity = tangential_relative_velocity;

        double radius_times_overlap_sqrt = sqrt(fiberRadius * normal_overlap);

        double sn = 2 * effectiveyoungsmodulus[walltype][fibertype] * radius_times_overlap_sqrt;

        double normal_contact_stiffness = 0.666666 * sn;

        double normal_damping_constant = 1.8257 * modelparameterbeta[walltype][fibertype] * sqrt(sn * fiberproperty->getNodemass());

        // tangential stiffness
        double tangential_contact_stiffness = 8 * effectiveshearmodulus[walltype][fibertype] * radius_times_overlap_sqrt;

        double tangential_damping_constant = 1.8257 * modelparameterbeta[walltype][fibertype] * sqrt(tangential_contact_stiffness * fiberproperty->getNodemass());

        Eigen::Vector3d normal_damping_force = (normal_damping_constant * normal_relative_velocity_value) * normal_unit_vector;
        // Calculation of normal force using spring and dashpot normal forces
        Eigen::Vector3d normal_force =
            ((normal_contact_stiffness * normal_overlap) * normal_unit_vector) + normal_damping_force;

        energyDissipation += std::fabs(normal_damping_force.dot(normal_relative_velocity * timeStep));

        // Calculation of tangential force using spring and dashpot tangential
        // forces. Since we need dashpot tangential force in the gross sliding again,
        // we define it as a separate variable
        Eigen::Vector3d dashpot_tangential_force = tangential_damping_constant * tangential_relative_velocity;
        Eigen::Vector3d tangential_force = (tangential_contact_stiffness * modified_tangential_displacement) + dashpot_tangential_force;

        double coulomb_threshold = effectiveslidingfriction[walltype][fibertype] * normal_force.norm();

        double tangential_force_value = tangential_force.norm();
        // Check for gross sliding
        if (tangential_force_value > coulomb_threshold)
        {
            // Gross sliding occurs and the tangential overlap and tangnetial
            // force are limited to Coulumb's criterion
            Eigen::Vector3d tangential_overlap = (coulomb_threshold * (tangential_force / (tangential_force_value + DBL_MIN)) - dashpot_tangential_force) /
                                                 (tangential_contact_stiffness + DBL_MIN);

            contact_info.tangentialDisplacement = tangential_overlap;
            tangential_force = (tangential_contact_stiffness * tangential_overlap) + dashpot_tangential_force;
        }
        energyDissipation += std::fabs(tangential_force.dot(tangential_relative_velocity * timeStep));

        // Calculation of torque
        // Torque caused by tangential force (tangential_torque)
        Eigen::Vector3d tangential_torque1 = (fiberRadius * normal_unit_vector).cross(tangential_force);

        Eigen::Vector3d results_force = normal_force + tangential_force;

        node1->addForce((weight_s - 1) * results_force);
        node2->addForce(-weight_s * results_force);
        node1->addTorque((weight_s - 1) * tangential_torque1);
        node2->addTorque(-weight_s * tangential_torque1);
        planewall->addForce(results_force);

        contact_info.normalForce = normal_force;
        contact_info.tangentialForce = tangential_force;
        nextplanewallfibercontactinformationlist[wallId][fiberId] = contact_info;
    }
}
void ContactForce::computeCylinderWallFiberForce(const std::unique_ptr<CylinderContainer> &cylinderWall, const std::unique_ptr<SphereCylinderBond> &fiberbond,
                                                 const std::unique_ptr<SphereParticle> &node1, const std::unique_ptr<SphereParticle> &node2, double timeStep)
{
    Eigen::Vector3d sphereCenter1 = node1->getPosition();
    Eigen::Vector3d sphereCenter2 = node2->getPosition();

    const auto &manager = fiberbond->getParticlePropertyManager();

    PropertyTypeID fibertype = fiberbond->getType();
    PropertyTypeID walltype = cylinderWall->getType();

    auto &fiberproperty = manager->getFiberProperties(fibertype);

    double fiberRadius = fiberproperty->getRadius();
    double cylinderRadius = cylinderWall->getRadius();
    // Retrieve the cylinder wall's point
    Eigen::Vector3d cylinderEnd = cylinderWall->getEndpoint();

    Eigen::Vector3d cylinderStart = cylinderWall->getStartpoint();

    // Calculate the axis direction of the cylinder
    Eigen::Vector3d cylinderAxis = (cylinderEnd - cylinderStart).normalized();

    // Calculate the vector from the start point of the cylinder to the sphere center
    Eigen::Vector3d startToCenter1 = sphereCenter1 - cylinderStart;

    // Project the vector onto the cylinder axis to find the closest point on the axis to the sphere center
    double projection1 = startToCenter1.dot(cylinderAxis);
    projection1 = std::clamp(projection1, 0.0, 1.0);
    Eigen::Vector3d closestPoint1 = cylinderStart + projection1 * cylinderAxis;
    // Eigen::Vector3d normal_unit_vector = sphereCenter1 - closestPoint;
    double distanceToAxis1 = (sphereCenter1 - closestPoint1).norm();

    Eigen::Vector3d startToCenter2 = sphereCenter2 - cylinderStart;

    // Project the vector onto the cylinder axis to find the closest point on the axis to the sphere center
    double projection2 = startToCenter2.dot(cylinderAxis);
    projection2 = std::clamp(projection2, 0.0, 1.0);
    Eigen::Vector3d closestPoint2 = cylinderStart + projection2 * cylinderAxis;
    double distanceToAxis2 = (sphereCenter2 - closestPoint2).norm();

    double normal_overlap1 = 0;
    double normal_overlap2 = 0;
    double weight_s = 0;
    double weight_t = 0;
    double weight = 0;

    if (distanceToAxis1 > cylinderRadius)
    {
        normal_overlap1 = fiberRadius + cylinderRadius - distanceToAxis1;
        normal_overlap2 = fiberRadius + cylinderRadius - distanceToAxis2;
    }
    else
    {
        normal_overlap1 = distanceToAxis1 - (cylinderRadius - fiberRadius);
        normal_overlap2 = distanceToAxis2 - (cylinderRadius - fiberRadius);
    }
    if (normal_overlap1 > 0)
    {
        weight = normal_overlap1;
        weight_t = projection1 * normal_overlap1;
    }
    if (normal_overlap2 > 0)
    {
        weight_s += normal_overlap2;
        weight += normal_overlap2;
        weight_t += projection2 * normal_overlap2;
    }

    double normal_overlap = normal_overlap1 > normal_overlap2 ? normal_overlap1 : normal_overlap2;

    if (normal_overlap > 0)
    {
        weight_s = weight_s / weight;
        weight_t = weight_t / weight;
        Eigen::Vector3d closestPoint = cylinderStart + weight_t * cylinderAxis;
        Eigen::Vector3d sphereCenter = (1 - weight_s) * sphereCenter1 + weight_s * sphereCenter2;
        Eigen::Vector3d normal_unit_vector = (sphereCenter - closestPoint).normalized();
        if (distanceToAxis1 > cylinderRadius)
        {
            normal_unit_vector = -normal_unit_vector;
        }
        const Eigen::Vector3d &node1Velocity = node1->getVelocity();
        const Eigen::Vector3d &node2Velocity = node2->getVelocity();

        const Eigen::Vector3d &node1Omega = node1->getOmega();
        const Eigen::Vector3d &node2Omega = node2->getOmega();

        Eigen::Vector3d fiber_velocity = weight_s * (node2Velocity - node1Velocity) + node1Velocity;
        Eigen::Vector3d fiber_omega = weight_s * (node2Omega - node1Omega) + node1Omega;
        // Calculation of contact relative velocity
        Eigen::Vector3d contact_relative_velocity = fiber_velocity + (fiberRadius * fiber_omega).cross(normal_unit_vector);

        double normal_relative_velocity_value = contact_relative_velocity.dot(normal_unit_vector);
        // Calculation of normal relative velocity
        Eigen::Vector3d normal_relative_velocity = normal_relative_velocity_value * normal_unit_vector;
        // Calculation of tangential relative velocity
        Eigen::Vector3d tangential_relative_velocity = contact_relative_velocity - normal_relative_velocity;
        // modify tangential_overlap
        int wallId = cylinderWall->getId();
        int fiberId = fiberbond->getId();

        auto &innerMap = cylinderwallfibercontactinformationlist[wallId]; // Access or create the inner map
        auto innerIt = innerMap.find(fiberId);

        if (innerIt == innerMap.end())
        {
            // Contact information for this pair does not exist, create a new one
            ContactInformation newContactInfo(wallId, fiberId);
            innerMap[fiberId] = newContactInfo; // Insert the new contact information
        }

        ContactInformation &contact_info = cylinderwallfibercontactinformationlist[wallId][fiberId];

        Eigen::Vector3d modified_tangential_displacement = contact_info.tangentialDisplacement + contact_info.previousTangentialVelocity * timeStep;
        // Updating the contact_info container based on the new calculated values

        contact_info.tangentialDisplacement = modified_tangential_displacement;
        contact_info.previousTangentialVelocity = tangential_relative_velocity;

        double radius_times_overlap_sqrt = sqrt(fiberRadius * normal_overlap);

        double sn = 2 * effectiveyoungsmodulus[walltype][fibertype] * radius_times_overlap_sqrt;

        double normal_contact_stiffness = 0.666666 * sn;

        double normal_damping_constant = 1.8257 * modelparameterbeta[walltype][fibertype] * sqrt(sn * fiberproperty->getNodemass());

        // tangential stiffness
        double tangential_contact_stiffness = 8 * effectiveshearmodulus[walltype][fibertype] * radius_times_overlap_sqrt;

        double tangential_damping_constant = 1.8257 * modelparameterbeta[walltype][fibertype] * sqrt(tangential_contact_stiffness * fiberproperty->getNodemass());

        Eigen::Vector3d normal_damping_force = (normal_damping_constant * normal_relative_velocity_value) * normal_unit_vector;
        // Calculation of normal force using spring and dashpot normal forces
        Eigen::Vector3d normal_force =
            ((normal_contact_stiffness * normal_overlap) * normal_unit_vector) + normal_damping_force;

        energyDissipation += std::fabs(normal_damping_force.dot(normal_relative_velocity * timeStep));

        // Calculation of tangential force using spring and dashpot tangential
        // forces. Since we need dashpot tangential force in the gross sliding again,
        // we define it as a separate variable
        Eigen::Vector3d dashpot_tangential_force = tangential_damping_constant * tangential_relative_velocity;
        Eigen::Vector3d tangential_force = (tangential_contact_stiffness * modified_tangential_displacement) + dashpot_tangential_force;

        double coulomb_threshold = effectiveslidingfriction[walltype][fibertype] * normal_force.norm();

        double tangential_force_value = tangential_force.norm();
        // Check for gross sliding
        if (tangential_force_value > coulomb_threshold)
        {
            // Gross sliding occurs and the tangential overlap and tangnetial
            // force are limited to Coulumb's criterion
            Eigen::Vector3d tangential_overlap = (coulomb_threshold * (tangential_force / (tangential_force_value + DBL_MIN)) - dashpot_tangential_force) /
                                                 (tangential_contact_stiffness + DBL_MIN);

            contact_info.tangentialDisplacement = tangential_overlap;
            tangential_force = (tangential_contact_stiffness * tangential_overlap) + dashpot_tangential_force;
        }
        energyDissipation += std::fabs(tangential_force.dot(tangential_relative_velocity * timeStep));

        // Calculation of torque
        // Torque caused by tangential force (tangential_torque)
        Eigen::Vector3d tangential_torque1 = (fiberRadius * normal_unit_vector).cross(tangential_force);

        Eigen::Vector3d results_force = normal_force + tangential_force;

        node1->addForce((weight_s - 1) * results_force);
        node2->addForce(-weight_s * results_force);
        node1->addTorque((weight_s - 1) * tangential_torque1);
        node2->addTorque(-weight_s * tangential_torque1);
        cylinderWall->addForce(results_force);

        contact_info.normalForce = normal_force;
        contact_info.tangentialForce = tangential_force;
        nextcylinderwallfibercontactinformationlist[wallId][fiberId] = contact_info;
    }
}
void ContactForce::updateContactInformation()
{
    spherespherecontactinformationlist.swap(nextspherespherecontactinformationlist);
    nextspherespherecontactinformationlist.clear();

    spherefibercontactinformationlist.swap(nextspherefibercontactinformationlist);
    nextspherefibercontactinformationlist.clear();

    planewallspherecontactinformationlist.swap(nextplanewallspherecontactinformationlist);
    nextplanewallspherecontactinformationlist.clear();

    planewallfibercontactinformationlist.swap(nextplanewallfibercontactinformationlist);
    nextplanewallfibercontactinformationlist.clear();

    fiberfibercontactinformationlist.swap(nextfiberfibercontactinformationlist);
    nextfiberfibercontactinformationlist.clear();

    cylinderwallspherecontactinformationlist.swap(nextcylinderwallspherecontactinformationlist);
    nextcylinderwallspherecontactinformationlist.clear();

    cylinderwallfibercontactinformationlist.swap(nextcylinderwallfibercontactinformationlist);
    nextcylinderwallfibercontactinformationlist.clear();
}
std::unordered_map<int, std::unordered_map<int, ContactInformation>> &ContactForce::getSphereSphereContactInformationList()
{
    return spherespherecontactinformationlist;
}
std::unordered_map<int, std::unordered_map<int, ContactInformation>> &ContactForce::getplanewallSphereContactInformationList()
{
    return planewallspherecontactinformationlist;
}
std::unordered_map<int, std::unordered_map<int, ContactInformation>> &ContactForce::getplanewallFiberContactInformationList()
{
    return planewallfibercontactinformationlist;
}
std::unordered_map<int, std::unordered_map<int, ContactInformation>> &ContactForce::getSphereFiberContactInformationList()
{
    return spherefibercontactinformationlist;
}
std::unordered_map<int, std::unordered_map<int, ContactInformation>> &ContactForce::getFiberFiberContactInformationList()
{
    return fiberfibercontactinformationlist;
}
std::unordered_map<int, std::unordered_map<int, ContactInformation>> &ContactForce::getcylinderwallSphereContactInformationList()
{
    return cylinderwallspherecontactinformationlist;
}
std::unordered_map<int, std::unordered_map<int, ContactInformation>> &ContactForce::getcylinderwallFiberContactInformationList()
{
    return cylinderwallfibercontactinformationlist;
}