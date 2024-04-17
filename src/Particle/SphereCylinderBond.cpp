#include "SphereCylinderBond.h"
SphereCylinderBond::SphereCylinderBond()
{
}
SphereCylinderBond::SphereCylinderBond(int id, PropertyTypeID type,
                                       const std::shared_ptr<ParticlePropertyManager> &manager,
                                       int fiberid,
                                       int node1, int node2, int neighborelement1, int neighborelement2,
                                       double energyDissipation,
                                       const Eigen::Vector3d &tangentialforce,
                                       const Eigen::Vector3d &twisttorque,
                                       const Eigen::Vector3d &bendtorque)
    : id(id), type(type), manager(manager), fiberid(fiberid), node1(node1), node2(node2),
      neighborelement1(neighborelement1), neighborelement2(neighborelement2), energyDissipation(energyDissipation),
      tangentialforce(tangentialforce), twisttorque(twisttorque), bendtorque(bendtorque)
{
}
int SphereCylinderBond::getNode1() const
{
    return node1;
}
int SphereCylinderBond::getNode2() const
{
    return node2;
}
int SphereCylinderBond::getNeighborelement1() const
{
    return neighborelement1;
}
int SphereCylinderBond::getNeighborelement2() const
{
    return neighborelement2;
}
int SphereCylinderBond::getId() const
{
    return this->id;
}
int SphereCylinderBond::getFiberId() const
{
    return this->fiberid;
}

const PropertyTypeID &SphereCylinderBond::getType() const
{
    return type;
}
void SphereCylinderBond::updateBond(std::unique_ptr<SphereParticle> &sphere1, std::unique_ptr<SphereParticle> &sphere2, double timeStep)
{
    const auto fiberproperties = manager->getFiberProperties(type);
    const Eigen::Vector3d &position1 = sphere1->getPosition();
    const Eigen::Vector3d &position2 = sphere2->getPosition();

    const Eigen::Vector3d &velocity1 = sphere1->getVelocity();
    const Eigen::Vector3d &velocity2 = sphere2->getVelocity();

    const Eigen::Vector3d &omega1 = sphere1->getOmega();
    const Eigen::Vector3d &omega2 = sphere2->getOmega();

    Eigen::Vector3d normal_length = position2 - position1;
    double length = normal_length.norm();
    Eigen::Vector3d normal = normal_length / length;

    Eigen::Vector3d normal_displacement = (length - fiberproperties->getElementlength()) * normal;

    Eigen::Vector3d r1 = 0.5 * (position1 + position2) - position1;
    Eigen::Vector3d r2 = -r1;

    Eigen::Vector3d contactvelocity1 = velocity1 + omega1.cross(r1);
    Eigen::Vector3d contactvelocity2 = velocity2 + omega2.cross(r2);

    Eigen::Vector3d relative_velocity = contactvelocity2 - contactvelocity1;

    Eigen::Vector3d normal_velocity = relative_velocity.dot(normal) * normal;

    Eigen::Vector3d tangential_velocity = relative_velocity - normal_velocity;

    Eigen::Vector3d normal_bond_force = fiberproperties->getNormalstiffnesses() * normal_displacement;
    Eigen::Vector3d normal_bond_damping_force = fiberproperties->getBonddampingcoefficient() *
                                                sqrt(2 * fiberproperties->getNodemass() * fiberproperties->getNormalstiffnesses()) * normal_velocity;
    energyDissipation += std::fabs(normal_bond_damping_force.dot(normal_velocity * timeStep));
    Eigen::Vector3d tangential_bond_damping_force = fiberproperties->getBonddampingcoefficient() *
                                                    sqrt(2 * fiberproperties->getNodemass() * fiberproperties->getShearstiffnesses()) * tangential_velocity;
    energyDissipation += std::fabs(tangential_bond_damping_force.dot(tangential_velocity * timeStep));

    tangentialforce += fiberproperties->getShearstiffnesses() * tangential_velocity * timeStep;

    Eigen::Vector3d tangential_torque1 = r1.cross(tangential_bond_damping_force + tangentialforce);
    // Eigen::Vector3d tangential_torque2 = r2.cross (tangential_bond_damping_force + tangentialforce);

    Eigen::Vector3d resultFoce = normal_bond_force + normal_bond_damping_force + tangentialforce + tangential_bond_damping_force;

    sphere1->addForce(resultFoce); // note the reference is sphere1, and these force is resist the motion
    sphere2->addForce(-resultFoce);

    sphere1->addTorque(tangential_torque1);
    sphere2->addTorque(tangential_torque1);

    Eigen::Vector3d relative_angle_velocity = omega2 - omega1;
    Eigen::Vector3d relative_twist_velocity = relative_angle_velocity.dot(normal) * normal;
    Eigen::Vector3d relative_bend_velocity = relative_angle_velocity - relative_twist_velocity;

    twisttorque += fiberproperties->getTwiststiffnesses() * relative_twist_velocity * timeStep;
    bendtorque += fiberproperties->getBendingstiffnesses() * relative_bend_velocity * timeStep;

    Eigen::Vector3d twist_bond_damping_torque = fiberproperties->getBonddampingcoefficient() *
                                                sqrt(2 * fiberproperties->getNodemomentofinertia() * fiberproperties->getTwiststiffnesses()) * relative_twist_velocity;
    energyDissipation += std::fabs(twist_bond_damping_torque.dot(relative_twist_velocity * timeStep));

    Eigen::Vector3d bend_bond_damping_torque = fiberproperties->getBonddampingcoefficient() *
                                               sqrt(2 * fiberproperties->getNodemomentofinertia() * fiberproperties->getBendingstiffnesses()) * relative_bend_velocity;
    energyDissipation += std::fabs(bend_bond_damping_torque.dot(relative_bend_velocity * timeStep));
    Eigen::Vector3d resultTorqur = twisttorque + bendtorque + twist_bond_damping_torque + bend_bond_damping_torque;
    sphere1->addTorque(resultTorqur);
    sphere2->addTorque(-resultTorqur);
}

void SphereCylinderBond::computeOverlap(const Eigen::Vector3d &node1position, const Eigen::Vector3d &node2position,
                                        const Eigen::Vector3d &sphereposition, double &t, Eigen::Vector3d &projection)
{
    Eigen::Vector3d bondVector = node2position - node1position;           // Vector along the bond
    Eigen::Vector3d sphereToNode1Vector = sphereposition - node1position; // Vector from node1 to sphere center

    // Projection of sphereToNode1Vector onto bondVector
    t = sphereToNode1Vector.dot(bondVector) / bondVector.squaredNorm();

    // Clamp t to the segment
    t = std::max(0.0, std::min(1.0, t));

    // Projection point on the bond segment
    projection = node1position + t * bondVector;
}
void SphereCylinderBond::computeOverlap(const Eigen::Vector3d &thisnode1position, const Eigen::Vector3d &thisnode2position, double &s,
                                        const Eigen::Vector3d &anothernode1position, const Eigen::Vector3d &anothernode2position, double &t)
{
    Eigen::Vector3d u = thisnode2position - thisnode1position;
    Eigen::Vector3d v = anothernode2position - anothernode1position;
    Eigen::Vector3d w = thisnode1position - anothernode1position;

    double a = u.dot(u); // Always >= 0
    double b = u.dot(v);
    double c = v.dot(v); // Always >= 0
    double d = u.dot(w);
    double e = v.dot(w);

    // The derivatives dR/ds(i,j) at the four corners of the domain.
    double f00 = d;
    double f10 = f00 + a;
    double f01 = f00 - b;
    double f11 = f10 - b;

    // The derivatives dR/dt(i,j) at the four corners of the domain.
    double g00 = -e;
    double g10 = g00 - b;
    double g01 = g00 + c;
    double g11 = g10 + c;

    double const zero = static_cast<double>(0);
    double const one = static_cast<double>(1);
    if (a > zero && c > zero)
    {
        // Compute the solutions to dR/ds(s0,0) = 0 and
        // dR/ds(s1,1) = 0.  The location of sI on the s-axis is
        // stored in classifyI (I = 0 or 1).  If sI <= 0, classifyI
        // is -1.  If sI >= 1, classifyI is 1.  If 0 < sI < 1,
        // classifyI is 0.  This information helps determine where to
        // search for the minimum point (s,t).  The fij values are
        // dR/ds(i,j) for i and j in {0,1}.

        std::array<double, 2> sValue{
            GetClampedRoot(a, f00, f10),
            GetClampedRoot(a, f01, f11)};

        std::array<int32_t, 2> classify{};
        for (size_t i = 0; i < 2; ++i)
        {
            if (sValue[i] <= zero)
            {
                classify[i] = -1;
            }
            else if (sValue[i] >= one)
            {
                classify[i] = +1;
            }
            else
            {
                classify[i] = 0;
            }
        }

        if (classify[0] == -1 && classify[1] == -1)
        {
            // The minimum must occur on s = 0 for 0 <= t <= 1.
            s = zero;
            t = GetClampedRoot(c, g00, g01);
        }
        else if (classify[0] == +1 && classify[1] == +1)
        {
            // The minimum must occur on s = 1 for 0 <= t <= 1.
            s = one;
            t = GetClampedRoot(c, g10, g11);
        }
        else
        {
            // The line dR/ds = 0 intersects the domain [0,1]^2 in a
            // nondegenerate segment. Compute the endpoints of that
            // segment, end[0] and end[1]. The edge[i] flag tells you
            // on which domain edge end[i] lives: 0 (s=0), 1 (s=1),
            // 2 (t=0), 3 (t=1).
            std::array<int32_t, 2> edge{0, 0};
            std::array<std::array<double, 2>, 2> end{};
            ComputeIntersection(sValue, classify, b, f00, f10, edge, end);

            // The directional derivative of R along the segment of
            // intersection is
            //   H(z) = (end[1][1]-end[1][0]) *
            //          dR/dt((1-z)*end[0] + z*end[1])
            // for z in [0,1]. The formula uses the fact that
            // dR/ds = 0 on the segment. Compute the minimum of
            // H on [0,1].
            ComputeMinimumParameters(edge, end, b, c, e, g00, g10,
                                     g01, g11, s, t);
        }
    }
    else
    {
        if (a > zero)
        {
            // The Q-segment is degenerate (Q0 and Q1 are the same
            // point) and the quadratic is R(s,0) = a*s^2 + 2*d*s + f
            // and has (half) first derivative F(t) = a*s + d.  The
            // closest P-point is interior to the P-segment when
            // F(0) < 0 and F(1) > 0.
            s = GetClampedRoot(a, f00, f10);
            t = zero;
        }
        else if (c > zero)
        {
            // The P-segment is degenerate (P0 and P1 are the same
            // point) and the quadratic is R(0,t) = c*t^2 - 2*e*t + f
            // and has (half) first derivative G(t) = c*t - e.  The
            // closest Q-point is interior to the Q-segment when
            // G(0) < 0 and G(1) > 0.
            s = zero;
            t = GetClampedRoot(c, g00, g01);
        }
        else
        {
            // P-segment and Q-segment are degenerate.
            s = zero;
            t = zero;
        }
    }
    // projection1 = (1 - s) * thisnode1position + s * thisnode2position;
    // projection2 = (1 - t) * anothernode1position + t * anothernode2position;
}

void SphereCylinderBond::computeOverlap(const Eigen::Vector3d &node1position, const Eigen::Vector3d &node2position,
                                        const std::unique_ptr<PlaneWall> &planewall, double &distance)
{

    // Retrieve the plane wall's point and normal
    Eigen::Vector3d planePoint = planewall->getCorner1();

    Eigen::Vector3d planeNormal = planewall->getNormal();

    // Compute the vector from a point on the plane to the sphere's center
    Eigen::Vector3d vecToSphereCenter1 = node1position - planePoint;
    // Compute the distance from the sphere's center to the plane
    double distanceToPlane1 = vecToSphereCenter1.dot(planeNormal);

    if (distanceToPlane1 < 0)
    {
        planeNormal = -planeNormal;
        distanceToPlane1 = -distanceToPlane1;
    }
    Eigen::Vector3d vecToSphereCenter2 = node2position - planePoint;
    double distanceToPlane2 = vecToSphereCenter2.dot(planeNormal);
    distance = distanceToPlane1 > distanceToPlane2 ? distanceToPlane2 : distanceToPlane1;
}

double SphereCylinderBond::GetClampedRoot(double const &slope, double const &h0, double const &h1)
{
    // Theoretically, r is in (0,1). However, when the slope is
    // nearly zero, then so are h0 and h1. Significant numerical
    // rounding problems can occur when using floating-point
    // arithmetic. If the rounding causes r to be outside the
    // interval, clamp it. It is possible that r is in (0,1) and has
    // rounding errors, but because h0 and h1 are both nearly zero,
    // the quadratic is nearly constant on (0,1). Any choice of p
    // should not cause undesirable accuracy problems for the final
    // distance computation.
    //
    // NOTE: You can use bisection to recompute the root or even use
    // bisection to compute the root and skip the division. This is
    // generally slower, which might be a problem for high-performance
    // applications.

    double const zero = static_cast<double>(0);
    double r;
    if (h0 < zero)
    {
        double const one = static_cast<double>(1);
        if (h1 > zero)
        {
            r = -h0 / slope;
            if (r > one)
            {
                r = static_cast<double>(0.5);
            }
            // The slope is positive and -h0 is positive, so there is
            // no need to test for a negative value and clamp it.
        }
        else
        {
            r = one;
        }
    }
    else
    {
        r = zero;
    }
    return r;
}

void SphereCylinderBond::ComputeIntersection(std::array<double, 2> const &sValue,
                                             std::array<int32_t, 2> const &classify, double const &b, double const &f00,
                                             double const &f10, std::array<int32_t, 2> &edge,
                                             std::array<std::array<double, 2>, 2> &end)
{
    // The divisions are theoretically numbers in [0,1]. Numerical
    // rounding errors might cause the result to be outside the
    // interval. When this happens, it must be that both numerator
    // and denominator are nearly zero. The denominator is nearly
    // zero when the segments are nearly perpendicular. The
    // numerator is nearly zero when the P-segment is nearly
    // degenerate (f00 = a is small). The choice of 0.5 should not
    // cause significant accuracy problems.
    //
    // NOTE: You can use bisection to recompute the root or even use
    // bisection to compute the root and skip the division. This is
    // generally slower, which might be a problem for high-performance
    // applications.

    double const zero = static_cast<double>(0);
    double const half = static_cast<double>(0.5);
    double const one = static_cast<double>(1);
    if (classify[0] < 0)
    {
        edge[0] = 0;
        end[0][0] = zero;
        end[0][1] = f00 / b;
        if (end[0][1] < zero || end[0][1] > one)
        {
            end[0][1] = half;
        }

        if (classify[1] == 0)
        {
            edge[1] = 3;
            end[1][0] = sValue[1];
            end[1][1] = one;
        }
        else // classify[1] > 0
        {
            edge[1] = 1;
            end[1][0] = one;
            end[1][1] = f10 / b;
            if (end[1][1] < zero || end[1][1] > one)
            {
                end[1][1] = half;
            }
        }
    }
    else if (classify[0] == 0)
    {
        edge[0] = 2;
        end[0][0] = sValue[0];
        end[0][1] = zero;

        if (classify[1] < 0)
        {
            edge[1] = 0;
            end[1][0] = zero;
            end[1][1] = f00 / b;
            if (end[1][1] < zero || end[1][1] > one)
            {
                end[1][1] = half;
            }
        }
        else if (classify[1] == 0)
        {
            edge[1] = 3;
            end[1][0] = sValue[1];
            end[1][1] = one;
        }
        else
        {
            edge[1] = 1;
            end[1][0] = one;
            end[1][1] = f10 / b;
            if (end[1][1] < zero || end[1][1] > one)
            {
                end[1][1] = half;
            }
        }
    }
    else // classify[0] > 0
    {
        edge[0] = 1;
        end[0][0] = one;
        end[0][1] = f10 / b;
        if (end[0][1] < zero || end[0][1] > one)
        {
            end[0][1] = half;
        }

        if (classify[1] == 0)
        {
            edge[1] = 3;
            end[1][0] = sValue[1];
            end[1][1] = one;
        }
        else
        {
            edge[1] = 0;
            end[1][0] = zero;
            end[1][1] = f00 / b;
            if (end[1][1] < zero || end[1][1] > one)
            {
                end[1][1] = half;
            }
        }
    }
}
void SphereCylinderBond::ComputeMinimumParameters(std::array<int32_t, 2> const &edge,
                                                  std::array<std::array<double, 2>, 2> const &end, double const &b, double const &c,
                                                  double const &e, double const &g00, double const &g10, double const &g01, double const &g11,
                                                  double &s, double &t)
{
    double const zero = static_cast<double>(0);
    double const one = static_cast<double>(1);
    double delta = end[1][1] - end[0][1];
    double h0 = delta * (-b * end[0][0] + c * end[0][1] - e);
    if (h0 >= zero)
    {
        if (edge[0] == 0)
        {
            s = zero;
            t = GetClampedRoot(c, g00, g01);
        }
        else if (edge[0] == 1)
        {
            s = one;
            t = GetClampedRoot(c, g10, g11);
        }
        else
        {
            s = end[0][0];
            t = end[0][1];
        }
    }
    else
    {
        double h1 = delta * (-b * end[1][0] + c * end[1][1] - e);
        if (h1 <= zero)
        {
            if (edge[1] == 0)
            {
                s = zero;
                t = GetClampedRoot(c, g00, g01);
            }
            else if (edge[1] == 1)
            {
                s = one;
                t = GetClampedRoot(c, g10, g11);
            }
            else
            {
                s = end[1][0];
                t = end[1][1];
            }
        }
        else // h0 < 0 and h1 > 0
        {
            double z = std::min(std::max(h0 / (h0 - h1), zero), one);
            double omz = one - z;
            s = omz * end[0][0] + z * end[1][0];
            t = omz * end[0][1] + z * end[1][1];
        }
    }
}
double SphereCylinderBond::getPotentialenergy(const std::unique_ptr<SphereParticle> &sphere1, const std::unique_ptr<SphereParticle> &sphere2)
{
    double potentialenergy = 0;
    const auto fiberproperties = manager->getFiberProperties(type);
    const Eigen::Vector3d &position1 = sphere1->getPosition();
    const Eigen::Vector3d &position2 = sphere2->getPosition();

    Eigen::Vector3d normal_length = position2 - position1;
    double length = normal_length.norm();
    double normal_displacement = (length - fiberproperties->getElementlength());

    potentialenergy += 0.5 * fiberproperties->getNormalstiffnesses() * normal_displacement * normal_displacement;

    potentialenergy += tangentialforce.dot(tangentialforce) / (2* fiberproperties->getShearstiffnesses());

    potentialenergy += bendtorque.dot(bendtorque)  / (2*fiberproperties->getBendingstiffnesses());

    potentialenergy += twisttorque.dot(twisttorque)  / (2*fiberproperties->getTwiststiffnesses());

    return potentialenergy;
}

std::string SphereCylinderBond::save_tostring() const
{
    std::ostringstream ss;
    ss.precision(std::numeric_limits<double>::digits10 + 1);

    ss << "Fiber, "
       << "Bond, " << id << ", " << type.getCategory() << ", " << type.getSubType() << ", " << fiberid
       << ", " << node1 << ", " << node2 << ", " << neighborelement1 << ", " << neighborelement2 << ", " << energyDissipation << ", "
       << tangentialforce.x() << ", " << tangentialforce.y() << ", " << tangentialforce.z() << ", "
       << twisttorque.x() << ", " << twisttorque.y() << ", " << twisttorque.z() << ", "
       << bendtorque.x() << ", " << bendtorque.y() << ", " << bendtorque.z();
    return ss.str();
}
