//
//  body.cpp
//  Two_Body_Problem
//
//  Created by Zackary Windham on 5/11/20.
//  Copyright © 2020 Zackary Windham. All rights reserved.
//

#include "body.h"

Body::Body(): mass(0)
{
    position[0] = 0;
    position[1] = 0;
    momentum[0] = 0;
    momentum[1] = 0;
}

Body::Body(const double& x, const double& y, const double& px, const double& py, const double& m):  mass(m)
{
    position[0] = x;
    position[1] = y;
    momentum[0] = px;
    momentum[1] = py;
}

double Body::getMass() const
{
    return mass;
}

double Body::getX() const
{
    return position[0];
}

double Body::getY() const
{
    return position[1];
}

double Body::getMomentumX() const
{
    return momentum[0];
}

double Body::getMomentumY() const
{
    return momentum[1];
}

double Body::getP() const
{
    return pow(pow(momentum[0], 2.0) + pow(momentum[1], 2.0), 0.5);
}

double Body::getKE() const
{
    return pow(getP(), 2.0) / (2 * mass);
}

void Body::setPosition(double x, double y)
{
    position[0] = x;
    position[1] = y;
}

void Body::setMomentum(double px, double py)
{
    momentum[0] = px;
    momentum[1] = py;
}

std::string Body::toString() const
{
    std::stringstream os;
    os << position[0] << ", " << position[1];
    
    return os.str();
}
