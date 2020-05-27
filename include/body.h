
//  body.hpp
//  Two_Body_Problem
//
//  Created by Zackary Windham on 5/11/20.
//  Copyright © 2020 Zackary Windham. All rights reserved.
//

#ifndef BODY_H
#define BODY_H

#include <string>
#include <sstream>

class Body
{
private:
    double mass;
    double position[2];
    double momentum[2];
    std::string toString() const;
public:
    //Constructor
    Body();
    Body(const double& x, const double& y, const double& px, const double& py, const double& m);
    
    //Getters and setters for different data members
    double getMass() const;
    double getX() const;
    double getY() const;
    double getMomentumX() const;
    double getMomentumY() const;
    void setPosition(double x, double y);
    void setMomentum(double px, double py);
    
    friend std::ostream& operator<< (std::ostream& os, const Body& body)
    {
        os << body.toString();
        return os;
    }
};

#endif /* BODY_H */
