#ifndef EQUATIONS_H
#define EQUATIONS_H

#include "body.h"

class Equations
{
protected:
    double qax;
    double qay;
    double qbx;
    double qby;
    double pax;
    double pay;
    double pbx;
    double pby;
    double ma;
    double mb;
    
    //An array which contains first the positions of both bodies, then the momentums
    double qdot0[2][2];
    double pdot0[2][2];
public:
    Equations(const Body& body1, const Body& body2):
    qax(body1.getX()),
    qay(body1.getY()),
    qbx(body2.getX()),
    qby(body2.getY()),
    pax(body1.getMomentumX()),
    pay(body1.getMomentumY()),
    pbx(body2.getMomentumX()),
    pby(body2.getMomentumY()),
    ma(body1.getMass()),
    mb(body2.getMass())
    {}
    
    double getdqax()
    {
        return qdot0[0][0];
    }
    
    double getdqay()
    {
        return qdot0[0][1];
    }
    
    double getdqbx()
    {
        return qdot0[1][0];
    }
    
    double getdqby()
    {
        return qdot0[1][1];
    }
    
    double getdpax()
    {
        return pdot0[0][0];
    }
    
    double getdpay()
    {
        return pdot0[0][1];
    }
    
    double getdpbx()
    {
        return pdot0[1][0];
    }
    
    double getdpby()
    {
        return pdot0[1][1];
    }
};

#endif /* EQUATIONS_H */
