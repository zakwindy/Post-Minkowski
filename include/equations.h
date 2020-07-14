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
    
    double dqax;
    double dqay;
    double dqbx;
    double dqby;
    double dpax;
    double dpay;
    double dpbx;
    double dpby;
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
        return dqax;
    }
    
    double getdqay()
    {
        return dqay;
    }
    
    double getdqbx()
    {
        return dqbx;
    }
    
    double getdqby()
    {
        return dqby;
    }
    
    double getdpax()
    {
        return dpax;
    }
    
    double getdpay()
    {
        return dpay;
    }
    
    double getdpbx()
    {
        return dpbx;
    }
    
    double getdpby()
    {
        return dpby;
    }
};

#endif /* EQUATIONS_H */
