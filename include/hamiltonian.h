#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include "equations.h"

class Hamiltonian : public Equations
{
private:
    double G;
    
    double difX = qbx - qax;
    double difY = qby - qay;
    double r = pow(pow(difX, 2.0) + pow(difY, 2.0), 1.5);
    double top = G * ma * mb;
public:
    Hamiltonian(const Body& body1, const Body& body2, const double& userG):
    G(userG),
    Equations(body1, body2)
    {
        dqax = pax / ma;
        dqay = pay / ma;
        dqbx = pbx / mb;
        dqby = pby / mb;
        dpax = top * difX / r;
        dpay = top * difY / r;
        dpbx = - top * difX / r;
        dpby = - top * difY / r;
    }
};

#endif /* HAMILTONIAN_H */
