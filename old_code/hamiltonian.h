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
        qdot0[0][0] = pax / ma;
        qdot0[0][1] = pay / ma;
        qdot0[1][0] = pbx / mb;
        qdot0[1][1] = pby / mb;
        pdot0[0][0] = top * difX / r;
        pdot0[0][1] = top * difY / r;
        pdot0[1][0] = - top * difX / r;
        pdot0[1][1] = - top * difY / r;
    }
};

#endif /* HAMILTONIAN_H */
