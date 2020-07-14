#ifndef PM_H
#define PM_H

#include "equations.h"

class PM : public Equations
{
private:
    double G;
public:
    PM(const Body& body1, const Body& body2, const double& userG);
};

#endif /* PM_H */
