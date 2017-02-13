//
// Created by vadim on 21.10.16.
//

#ifndef PIC_DICHOTOMOUSSOLVER_H
#define PIC_DICHOTOMOUSSOLVER_H

class BaseDichotomousSolver{
public:
    double maxErrorX;
    double maxErrorY;
    BaseDichotomousSolver();
    virtual double function(double x);
    double solve(double minX, double maxX);
};

class TemperatureRelativisticMaxwellSolver : public BaseDichotomousSolver{
public:
    double alphaNormal;
    double rightPart;

    TemperatureRelativisticMaxwellSolver(double alphaNormalValue, double rightPartValue);
    double function(double x);
};

#endif //PIC_DICHOTOMOUSSOLVER_H
