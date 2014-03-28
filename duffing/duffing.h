#ifndef DUFFING_H_
#define DUFFING_H_

#include <cmath>

class Duffing
{
public:
    Duffing(double dt);
    ~Duffing();

    double time();
    double pos();
    double vel();

    double force();

    void step();

private:
    double F0;
    double GAMMA;
    double OMEGA;

    double t;
    double dt;

    double x;
    double v;
    double xmh;
    double xph;

};

#endif
