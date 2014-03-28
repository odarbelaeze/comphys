#include "duffing.h"

Duffing::Duffing(double dt)
{
	F0 = 2.0;
	GAMMA = 0.1;
	OMEGA = 2.4;

	t = 0.0;
	this -> dt = dt;

	x = 0.5;
	v = 0.0;

	xph = x + v * dt + 0.5 * dt * dt * force();
}


Duffing::~Duffing()
{

}


double Duffing::time()
{
    return t;
}


double Duffing::pos()
{
	return x;
}


double Duffing::vel()
{
	return v;
}


double Duffing::force()
{
	return - GAMMA * v + x - x * x * x + F0 * std::cos(OMEGA * t);
}


void Duffing::step()
{
	t = t + dt;
	xmh = x;
	x = xph;
	xph = 2.0 * x - xmh + dt * dt * force();
	v = (xph - xmh) / (2.0 * dt);
}
