#include <iostream>
#include <iomanip>

#include "duffing.h"

int main(int argc, char const *argv[])
{
	double dt = 2.0 * M_PI / 2.4 / 10000;

	Duffing duff(dt);

	for (int i = 0; i < 50000; ++i)
	{
		for (int j = 0; j < 10000; ++j)
		{
			duff.step();
		}
		std::cout << std::setw(15) << duff.pos()
		          << std::setw(15) << duff.vel()
		          << std::endl;
	}

	return 0;
}
