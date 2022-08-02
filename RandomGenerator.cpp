#include "RandomGenerator.h"

RandomGenerator::RandomGenerator(int seed) : seed(seed)
{
    srand(seed);
}

bool RandomGenerator::randomBool()
{
    if (rand()%2 == 0)
    {
        return false;
    }
    else
    {
        return true;
    }
}

// ;)