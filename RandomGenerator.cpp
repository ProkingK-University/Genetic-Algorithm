#include "RandomGenerator.h"

RandomGenerator::RandomGenerator(int seed) : seed(seed)
{
    srand(seed);
}

//Generates a true or false value depending on if the random number is odd or even
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