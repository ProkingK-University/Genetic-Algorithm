#include "FitnessFunction.h"

double FitnessFunction::calculateFitness(Chromosome* chromosome, int numGenes)
{
    double m = 0;
    double fitness = 0;

    for (int i = 0; i < numGenes; i++)
    {
        if (chromosome->getGenes()[i] == true)
        {
            m++;
        }
    }

    fitness = m / numGenes;

    return fitness;
}

// ;)