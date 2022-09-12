#include "FitnessFunction.h"

//Calculates the fittness of the chromosome by counting the number of true and dividing by the number of genes
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