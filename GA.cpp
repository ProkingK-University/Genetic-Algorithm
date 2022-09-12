#include "GA.h"
#include <iostream>

//Creates an array of chromosomes
GA::GA(int populationSize, RandomGenerator* rand, int numGenes, int selectionSize) : populationSize(populationSize), selectionSize(selectionSize)
{
    population = new Chromosome*[populationSize];

    for (int i = 0; i < populationSize; i++)
    {
        population[i] = new Chromosome(numGenes, rand);
    }
}

GA::GA(GA* geneticAlgorithm)
{
    populationSize = geneticAlgorithm->populationSize;
    selectionSize = geneticAlgorithm->selectionSize;

    population = new Chromosome*[populationSize];

    for (int i = 0; i < populationSize; i++)
    {
        population[i] = new Chromosome(geneticAlgorithm->population[i]);
    }
}

GA::~GA()
{
    for (int i = 0; i < populationSize; i++)
    {
        population[i] = NULL;
        delete population[i];
    }

    population = NULL;
    delete [] population;
}

// Runs the genetic algorithm
Chromosome** GA::run(FitnessFunction* fitnessFunction)
{
    Chromosome** winners = selection(fitnessFunction);
    Chromosome** losers = inverseSelection(fitnessFunction);

    Chromosome** offspring = new Chromosome*[3*selectionSize];
    Chromosome** P = new Chromosome*[populationSize];

    for (int i = 0; i < 2*selectionSize; i++)
    {
        Chromosome** nChromosomes = crossOver(winners[i], winners[i+1]);

        offspring[i] = nChromosomes[0];
        offspring[i+1] = nChromosomes[1];
        i++;

        for (int j = 0; j < 2; j++)
        {
            nChromosomes[j] = NULL;
            delete nChromosomes[j];
        }

        nChromosomes = NULL;
        delete [] nChromosomes;
    }

    for (int i = 0; i < selectionSize; i++)
    {
        offspring[i+2*selectionSize] = new Chromosome(mutate(winners[i+2*selectionSize]));
    }
    
    P = population;

    for (int i = 0; i < 3*selectionSize; i++)
    {
        Chromosome* dyingChromosome = losers[i];

        for (int u = 0; u < populationSize; u++)
        {
            if (P[u] == dyingChromosome)
            {
                P[u] = new Chromosome(offspring[i]);
                break;
            }
        }

        delete dyingChromosome;
    }
    
    for (int i = 0; i < populationSize; i++)
    {
        winners[i] = NULL;
        losers[i] = NULL;

        delete winners[i];
        delete losers[i];
    }

    winners = NULL;
    losers = NULL;

    delete [] winners;
    delete [] losers;

    
    for (int i = 0; i < 3*selectionSize; i++)
    {
        offspring[i] = NULL;

        delete offspring[i];
    }

    offspring = NULL;
    delete [] offspring;

    return P;
}

//Runs multiple generations
double** GA::run(FitnessFunction* fitnessFunction, int numGenerations)
{
    double** results = new double*[numGenerations];

    for (int i = 0; i < numGenerations; i++)
    {
        results[i] = new double[3];
    }

    for (int i = 0; i < numGenerations; i++)
    {
        Chromosome** P = run(fitnessFunction);

        population = P;

        results[i][0] = calculateAvgAccuracy(fitnessFunction);
        results[i][1] = calculateStd(fitnessFunction);
        results[i][2] = calculateVariance();
    }
    
    return results;
}

//Sorts the population from strongest to weakest
Chromosome** GA::selection(FitnessFunction* fitnessFunction)
{
    int j = 0;
    double k1 = 0;
    Chromosome* k2;;

    double fitness[populationSize];
    
    Chromosome** p = new Chromosome*[populationSize];

    for (int i = 0; i < populationSize; i++)
    {
        p[i] = population[i];
    }

    for (int i = 0; i < populationSize; i++)
    {
        fitness[i] = population[i]->fitness(fitnessFunction, population[i], population[i]->getNumGenes());
    }

    for(int i = 1; i < populationSize; i++)
    {
        k1 = fitness[i];
        k2 = p[i];

        j = i-1;

        while(j >= 0 && fitness[j] < k1)
        {
            fitness[j+1] = fitness[j];
            p[j+1] = p[j];

            j--;
        }

        fitness[j+1] = k1;
        p[j+1] = k2;
    }

    return p;
}

//Sorts the population from weakest to strongest
Chromosome** GA::inverseSelection(FitnessFunction* fitnessFunction)
{
    int j = 0;

    Chromosome** p = new Chromosome*[populationSize];

    for (int i = populationSize-1; i >= 0; i--)
    {
        p[j] = selection(fitnessFunction)[i];
        j++;
    }
    
    return p;
}

//Mixes two chromosomes
Chromosome** GA::crossOver(Chromosome* c1, Chromosome* c2)
{
    Chromosome** c = new Chromosome*[2];

    c[0] = new Chromosome(c1->crossOver(c2));
    c[1] = new Chromosome(c2->crossOver(c1));

    return c;
}

//Inverts the passed in chromosomes
Chromosome* GA::mutate(Chromosome* c1)
{
    Chromosome* c = new Chromosome(c1->mutate());
    
    return c;
}

//Calculates avarage fitness
double GA::calculateAvgAccuracy(FitnessFunction* fitnessFunction)
{
    double fitness = 0;

    for (int i = 0; i < populationSize; i++)
    {
        fitness += population[i]->fitness(fitnessFunction, population[i], population[i]->getNumGenes());
    }

    fitness /= populationSize;

    return fitness;
}

//Calculates the standard deviation
double GA::calculateStd(FitnessFunction* fitnessFunction)
{
    int i = 0;
    double sum = 0;
    double mean = 0;
    double stD = 0;

  for(i = 0; i < populationSize; i++)
  {
    sum += population[i]->fitness(fitnessFunction, population[i], population[i]->getNumGenes());
  }

  mean = sum / populationSize;

  for(i = 0; i < populationSize; ++i) 
  {
    stD += pow(population[i]->fitness(fitnessFunction, population[i], population[i]->getNumGenes()) - mean, 2);
  }

  return sqrt(stD / populationSize);
}

//Calculates variance
double GA::calculateVariance()
{
    double var = 0;
    double numChromosomes = 0;

    for (int i = 0; i < populationSize; i++)
    {
        int j = 0;
        for (j = 0; j < i; j++)
        {
            if (population[i]->toString() == population[j]->toString())
            {
                break;
            }
        }

        if (i == j)
        {
            numChromosomes++;
        }
    }

    var = numChromosomes / populationSize;

    return var;
}

//Setter for population array
void GA::setPopulation(Chromosome** p)
{
    for (int i = 0; i < populationSize; i++)
    {
        population[i] = new Chromosome(p[i]);
    }
}

// ;)