#include "GA.h"
#include <sstream>
#include <iostream>

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
        Chromosome* dyingChromosome = new Chromosome(losers[i]);

        for (int u = 0; u < populationSize; u++)
        {
            if (P[u]->toString() == dyingChromosome->toString())
            {
                P[u] = offspring[i];
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

Chromosome** GA::selection(FitnessFunction* fitnessFunction)
{
    int j = 0;
    double k1 = 0;
    std::string k2 = "";

    double fitness[populationSize];
    std::string id[populationSize];
    
    Chromosome** p = new Chromosome*[populationSize];

    for (int i = 0; i < populationSize; i++)
    {
        fitness[i] = population[i]->fitness(fitnessFunction, population[i], population[i]->getNumGenes());

        id[i] = population[i]->toString();
    }

    for(int i = 1; i < populationSize; i++)
    {
        k1 = fitness[i];
        k2 = id[i];

        j = i-1;

        while(j >= 0 && fitness[j] < k1)
        {
            fitness[j+1] = fitness[j];
            id[j+1] = id[j];

            j--;
        }

        fitness[j+1] = k1;
        id[j+1] = k2;
    }

    for (int i = 0; i < populationSize; i++)
    {
        for (int j = 0; j < populationSize; j++)
        {
            if (fitness[i] == population[j]->fitness(fitnessFunction, population[j], population[j]->getNumGenes()) && id[i] == population[j]->toString())
            {
                p[i] = new Chromosome (population[j]);
                break;
            }
        }  
    }

    return p;
}

Chromosome** GA::inverseSelection(FitnessFunction* fitnessFunction)
{
    int j = 0;

    Chromosome** p = new Chromosome*[populationSize];

    for (int i = populationSize-1; i >= 0; i--)
    {
        p[j] = new Chromosome(selection(fitnessFunction)[i]);
        j++;
    }
    
    return p;
}

Chromosome** GA::crossOver(Chromosome* c1, Chromosome* c2)
{
    Chromosome** c = new Chromosome*[2];

    c[0] = new Chromosome(c1->crossOver(c2));
    c[1] = new Chromosome(c2->crossOver(c1));

    return c;
}

Chromosome* GA::mutate(Chromosome* c1)
{
    Chromosome* c = new Chromosome(c1->mutate());
    
    return c;
}

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

void GA::setPopulation(Chromosome** p)
{
    for (int i = 0; i < populationSize; i++)
    {
        population[i] = new Chromosome(p[i]);
    }
}

// ;)