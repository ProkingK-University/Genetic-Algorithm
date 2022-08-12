#include "Chromosome.h"

Chromosome::Chromosome(int numGenes, RandomGenerator* rand)
{
    if (numGenes < 0)
    {
        this->numGenes = 0;
    }
    else
    {
        this->numGenes = numGenes;
    }

    genes = new bool[numGenes];

    if (rand == NULL)
    {
        for (int i = 0; i < numGenes; i++)
        {
            genes[i] = false;
        }
    }
    else
    {
        for (int i = 0; i < numGenes; i++)
        {
            genes[i] = rand->randomBool();
        }
    }
}

Chromosome::Chromosome(Chromosome* chromosome)
{
    if (chromosome == NULL)
    {
        numGenes = 0;

        genes = new bool[0];
    }
    else
    {
        numGenes = chromosome->getNumGenes();

        genes = new bool[numGenes];

        for (int i = 0; i < numGenes; i++)
        {
            genes[i] = chromosome->getGenes()[i];
        }
    }  
}

Chromosome::Chromosome(bool* genes, int numGenes)
{
    if (numGenes < 0)
    {
        this->numGenes = 0;

        genes = new bool[0];
    }
    else
    {
        this->numGenes = numGenes;

        this->genes = new bool[this->numGenes];

        if (genes == NULL)
        {
            for (int i = 0; i < this->numGenes; i++)
            {
                genes[i] = false;
            }
        }
        else
        {
            for (int i = 0; i < this->numGenes; i++)
            {
                this->genes[i] = genes[i];
            }
        }
    }
}

Chromosome::~Chromosome()
{
    genes = NULL;
    delete [] genes;
}

double Chromosome::fitness(FitnessFunction* fitnessFunction, Chromosome* chromosome, int numGenes)
{
    if (fitnessFunction == NULL || chromosome == NULL || numGenes <= 0)
    {
        return 0;
    }
    else
    {
        return fitnessFunction->calculateFitness(chromosome, numGenes);
    }
}

int Chromosome::getNumGenes() {return numGenes;}

Chromosome* Chromosome::crossOver(Chromosome* c2)
{
    Chromosome* c;

    if (c2 == NULL)
    {
        c = this;

        return c;
    }
    else
    {
        int crossOverpoint = floor(c2->getNumGenes() / 2);

        bool* nGenes = new bool[numGenes];

        for (int i = 0; i < crossOverpoint; i++)
        {
            nGenes[i] = this->genes[i];
        }

        for (int i = crossOverpoint; i < numGenes; i++)
        {
            nGenes[i] = c2->getGenes()[i];
        }

        c = new Chromosome(nGenes, numGenes);

        nGenes = NULL;
        delete [] nGenes;

        return c;
    }
    
}

Chromosome* Chromosome::mutate()
{
    bool* nGenes = new bool[numGenes];

    for (int i = numGenes-1; i >= 0; i--)
    {
        if (genes[i] == true)
        {
            nGenes[i] = false;
        }
        else
        {
            nGenes[i] = true;
        }
    }

    Chromosome* c = new Chromosome(nGenes, numGenes);

    nGenes = NULL;
    delete [] nGenes;

    return c;
}

std::string Chromosome::toString()
{
    if (numGenes == 0)
    {
        return "";
    }
    else
    {
        std::string str = "";

        for (int i = 0; i < numGenes; i++)
        {
            if (genes[i] == true)
            {
                str += "1";
            }
            else if(genes[i] == false)
            {
                str += "0";
            }
        }

        return str;
    }
}

bool* Chromosome::getGenes() {return genes;}

// ;)