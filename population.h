#ifndef POPULATION_H_

#define POPULATION_H_

/*
 * =====================================================================================
 *
 *       Filename:  population.h
 *
 *    Description:  Header file containing the struct Population which encapsulates
 *                  the information required to carry out a genetic algorithm.
 *
 *        Version:  1.0
 *        Created:  02/01/16 14:31:36
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Michael Tierney (MT), tiernemi@tcd.ie
 *   Organization:  
 *
 * =====================================================================================
 */

/* 
 * ===  STRUCT  ======================================================================
 *         Name:  Population
 *       Fields:  struct Chromome * organisms - Array of all the organisms in the
 *                population.
 *                struct Chromome * offspring - Array used to store offspring.
 *                unsigned int * genePool - Array containing all the genes in the system.
 *                unsigned int * offGenePool - Array containing all the genes of offspring.
 *                double * fitness - Array containing the fitnesses of the population. 
 *                double * fitnessCDF - Cumalitive fitness of population where the 
 *                last element is the total fitness.
 *                int popSize - The size of the population i.e number of chromosomes.
 *                double mutationRate - The rate of mutation as a % occurence.
 *                double crossRate - The rate of crossover as a % occurence.
 *                int chromoArSize - The number of ints required to reprensent the genes.
 *                int chromolength - Length of the chromosome.
 *                int numGames - Number of games of prisoners dilemma to play.
 *                int bitsizeint - The number of bits in an int.
 *                int nproc - Number of processes.
 *                int numSends - The number of sends/recv pairs that need to be sent.
 *                int ** sendRecvIndexPairs - Indices of the pairs of chromosomes that
 *                that need to be tested against each other.
 *                int * indexPairs - Array storing the actual indices.
 *                int * prisonerRes[2] - Array storing player outcomes.
 *                int outcomes[8] - Array storing all outcomes of prisoners dilemma.
 *  Description:  Population struct containing all the information required to
 *                perform a genetic algorithm.
 * =====================================================================================
 */

struct Chromosome ;

typedef struct Population {
	struct Chromosome * organisms ;
	struct Chromosome * offspring ;
	int * fitness ;
	double * fitnessCDF ;
	unsigned int * genePool ; 
	unsigned int * offGenePool ;
	int popSize ;
	double mutationRate ;
	double crossRate ;
	int chromoArSize ;
	int chromoLength ;
	int numGames ;
	int bitSizeInt ;
	int nproc ;
	int numSends ;
	int ** sendRecvIndexPairs ;
	int * indexPairs ;
	int * prisonerRes[2] ;
	int outcomes[8] ;
} Population ;		/* ----------  end of struct Population  ---------- */

Population * makePopulation(int popSize, double mutationRate, double crossRate, int chromoLength, int numGames) ;
void simulateM(Population * pop, int numGenerations, double * fitnessResults) ;
void simulateS(Population * pop) ;
void calcFitnessPop(Population * pop) ;
void mutatePop(Population * pop) ;
void selectPop(Population * pop) ;
void printPopulation(Population * pop) ;
void freePopulation(Population * pop) ;
unsigned int getNumSends(Population * pop) ;

#endif
