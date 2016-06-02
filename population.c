/*
 * =====================================================================================
 *
 *       Filename:  population.c
 *
 *    Description:  This file contains the functions that involve population.
 *
 *        Version:  1.0
 *        Created:  09/01/16 18:39:44
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Michael Tierney (MT), tiernemi@tcd.ie
 *   Organization:  
 *
 * =====================================================================================
 */

#include <stdlib.h>
#include <stdio.h>
#include "population.h"
#include "chromosome.h"
#include <stdbool.h>
#include <string.h>
#include <mpi/mpi.h>

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  makePopulation
 *    Arguments:  int popSize - The size of the population.
 *                double mutationRate - The mutation rate of the population.
 *                double crossRate - The crossover rate.
 *                int chromoLength - Length of a chromosome.
 *                int numGames - The number of iterations of the prisoners dilemma to play.
 *  Description:  Makes the population object which is used to carry out the genetic
 *                simulation.
 * =====================================================================================
 */

Population * makePopulation (int popSize, double mutationRate, double crossRate, int chromoLength, int numGames) {
	int rank, nproc ;
	int i,j ;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank) ;
	MPI_Comm_size(MPI_COMM_WORLD, &nproc) ;
	Population * newPop = malloc(sizeof(Population)) ;
	newPop->numGames = numGames ;
	newPop->popSize = popSize ;
	newPop->mutationRate = mutationRate ;
	newPop->crossRate = crossRate ;
	newPop->nproc = nproc ;
	// Number of send/Recvs that need to be carried out. //
	newPop->numSends = getNumSends(newPop) ;
	newPop->bitSizeInt = sizeof(unsigned int)*8 ;
	newPop->chromoArSize = (chromoLength / newPop->bitSizeInt) + ((chromoLength % newPop->bitSizeInt) != 0) ;
	newPop->chromoLength = chromoLength ;
	// Program prisoner dilemma outcomes. //
	newPop->prisonerRes[0] = &newPop->outcomes[0] ;
	newPop->prisonerRes[1] = &newPop->outcomes[4] ;
	newPop->prisonerRes[0][0] = 3 ; // Co-op
	newPop->prisonerRes[1][0] = 3 ; // Co-op
	newPop->prisonerRes[0][1] = 0 ; // Prisoner 2 defects
	newPop->prisonerRes[1][1] = 5 ; // Prisoner 2 defects
	newPop->prisonerRes[0][2] = 5 ; // Prisoner 1 defects
	newPop->prisonerRes[1][2] = 0 ; // Prisoner 1 defects
	newPop->prisonerRes[0][3] = 1 ; // Defect
	newPop->prisonerRes[1][3] = 1 ; // Defect

	// Master is the only thread that requires this storage. //
	if (rank == 0) {
		// Array storing the chromosomes for the current population. //
		newPop->organisms = malloc(sizeof(Chromosome)*popSize) ;
		// Array storing the chromosomes for the next population. //
		newPop->offspring = malloc(sizeof(Chromosome)*popSize) ;
		// CDF for choosing the next population. //
		newPop->fitnessCDF = malloc(sizeof(double)*popSize) ;
		newPop->fitness = malloc(sizeof(int)*popSize) ;
		newPop->genePool = malloc(sizeof(unsigned int)*popSize*newPop->chromoArSize) ;
		newPop->offGenePool = malloc(sizeof(unsigned int)*popSize*newPop->chromoArSize) ;
		// Initialise the chromosomes. //
		for (i = 0 ; i < popSize ; ++i) {
			initialiseChromosome(&newPop->organisms[i], newPop, &newPop->genePool[i*newPop->chromoArSize]) ;
			initialiseChromosome(&newPop->offspring[i], newPop, &newPop->offGenePool[i*newPop->chromoArSize]) ;
		}
		// Create array containing the pair indices of chromosomes that need to be evaluated. //
		newPop->sendRecvIndexPairs = malloc(sizeof(int*)*newPop->numSends) ;
		newPop->indexPairs = malloc(sizeof(int)*2*newPop->numSends) ;
		for (i = 0 ; i < newPop->numSends ; ++i) {
			newPop->sendRecvIndexPairs[i] = &newPop->indexPairs[2*i] ;
		}
		int index = 0 ;
		for (i = 0 ; i < newPop->popSize ; ++i) {
			for (j = i+1 ; j < newPop->popSize ; ++j) {
				newPop->sendRecvIndexPairs[index][0] = i ;
				newPop->sendRecvIndexPairs[index][1] = j ;
				++index ;
			}
		}
	}
	return newPop ;
}	/* -----  end of function makePopulation  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  simulateM
 *    Arguments:  Population * pop - Population to simulate.
 *                int numGenerations - Number of generations to run natural selection
 *                  simulations.
 *                double * fitnessResults - The array storing the population fitness
 *                  for each generation.
 *  Description:  This function performs the natural selection genetic algorithm by
 *                selecting a population based on some fitness function, breeding the
 *                winners and mutating the resultant population. Only used by the master
 *                thread. This thread delegates the fitness evaluation to slaves.
 * =====================================================================================
 */

void simulateM (Population * pop, int numGenerations, double * fitnessResults) {
	int i = 0 ;
	for (i = 0 ; i < numGenerations ; ++i) {
		// Get the slaves to calculate the fitness of each member of pop. //
		calcFitnessPop(pop) ;
		// Save the result for the population fitness. //
		fitnessResults[i] = pop->fitnessCDF[pop->popSize-1] ;
		// Use the roulette method to select the population + cross population. //
		selectPop(pop) ;
		// Mutate the new population. //
		mutatePop(pop) ;
	}
	// Terminate the slaves. //
	int termSignal = 0 ;
	MPI_Request requestTerm ;
	MPI_Request requestDummy ;
	for (i = 1 ; i < pop->nproc ; ++i) {
		// Send the termination signal. //
		MPI_Isend(&termSignal, 1, MPI_INT, i, 1, MPI_COMM_WORLD, &requestTerm) ;
		// Dummy signal to get slave to check state. //
		MPI_Isend(&termSignal, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &requestDummy) ;
		MPI_Request_free(&requestTerm) ;
		MPI_Request_free(&requestDummy) ;
	}
}		/* -----  end of function mutate  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  simulateS
 *    Arguments:  Population * pop - Struct containing info about chromosomes.
 *  Description:  This slave thread waits for chromonsomes to be sent to it and then
 *                evalauates the given pair. Sends back an array containing the fitness.
 * =====================================================================================
 */

void simulateS (Population * pop) {
	// Two dummy chromosomes. //
	Chromosome chromo1 ;
	Chromosome chromo2 ;
	int i ;
	int hist = 0 ;
	int threadOpen = 1 ;
	// Genepool used to receive genetic data. //
	unsigned int * genePool = malloc(sizeof(unsigned int)*2*pop->chromoArSize) ;
	int * fitness = malloc(sizeof(int)*2) ;
	MPI_Status stat ;
	MPI_Request requestTerm ;
	MPI_Request requestSend ;
	// Posts receive for the termination signal. //
	MPI_Irecv(&threadOpen, 1 , MPI_INT, 0, 1, MPI_COMM_WORLD, &requestTerm) ;

	// Waits to receive genes, initialises chromosomes and then evaluates fitness.
	while (1) {
		// Receive genetic info. //
		MPI_Recv(genePool, 2*pop->chromoArSize , MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD, &stat) ;
		if (!threadOpen) {
			break ;
		}
		chromo1.genes = &genePool[0] ;
		chromo2.genes = &(genePool[pop->chromoArSize]) ;
		memset(fitness, 0, sizeof(int)*2) ;
		// Play PD game multiple times. //
		hist = rand() % pop->chromoLength ;
		for (i = 0 ; i < pop->numGames ; ++i) {
			calcFitness(&chromo1, &chromo2, pop, &hist, &fitness[0], &fitness[1]) ;
		}
		// Send the fitness back to the master. //
		MPI_Isend(fitness, 2, MPI_INT, 0, 0, MPI_COMM_WORLD, &requestSend) ;
		// Free request to avoid memory leak. //
		MPI_Request_free(&requestSend) ;
	}
	MPI_Request_free(&requestTerm) ;
	free(genePool) ;
	free(fitness) ;
}		/* -----  end of function simulateS  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  mutatePop
 *    Arguments:  Population * pop - Struct containing the chromosomes to be mutated.
 *  Description:  Mutates each member of the current population.
 * =====================================================================================
 */

void mutatePop (Population * pop) {
	int i = 0 ;
	for (i = 0 ; i < pop->popSize ; ++i) {
		mutate(&pop->organisms[i], pop->mutationRate, pop) ;
	}
}	/* -----  end of function mutate  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  calcFitnessPop
 *    Arguments:  Population * pop - Population for which we are calculating the fitness.
 *  Description:  Calculates the fitness of the population by delegating the evaluation
 *                to slave threads. Evaluates fitnesses in chunks of nproc-1.
 * =====================================================================================
 */

void calcFitnessPop (Population * pop) {
	int i = 0 ;

	memset(pop->fitness, 0, sizeof(int)*pop->popSize) ;
	unsigned int * geneticPackage = malloc(sizeof(unsigned int)*pop->chromoArSize*2) ;
	int * fitnessPackage = malloc(sizeof(unsigned int)*2*pop->nproc) ;
	MPI_Request * requestsR = malloc(sizeof(MPI_Request)*(pop->nproc-1)) ;
	MPI_Request * requestsS = malloc(sizeof(MPI_Request)*(pop->nproc-1)) ;
	int * sourceIndices = malloc(sizeof(int)*(pop->nproc-1)) ;

	// Make the organisms compete with each other in succesive PD games. //
	int dest = 0 ;

	int jobsSent = 0 ;
	int jobsRecv = 0 ;
	// Send initial data. //
	for (i = 0 ; i < pop->nproc-1 ; ++i) {
		// Extract genetic info in preparation for sending. //
		memcpy(&geneticPackage[0], pop->organisms[pop->sendRecvIndexPairs[i][0]].genes, sizeof(unsigned int)*pop->chromoArSize) ;
		memcpy(&geneticPackage[pop->chromoArSize], pop->organisms[pop->sendRecvIndexPairs[i][1]].genes, sizeof(unsigned int)*pop->chromoArSize) ;
		// Determine the destination thread. //
		dest = (i%(pop->nproc-1))+1 ;
		sourceIndices[dest-1] = i ;
		// Send the genetic info. //
		MPI_Isend(geneticPackage, 2*pop->chromoArSize, MPI_UNSIGNED, dest, 0, MPI_COMM_WORLD, &requestsS[dest-1]) ;
		// Free request to avoid memory leak. //
		MPI_Request_free(&requestsS[dest-1]) ;
		// Receive the resulting fitness. //
		MPI_Irecv(&fitnessPackage[2*(dest-1)], 2, MPI_INT, dest, 0, MPI_COMM_WORLD, &requestsR[dest-1]) ;
		++jobsSent ;
	}

	int indexComp = 0 ;
	int jobCompleted = 0 ;
	MPI_Status recvStat ;
	// If jobs are finished send a new job. //
	while (jobsSent < pop->numSends) {
		MPI_Testany(pop->nproc-1, requestsR, &indexComp, &jobCompleted, &recvStat) ;
		if (jobCompleted) {
			pop->fitness[pop->sendRecvIndexPairs[sourceIndices[indexComp]][0]] += fitnessPackage[2*indexComp] ;
			pop->fitness[pop->sendRecvIndexPairs[sourceIndices[indexComp]][1]] += fitnessPackage[(2*indexComp)+1] ;
			++jobsRecv ;

			// Extract genetic info in preparation for sending. //
			memcpy(&geneticPackage[0], pop->organisms[pop->sendRecvIndexPairs[jobsSent][0]].genes, sizeof(unsigned int)*pop->chromoArSize) ;
			memcpy(&geneticPackage[pop->chromoArSize], pop->organisms[pop->sendRecvIndexPairs[jobsSent][1]].genes, sizeof(unsigned int)*pop->chromoArSize) ;
			// Determine the destination thread. //
			dest = indexComp+1 ;
			sourceIndices[indexComp] = jobsSent ;
			// Send the genetic info. //
			MPI_Isend(geneticPackage, 2*pop->chromoArSize, MPI_UNSIGNED, dest, 0, MPI_COMM_WORLD, &requestsS[indexComp]) ;
			// Free request to avoid memory leak. //
			MPI_Request_free(&requestsS[indexComp]) ;
			// Receive the resulting fitness. //
			MPI_Irecv(&fitnessPackage[2*(indexComp)], 2, MPI_INT, dest, 0, MPI_COMM_WORLD, &requestsR[indexComp]) ;
			++jobsSent ;
			jobCompleted = 0 ;
		}
	}

	// Receive final data sends. //
	while (jobsRecv < pop->numSends) {
		MPI_Testany(pop->nproc-1, requestsR, &indexComp, &jobCompleted, &recvStat) ;
		if (jobCompleted) {
			pop->fitness[pop->sendRecvIndexPairs[sourceIndices[indexComp]][0]] += fitnessPackage[2*indexComp] ;
			pop->fitness[pop->sendRecvIndexPairs[sourceIndices[indexComp]][1]] += fitnessPackage[(2*indexComp)+1] ;
			++jobsRecv ;
			jobCompleted = 0 ;
		}
	}

	free(geneticPackage) ;
	free(fitnessPackage) ;
	free(requestsR) ;
	free(requestsS) ;
	free(sourceIndices) ;

	// Generate the CDF. //
	pop->fitnessCDF[0] = pop->fitness[0] ;
	for (i = 1 ; i < pop->popSize ; ++i) {
		pop->fitnessCDF[i] = pop->fitness[i] + pop->fitnessCDF[i-1] ;
	}
}		/* -----  end of function calcFitnessPop  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  selectPop
 *    Arguments:  Population * pop - Population struct containing populations.
 *  Description:  Normalises the fitness cdf, and uses this to generate a random variable.
 *                This RV is used to select more fit organisms with a higher probability.
 *                These "parents" have a chance to cross and produce a new organism.
 *                Finally this new population replaces the old population.
 * =====================================================================================
 */

void selectPop (Population * pop) {
	int i ;
	// Normalise the fitness cdf. //
	for (i = 0 ; i < pop->popSize ; ++i) {
		pop->fitnessCDF[i] /= pop->fitnessCDF[pop->popSize-1] ;
	}

	// Generate new population in groups of two. //
	for (i = 0 ; i < pop->popSize-1 ; i+=2) {
		double randomNum = (double) (rand() / (double) RAND_MAX) ;
		// Parent 1. //
		int indexPar1 = 0 ;
		while (pop->fitnessCDF[indexPar1] <= randomNum) {
			++indexPar1 ;
		}
		
		bool notSame = false ;
		// Parent 2. //
		int indexPar2 = 0 ;
		while (!notSame) {
			randomNum = (double) (rand() / (double) RAND_MAX) ;
			indexPar2 = 0 ;
			while (pop->fitnessCDF[indexPar2] <= randomNum) {
				++indexPar2 ;
			}
			notSame = (indexPar2 != indexPar1) ;
		}

		// Chance of breeding to produce new organisms. //
		randomNum = (double) (rand() / (double) RAND_MAX) ;
		if (randomNum <= pop->crossRate) {
			doubleCrossOver(&pop->organisms[indexPar1], &pop->organisms[indexPar2], &pop->offspring[i], &pop->offspring[i+1], pop) ;	
		} else {
			memcpy(pop->offspring[i].genes, pop->organisms[indexPar1].genes, sizeof(unsigned int)*pop->chromoArSize) ;
			memcpy(pop->offspring[i+1].genes, pop->organisms[indexPar2].genes, sizeof(unsigned int)*pop->chromoArSize) ;
		}
	}
	// If odd number of organisms there will be a hanging index. This fills that index. //
	if (pop->popSize % 2 != 0) {
		double randomNum = (double) (rand() / (double) RAND_MAX) ;
		// Parent 1. //
		int indexPar1 = 0 ;
		while (pop->fitnessCDF[indexPar1] <= randomNum) {
			++indexPar1 ;
		}
		
		bool notSame = false ;
		// Parent 2. //
		int indexPar2 = 0 ;
		while (!notSame) {
			randomNum = (double) (rand() / (double) RAND_MAX) ;
			indexPar2 = 0 ;
			while (pop->fitnessCDF[indexPar2] <= randomNum) {
				++indexPar2 ;
			}
			notSame = (indexPar2 != indexPar1) ;
		}

		// Chance of breeding to produce new organism. //
		randomNum = (double) (rand() / (double) RAND_MAX) ;
		if (randomNum <= pop->crossRate) {
			singleCrossOver(&pop->organisms[indexPar1], &pop->organisms[indexPar2], &pop->offspring[pop->popSize-1], pop) ;
		} else {
			memcpy(pop->offspring[pop->popSize-1].genes, pop->organisms[indexPar1].genes, sizeof(unsigned int)*pop->chromoArSize) ;
		}
	}
	
	// Swap offspring and parents. //
  	Chromosome * temp = pop->organisms ;
	pop->organisms = pop->offspring ;
	pop->offspring = temp ;
}		/* -----  end of function selectPop  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  printPopulation
 *    Arguments:  Population * pop - Population to print. 
 *  Description:  Prints all the chromosomes of the population.
 * =====================================================================================
 */

void printPopulation (Population * pop) {
	int i ;
	for (i = 0 ; i < pop->popSize ; ++i) {
		printChromosome(&pop->organisms[i], pop) ;
	}
}		/* -----  end of function printPopulation  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  freePopulation
 *    Arguments:  Population * pop - Population we want to free.
 *  Description:  Frees the dynamic memory associated with pop.
 * =====================================================================================
 */

void freePopulation (Population * pop) {
	free(pop->organisms) ;
	free(pop->offspring) ;
	free(pop->genePool) ;
	free(pop->offGenePool) ;
	free(pop->fitnessCDF) ;
	free(pop->fitness) ;
	free(pop->sendRecvIndexPairs) ;
	free(pop->indexPairs) ;
	free(pop) ;
}		/* -----  end of function freePopulation  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  getNumSends
 *    Arguments:  
 *  Description:  
 * =====================================================================================
 */

unsigned int getNumSends (Population * pop) {
	unsigned int ans = ((pop->popSize-1)*pop->popSize/2) ;
	return ans ;
}		/* -----  end of function getNumSends  ----- */
