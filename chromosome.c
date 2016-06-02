/*
 * =====================================================================================
 *
 *       Filename:  chromosome.c
 *
 *    Description:  File containing the functions which involve Chromosome.
 *
 *        Version:  1.0
 *        Created:  02/01/16 15:51:13
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
#include <math.h>
#include <string.h>
#include "chromosome.h"
#include "population.h"

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  initialiseChromosome
 *    Arguments:  Population * pop - Population struct containing all info for chromosome
 *                dimensions.
 *  Description:  Makes a chromsome object. Genes are stored as bits within the array
 *                genes.
 * =====================================================================================
 */

void initialiseChromosome(Chromosome * chromo, Population * pop, unsigned int * arrayStart) {
	chromo->genes = arrayStart ;
	mutate(chromo, 0.5, pop) ;
}		/* -----  end of function makeChromosome  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  mutate
 *    Arguments:  Chromosome * chromo - Chromosome to mutate.
 *                double rate - The mututation rate per gene.
 *                Population * pop - Population struct containing all info for chromosome
 *                dimensions.
 *  Description:  Mutates the chromosome at some rate by generating a mutation bitmask
 *                and logically ^ it with the chromosome.
 * =====================================================================================
 */

void mutate(Chromosome * chromo, double rate, Population * pop) {
	unsigned int * mutationBits = calloc(pop->chromoArSize, sizeof(unsigned int)) ;
	int overflow = (pop->chromoLength % (pop->bitSizeInt)) + pop->bitSizeInt*((pop->chromoLength % (pop->bitSizeInt)) == 0) ; 
	int i,j ;
	for (i = 0 ; i < pop->chromoArSize-1 ; ++i) {
		for (j = 0 ; j < pop->bitSizeInt ; ++j) {
			mutationBits[i] <<= 1 ;
			mutationBits[i] += ((double) rand()/ (double) RAND_MAX) < rate ;
		}
	}
	// Truncates the final part of the array to the right size. //
	for (j = 0 ; j < overflow ; ++j) {
		mutationBits[pop->chromoArSize-1] <<= 1 ;	
		mutationBits[pop->chromoArSize-1] += ((double) rand()/ (double) RAND_MAX) < rate ;
	}
	for (i = 0 ; i < pop->chromoArSize ; ++i) {
		chromo->genes[i] ^= mutationBits[i] ;
	}
	free(mutationBits) ;
}		/* -----  end of function mutate  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  doubleCrossOver
 *    Arguments:  Chromosome * par1 - The chromosome of the first parent.
 *                Chromosome * par2 - The chromosome of the second parent.
 *                Chromosome * off1 - The chromosome of the first offspring.
 *                Chromosome * off2 - The chromosome of the second offspring.
 *                Population * pop - Population struct containing all info for chromosome
 *                dimensions.
 *  Description:  Crosses over two chromosomes at a random bit and overwrites to form
 *                new offspring.
 * =====================================================================================
 */

void doubleCrossOver(Chromosome * par1, Chromosome * par2, Chromosome * off1, Chromosome * off2, Population * pop) {
	// Generate to the location of the crossover. //
	int crossOverPoint = rand() % pop->chromoLength ;
	// Index of the array holding that group of 32 bits. //
	int arrayIndex = crossOverPoint / pop->bitSizeInt ;
	int locWithinArray = crossOverPoint % pop->bitSizeInt ;

	memcpy(off1->genes, par1->genes, sizeof(unsigned int)*(pop->chromoArSize)) ;
	memcpy(off2->genes, par2->genes, sizeof(unsigned int)*(pop->chromoArSize)) ;
	memcpy(off1->genes, par2->genes, sizeof(unsigned int)*(arrayIndex)) ;
	memcpy(off2->genes, par1->genes, sizeof(unsigned int)*(arrayIndex)) ;
	// At the array index, splits the array into a left/right side and swaps.//
	if (locWithinArray != 0) {
		int rShift = pop->bitSizeInt-locWithinArray ;
		unsigned int rightBits1 = (par1->genes[arrayIndex] << (rShift)) ;
		rightBits1 >>= rShift ;
		unsigned int leftBits1 = (par1->genes[arrayIndex] >> locWithinArray) ;
		leftBits1 <<= locWithinArray ;
		unsigned int rightBits2 = (par2->genes[arrayIndex] << (rShift)) ;
		rightBits2 >>= rShift ;
		unsigned int leftBits2 = (par2->genes[arrayIndex] >> locWithinArray) ;
		leftBits2 <<= locWithinArray ;
		off1->genes[arrayIndex] = rightBits2 + leftBits1 ;
		off2->genes[arrayIndex] = leftBits2 + rightBits1 ;
	}	
}		/* -----  end of function crossover  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  singleCrossOver
 *    Arguments:  Chromosome * par1 - The chromosome of the first parent.
 *                Chromosome * par2 - The chromosome of the second parent.
 *                Chromosome * off1 - The chromosome of the offspring.
 *                Population * pop - Population struct containing all info for chromosome
 *                dimensions.
 *  Description:  Crosses over two chromosomes at a random bit and overwrites to form
 *                a single new offspring.
 * =====================================================================================
 */

void singleCrossOver(Chromosome * par1, Chromosome * par2, Chromosome * off1, Population * pop) {
	// Generate to the location of the crossover. //
	int crossOverPoint = rand() % pop->chromoLength ;
	// Index of the array holding that group of 32 bits. //
	int arrayIndex = crossOverPoint / pop->bitSizeInt ;
	int locWithinArray = crossOverPoint % pop->bitSizeInt ;
	
	memcpy(off1->genes, par1->genes, sizeof(unsigned int)*(pop->chromoArSize)) ;
	memcpy(off1->genes, par2->genes, sizeof(unsigned int)*(arrayIndex)) ;
	// At the array index, splits the array into a left/right side and swaps.//
	if (locWithinArray != 0) {
		int rShift = pop->bitSizeInt-locWithinArray ;
		unsigned int leftBits1 = (par1->genes[arrayIndex] >> locWithinArray) ;
		leftBits1 <<= locWithinArray ;
		unsigned int rightBits2 = (par2->genes[arrayIndex] << (rShift)) ;
		rightBits2 >>= rShift ;
		off1->genes[arrayIndex] = rightBits2 + leftBits1 ;
	}
}		/* -----  end of function crossover  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  calcFitness
 *    Arguments:  Chromosome * chromo1 - Chromosome 1 from which evaluate the fitness.
 *                Chromosome * chromo2 - Chromosome 2 from which evaluate the fitness.
 *                Population * pop - Population containing info about chromosome dimensions.
 *                int * hist - The previous outcomes.
 *                int * fitness1 - Pointer to where the fitness is being stored for 1.
 *                int * fitness2 - Pointer to where the fitness is being stored for 2.
 *  Description:  Evaluates the fitness of the chromosomes versus each other in a game
 *                of prisoners dilemma given a memory of the two previous games outcomes.
 * =====================================================================================
 */

void calcFitness (Chromosome * chromo1, Chromosome * chromo2, Population * pop, int * hist, int * fitness1, int * fitness2) {
	// Find the relevant strategy in the genome given the history. //
	int arrayIndex = *hist / pop->bitSizeInt ;
	int locWithinArray = *hist % pop->bitSizeInt ;

	// Obtain the courses of action. //
	int action1 = (chromo1->genes[arrayIndex] >> locWithinArray) ;
	action1 &= 1 ;
	int action2 = (chromo2->genes[arrayIndex] >> locWithinArray) ;
	action2 &= 1 ;
	// Update the history by concatinating bits and chopping off the left side.//
	int newEvent = action1 << 1 ;
	newEvent |= action2 ;
	int newHist = *hist ;
	newHist <<= 2 ;
	newHist |= newEvent ;
	*hist = newHist & (pop->chromoLength-1) ;

	// Play the prisoners dilemma using an array of results to avoid branching. //
	*fitness1 += pop->prisonerRes[0][newEvent] ;
	*fitness2 += pop->prisonerRes[1][newEvent] ;
	
}		/* -----  end of function calcFitness  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  convertToBin
 *    Arguments:  unsigned int num - Number to convert.
 *                char * out - String to which we save the converted number.
 *  Description:  Converts a number to binary and saves the string in out.
 * =====================================================================================
 */

void convertToBin (unsigned int num, char * out) {
	int i = 0 ;
	unsigned int mask = (unsigned int) pow(2,sizeof(int)*8-1) ;
	for (i = 0 ; i < sizeof(int)*8 ; ++i) {
		*out = (num & mask) ? '1' : '0' ;
		num <<= 1 ;
		++out ;
	}
	*out = '\0' ;
}		/* -----  end of function convertToBin  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  printChromosome
 *    Arguments:  Chromosome * chromo - Chromosome to print.
 *                Population * pop - Population struct containing info about chromosome.
 *  Description:  Prints the genetic information in the chromosome.
 * =====================================================================================
 */

void printChromosome (Chromosome * chromo, Population * pop) {
	int i ;
	char buffer[pop->bitSizeInt] ;
	for (i = pop->chromoArSize-1 ; i >= 0 ; --i) {
		convertToBin(chromo->genes[i], buffer) ;
		printf("%s", buffer) ;
	}
	printf("\n") ;
}		/* -----  end of function printChromosome  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  freeChromosome
 *    Arguments:  Chromosome * chromo - Chromosome to free.
 *  Description:  Frees memory associated with chromosome.
 * =====================================================================================
 */

void freeChromosome (Chromosome * chromo) {
	free(chromo) ;
}		/* -----  end of function freeChromosome  ----- */
