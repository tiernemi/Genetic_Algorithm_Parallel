#ifndef  CHROMOSOME_INC

#define  CHROMOSOME_INC

/*
 * =====================================================================================
 *
 *       Filename:  chromosome.h
 *
 *    Description:  Chromome.h contains the struct Chromosome which contains the inform-
 *                  ation about the genes of a population member in a genetic algorithm.
 *                  Information about size of chromosomes stored in population.
 *
 *        Version:  1.0
 *        Created:  02/01/16 15:43:59
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
 *         Name:  Chromosome
 *       Fields:  unsigned int * genes - Pointer to array containing 0-1 bits 
 *                reprensenting the genes.
 *  Description:  Chromosome contains the genetic information of a population member
 *                in a genetic algorithm.
 * =====================================================================================
 */

struct Population ;

typedef struct Chromosome {
	unsigned int * genes ;
} Chromosome;				/* ----------  end of struct Chromosome  ---------- */

void initialiseChromosome(Chromosome * chromo, struct Population * pop, unsigned int * arrayStart) ;
void mutate(Chromosome * chromo, double rate, struct Population * pop) ;
void doubleCrossOver(Chromosome * par1, Chromosome * par2, Chromosome * off1, Chromosome * off2, struct Population * pop) ;
void singleCrossOver(Chromosome * par1, Chromosome * par2, Chromosome * off1, struct Population * pop) ;
void calcFitness(Chromosome * chromo1, Chromosome * chromo2, struct Population * pop, int * hist, int * fitness1, int *fitness2) ;
void convertToBin (unsigned int num, char * out) ;
void printChromosome (Chromosome * chromo, struct Population * pop) ;
void freeChromosome (Chromosome * chromo) ;

#endif
