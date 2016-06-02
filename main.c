#include <stdlib.h>
#include "population.h"
#include "chromosome.h"
#include <stdio.h>
#include <string.h>
#include <mpi/mpi.h>
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  pdPara
 *  Description:  Main function for genetic algorithm, randomly generates chromosomes
 *                and uses fit,select,cross,mutate steps to generate a new population
 *                which after many iterations should converge to a high fitness. This
 *                particular algorithm finds the optimal strategy of the iterative
 *                prisoners dilemma. This particular algorithm uses a master slave
 *                scheme to evaluate the fitness in parallel.
 *    Arguments:  -ps POPULATION SIZE 
 *                -ng NUMBER OF GAMES
 *                -c CROSSOVER RATE
 *                -m MUTATION RATE
 *                -ng NUMBER OF GENERATIONS
 * =====================================================================================
 */

int main ( int argc, char *argv[] ) {
	MPI_Init(&argc, &argv) ;
	int rank, nproc ;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank) ;
	MPI_Comm_size(MPI_COMM_WORLD, &nproc) ;
	srand(15071992+rank*200) ;
	if (argc != 11) {
		fprintf(stderr, "Insufficient Command Line Arguments \n");
		exit(-1) ;
	}
	int i = 0 ;
	int numDetected = 0 ;

	int popSize ;
	int numGenerations ;
	int numGames ;
	double mutationRate ;
	double crossRate ;
	for (i = 1 ; i < argc ; i+=2) {
		if (strcmp(argv[i], "-ps") == 0) {
			++numDetected ;
			sscanf(argv[i+1], "%d", &popSize) ;
		} else if(strcmp(argv[i], "-ni") == 0) {
			++numDetected ;
			sscanf(argv[i+1], "%d", &numGames) ;
		} else if(strcmp(argv[i], "-ng") == 0) {
			++numDetected ;
			sscanf(argv[i+1], "%d", &numGenerations) ;
		} else if(strcmp(argv[i], "-c") == 0) {
			++numDetected ;
			sscanf(argv[i+1], "%lf", &crossRate) ;
		} else if(strcmp(argv[i], "-m") == 0) {
			++numDetected ;
			sscanf(argv[i+1], "%lf", &mutationRate) ;
		}		
	}
	if (numDetected != 5) {
		fprintf(stderr, "Insufficient Command Line Arguments \n");
		exit(-1) ;
	}
	if (popSize < nproc-1) {
		fprintf(stderr, "Population below the number of processes, will causes halting.\n");
		exit(-2) ;
	}


	int chromoLength = 16 ; // 2 - Game memory. //
	// Generate the population for each rank, Only the master contains heap memory. //
	Population * pop = makePopulation(popSize, mutationRate, crossRate, chromoLength, numGames) ;
	// Master process. //
	if (rank == 0) {
		char * outputFile = "fitness.dat" ;
		double * results = malloc(sizeof(double)*numGenerations) ;
		// Run numgenerations amount of steps of the genetic algorithm (MASTER). //
		simulateM(pop, numGenerations, results) ;
	
		FILE * output = fopen(outputFile, "w") ;
		for (i = 0 ; i < numGenerations ; ++i) {
			fprintf(output, "%d\t%lf\n", i, results[i]);
		}
		fclose(output) ;
		free(results) ;
		freePopulation(pop) ;
	} else {
		// Open a thread that evaluates fitnesses of chromosomes. //
		simulateS(pop) ;
		free(pop) ;
	}

	MPI_Finalize() ;
	return EXIT_SUCCESS ;
}				/* ----------  end of function main  ---------- */

