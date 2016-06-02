# compiler variables
CC = mpicc
CFLAGS = -Wall -O2
LDLIBS = -lm

# custom variables
objects = population.o main.o chromosome.o

pdPara : $(objects)
	$(CC) -o $@ $(objects) $(LDLIBS) $(CFLAGS) 

population.o : population.c population.h
	$(CC) -c $< $(CFLAGS) 

chromosome.o : chromosome.c chromosome.h
	$(CC) -c $< $(CFLAGS) 

main.o : main.c chromosome.h population.h
	$(CC) -c $< $(CFLAGS) 
# test target

.PHONY: clean
clean:
	rm -f pdPara $(objects)

debug : $(objects)
	$(CC) -o $@ $(objects) $(LDLIBS) $(CFLAGS) -g
