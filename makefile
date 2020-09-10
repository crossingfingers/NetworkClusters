FLAGS =
LIBS = -lm

all: cluster
clean:
	rm -rf *.o cluster

spmat.o: spmat.c
	gcc $(FLAGS) -c spmat.c
	
utils.o: utils.c
	gcc $(FLAGS) -c utils.c

divider.o: divider.c
	gcc $(FLAGS) -c divider.c
	
main.o: main.c
	gcc $(FLAGS) -c main.c

cluster: spmat.o main.o divider.o utils.o
	gcc main.o spmat.o utils.o divider.o -o cluster $(LIBS)
