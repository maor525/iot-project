FLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors
LIBS = -lm

all: main.o spmat.o group.o
	gcc main.o spmat.o group.o -o cluster $(LIBS)
clean:
	rm -rf *.o cluster

main.o: main.c spmat.c group.c spmat.h group.h
	gcc $(FLAGS) -c main.c

spmat.o: spmat.c spmat.h main.c
	gcc $(FLAGS) -c spmat.c

group.o: group.c group.h main.c 
	gcc $(FLAGS) -c group.c