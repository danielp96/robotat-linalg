

CC = gcc


all: matf32 linsolve control quadprog

matf32:
	$(CC) -c matf32*.c math_util.c -lm

linsolve:
	$(CC) -c linsolve.c

quadprog:
	$(CC) -c quadprog.c

control:
	$(CC) -c robotat_control.c


clean:
	rm -f *.o


.PHONY: clean matf32 linsolve quadprog