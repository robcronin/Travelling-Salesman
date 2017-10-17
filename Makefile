CC = gcc
CFLAGS = -W -Wall
LDFLAGS = -lm -g

objects = tsp.c

.PHONY: all test clean fullclean
all: open.out closed.out plot.out

open.out: $(objects)
	$(CC) $(objects) -o $@ $(LDFLAGS)

closed.out: $(objects)
	$(CC) $(objects) -o $@ -DCLOSED $(LDFLAGS)

plot.out: $(objects)
	$(CC) $(objects) -o $@ -DCLOSED -DPLOT $(LDFLAGS)

test: open.out closed.out path.gnu time.gnu
	./closed.out -f inputs/berlin52.tsp
	gnuplot path.gnu
	./open.out -b -t 11 -o timing.dat
	gnuplot time.gnu

path: path.gnu
	gnuplot path.gnu

time: time.gnu
	gnuplot time.gnu

clean :
	$(RM) *.o open.out closed.out plot.out

fullclean: clean
	$(RM) *.dat *.png *.ps
