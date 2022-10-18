CC = gcc
SOURCE = $(wildcard *.c)
OBJS = $(SOURCE:.c=.o)
HEADERS = $(wildcard *.h)

OUT = main

FLAGS = -g -Wall -O2 -pg -lgmp -lmpfr

all: $(OBJS)
	$(CC) $(FLAGS) -o $(OUT) $^

%.o: %.c $(HEADERS)
	$(CC) $(FLAGS) -o $@ -c $< 

clean:
	rm -f $(OBJS) $(OUT)
