CC = gcc
SOURCE = $(wildcard *.c)
OBJS = $(SOURCE:.c=.o)
HEADERS = $(wildcard *.h)

OUT = main
FLAGS = -g -Wall -O2 -pg -lm
FLAGS += -lgmp -lmpfr
FLAGS += -Wextra -Werror

all: $(OBJS)
	$(CC) $(FLAGS) -o $(OUT) $^ -lm -lmpfr -lgmp
# -lm at the end is for math.h


%.o: %.c $(HEADERS)
	$(CC) $(FLAGS) -o $@ -c $< 


clean:
	rm -f $(OBJS) $(OUT)
