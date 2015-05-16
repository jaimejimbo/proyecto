NAME = proyecto
CC = g++
CFLAGS = -g #-Wall
SOURCES = main.cc
OBJS = main.o

$(OBJS): $(SOURCES)
	$(CC) $(CFLAGS) -c -o $@ $<

$(NAME): $(OBJS)
	$(CC) -o programa $(OBJS) -lm

all: $(NAME)

clean: 
	rm *.o
	rm *.a
	rm programa
