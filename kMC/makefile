OBJS = AG.o
INC = /usr/local/include
DIR = /usr/local/lib
COMPILER = g++
FLAGS = -std=c++11 -g -lm
AG : $(OBJS)
	$(COMPILER) -L$(DIR) $(FLAGS) $(OBJS) -o AG

$(OBJS): AG.cpp
	$(COMPILER) -I$(INC) -c AG.cpp

clean:
	rm *~ *.o *.dat AG
