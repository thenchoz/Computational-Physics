CC = g++
CFLAGS = -Wall  -pedantic -ansi -g -std=c++11
EXEC_NAME = Ncorps_OSullivan_Henchoz.e
INCLUDES =
LIBS = -lm
OBJ_FILES = Ncorps_OSullivan_Henchoz.o

all : $(EXEC_NAME)

clean :
	rm $(EXEC_NAME) $(OBJ_FILES) *.out

$(EXEC_NAME) : $(OBJ_FILES)
	$(CC) -o $(EXEC_NAME) $(OBJ_FILES) $(LIBS)

%.o: %.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

Exercice4.o: ConfigFile.tpp
