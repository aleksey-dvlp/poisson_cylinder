

CC = g++

OBJ = pcylind.o main.o

run: $(OBJ)
	$(CC) -o $@ $^
clean:
	rm -f run *.o
