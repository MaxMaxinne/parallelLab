TARGET=swsse2

CC=icc
INCLUDE=./include

FLAGS= -g -O3 -fopenmp
SOURCES = swsse2.c swstriped.c fastalib.c matrix.c
OBJECTS = $(SOURCES:.c=.o)

all:$(TARGET)
.PHONY:clean

$(shell if [ ! -e ./obj ];then mkdir -p ./obj; fi)

%.o : %.c
	$(CC) -c $< -o $@ $(FLAGS) -I$(INCLUDE)

$(TARGET) : $(OBJECTS)
	$(CC) -o $(TARGET) $^ $(FLAGS) -I$(INCLUDE)
	
clean:
	rm -f $(OBJECTS) $(TARGET)
