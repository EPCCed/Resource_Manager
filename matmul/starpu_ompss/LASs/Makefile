include make.inc

CC		 =	mcc
SRC		 =	src
BIN		 =	bin
INC		 +=	-Iinclude
#FLAGS
#FLAGS   = -O3 -Wall --ompss --instrumentation 
FLAGS    += -O3 -Wall --ompss-2

all: liblass.a

liblass.a:
	mkdir -p $(BIN)

	#TILE TEST
	$(CC) $(INC) $(FLAGS) -c $(SRC)/ddss_tile.c -o $(BIN)/ddss_tile.o $(LIBS)
	$(CC) $(INC) $(FLAGS) -c $(SRC)/ddss_flat2tiled.c -o \
          $(BIN)/ddss_flat2tiled.o $(LIBS)
	$(CC) $(INC) $(FLAGS) -c $(SRC)/ddss_tiled2flat.c -o \
          $(BIN)/ddss_tiled2flat.o $(LIBS)
	
	#DGEMM
	$(CC) $(INC) $(FLAGS) -c $(SRC)/ddss_dgemm.c -o \
	      $(BIN)/ddss_dgemm.o $(LIBS)  
	$(CC) $(INC) $(FLAGS) -c $(SRC)/kdgemm.c -o $(BIN)/kdgemm.o $(LIBS) 

	# DGEMM BIN
	$(AR) rcs liblass.a $(BIN)/ddss_tile.o \
          $(BIN)/ddss_flat2tiled.o $(BIN)/ddss_tiled2flat.o \
		  $(BIN)/ddss_dgemm.o $(BIN)/kdgemm.o

clean:
	rm -rf $(BIN)
	rm liblass.a
