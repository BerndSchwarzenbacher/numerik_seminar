NGSCXX = ${NETGENDIR}/ngscxx
FLAGS = -lngstd -lsolve -lngcomp -lngfem -linterface -lngla

all: main.o

main.o : main.cpp
	$(NGSCXX) $? $(FLAGS) -o $@

clean:
	rm main.o
