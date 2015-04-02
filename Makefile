NGSCXX = ${NETGENDIR}/ngscxx
FLAGS = -lngstd -lsolve -lngcomp -lngfem -linterface -lngla

all: main

main : main.cpp
	$(NGSCXX) $? $(FLAGS) -o $@

clean:
	rm main
