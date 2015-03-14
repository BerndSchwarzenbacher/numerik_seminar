NGSCXX = ${NETGENDIR}/ngscxx -lngstd -lsolve -lngcomp -lngfem -linterface -lngla

all: main.o

main.o : main.cpp
	$(NGSCXX) $? -o $@

