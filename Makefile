NGSCXX = ${NETGENDIR}/ngscxx
FLAGS = -lngstd -lsolve -lngcomp -lngfem -linterface -lngla

all: main init_par_for

main : main.cpp
	$(NGSCXX) $? $(FLAGS) -o $@

init_par_for : init_par_for.cpp
	$(NGSCXX) $? $(FLAGS) -o $@

clean:
	rm main init_par_for
