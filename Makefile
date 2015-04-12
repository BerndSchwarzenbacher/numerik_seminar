NGSCXX = ${NETGENDIR}/ngscxx
FLAGS = -lngstd -lsolve -lngcomp -lngfem -linterface -lngla
TARGETS = main init_par_for init_balancing

all: $(TARGETS)

% : %.cpp
	$(NGSCXX) $? $(FLAGS) -o $@

clean:
	rm $(TARGETS)
