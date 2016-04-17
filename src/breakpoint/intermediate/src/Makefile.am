bin_PROGRAMS=intermediate

# GUROBI_HOME, BOOST_HOME required as environmental variables

GBASEDIR = ../../../../lib/gbase/src/

intermediate_CPPFLAGS = -I$(GBASEDIR) -I$(BOOST_HOME) -I$(GUROBI_HOME)/include
intermediate_LDFLAGS = -L$(GBASEDIR) -L$(GUROBI_HOME)/lib
intermediate_LDADD = -lgbase -lgurobi_c++ -lgurobi

intermediate_SOURCES= config.h config.cc \
					  common.h common.cc \
					  genome.h genome.cc \
					  simulator.h simulator.cc \
					  pbase.h pbase.cc \
					  candidate.h candidate.cc \
					  shadj.h shadj.cc \
					  ilp.h ilp.cc \
					  pbase0.h pbase0.cc \
					  ilp0.h ilp0.cc \
					  mygraph.h mygraph.cc \
					  psolver.h psolver.cc \
					  verifier.h verifier.cc \
					  draw.h draw.cc \
					  intermediate.h intermediate.cc \
					  trier.h trier.cc \
					  main.cc
