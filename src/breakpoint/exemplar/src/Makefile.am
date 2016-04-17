bin_PROGRAMS=exemplar

# GUROBI_HOME, BOOST_HOME required as environmental variables

GBASEDIR = ../../../../lib/gbase/src/

exemplar_CPPFLAGS = -I$(GBASEDIR) -I$(BOOST_HOME) -I$(GUROBI_HOME)/include
exemplar_LDFLAGS = -L$(GBASEDIR) -L$(GUROBI_HOME)/lib
exemplar_LDADD = -lgbase -lgurobi_c++ -lgurobi

exemplar_SOURCES= config.h config.cc \
				  genome.h genome.cc \
				  simulator.h simulator.cc \
				  pbase.h pbase.cc \
				  spliter.h spliter.cc \
				  ilp.h ilp.cc \
				  exemplar.h exemplar.cc \
				  mygraph.h mygraph.cc \
				  trier.h trier.cc \
				  shadj.h shadj.cc \
				  psolver.h psolver.cc \
				  decomposer.h decomposer.cc \
				  gfamily.h gfamily.cc \
				  candidate.h candidate.cc \
				  isolver.h isolver.cc \
				  bbstate.h bbstate.cc \
				  bbsearch.h bbsearch.cc \
				  main.cc
