bin_PROGRAMS = dcj

dcj_CPPFLAGS = -I$(BOOST_HOME) -I$(GUROBI_HOME)/include
dcj_LDFLAGS = -L$(GUROBI_HOME)/lib
dcj_LDADD = -lgurobi_c++ -lgurobi

dcj_SOURCES = gene.h gene.cc \
			  operation.h operation.cc \
			  chrm.h chrm.cc \
			  linear_chrm.h linear_chrm.cc \
			  circular_chrm.h circular_chrm.cc \
			  common.h common.cc \
			  genome_base.h genome_base.cc \
			  extremity.h extremity.cc \
			  adjacency.h adjacency.cc \
			  config.h config.cc \
			  genome.h genome.cc \
			  presolver.h presolver.cc \
			  ilp_base.h ilp_base.cc \
			  ilp6.h ilp6.cc \
			  distoptimizer.h distoptimizer.cc \
			  mygraph.h mygraph.cc \
			  gredo.h gredo.cc \
			  main.cc
