bin_PROGRAMS = segdcj

segdcj_CPPFLAGS = -I$(BOOST_HOME) -I$(GUROBI_HOME)/include
segdcj_LDFLAGS = -L$(GUROBI_HOME)/lib
segdcj_LDADD = -lgurobi_c++ -lgurobi

segdcj_SOURCES = main.cc \
				 trier.h trier.cc \
				 mygraph.h mygraph.cc \
				 simulator.h simulator.cc \
				 config.h config.cc \
				 genome.h genome.cc \
				 draw.h draw.cc \
				 shseg.h shseg.cc \
				 pbase.h pbase.cc \
				 asolver.h asolver.cc \
				 bsolver.h bsolver.cc \
				 psolver.h psolver.cc \
				 lpsolver.h lpsolver.cc \
				 lpsolver2.h lpsolver2.cc \
				 adjacency.h adjacency.cc \
				 chrm.h chrm.cc \
				 circular_chrm.h circular_chrm.cc \
				 common.h common.cc \
				 extremity.h extremity.cc \
				 gene.h gene.cc \
				 genome_base.h genome_base.cc \
				 item.h item.cc \
				 linear_chrm.h linear_chrm.cc \
				 operation.h operation.cc
