#$Id: Makefile,v 1.11 2013/10/16 10:32:02 jiayi Exp $
#!/usr/bin/make

CC = gcc
CXX = g++
CFLAGS = -Wall -I../include -I/home/moon/jiayiliu/local/include
CXXFLAGS = $(CFLAGS)
LOADLIBES = -lm -lpress2ndkr -lpress2ndkrD -lcfitsio -lfftw3
LDFLAGS = -L/home/moon/jiayiliu/local/lib
RUN_DIR = ../run/

EXEC = test mknull mkbeam mkszmap genfilter filtermap szpeakfinder szpeakcmp formatcluster \
      genfilterM filtermapM mkoneszmap checkszsn

vpath %.c ../src
vpath %.h ../include
vpath %.cpp ../src

csources = test.c main_mknull.c main_mkbeam.c main_mkszmap.c main_genfilter.c \
		main_filtermap.c main_genfilterM.c main_mkszonemap.c checkszsn.c\
		main_filtermapM.c sub_clinput.c \
		sub_fft.c sub_fft2image.c sub_fftsave.c sub_fft2image.c sub_filter.c \
		sub_mkbolo.c sub_mkcmb.c sub_mkdirtmap.c sub_setup.c sub_stdroutines.c sub_tauvalue.c 

cppsources = main_szpeakfinder.cpp main_szpeakcmp.cpp aux_formatcluster.cpp sub_cluster.cpp sub_peak.cpp

all: main_mkbeam main_mknull main_mkszmap main_mkoneszmap main_genfilter main_genfilterM\
	 main_filtermap main_filtermapM main_szpeakfinder main_szpeakcmp checkszsn
	cp main_mkbeam $(RUN_DIR)mkbeam
	cp main_mknull $(RUN_DIR)mknull
	cp main_mkszmap $(RUN_DIR)mkszmap
	cp main_mkoneszmap $(RUN_DIR)mkoneszmap
	cp main_genfilter $(RUN_DIR)genfilter
	cp main_genfilterM $(RUN_DIR)genfilterM
	cp main_filtermap $(RUN_DIR)filtermap
	cp main_filtermapM $(RUN_DIR)filtermapM
	cp main_szpeakfinder $(RUN_DIR)szpeakfinder
	cp main_szpeakcmp $(RUN_DIR)szpeakcmp
	cp test $(RUN_DIR)test

.PHONY: clean

clean:
	rm -f *.o *.d
	rm -f main_* aux_*
	cd $(RUN_DIR); 	rm -f $(EXEC)

# test files
test: test.o sub_stdroutines.o sub_clinput.o sub_mkcmb.o sub_setup.o 

# main program
main_mkbeam: main_mkbeam.o sub_fftsave.o sub_setup.o sub_stdroutines.o sub_fft.o

main_mknull: main_mknull.o sub_stdroutines.o

main_mkszmap: main_mkszmap.o sub_clinput.o sub_fft.o sub_fft2image.o sub_fftsave.o \
		sub_mkbolo.o sub_mkcmb.o sub_mkdirtmap.o sub_setup.o sub_stdroutines.o

main_mkoneszmap: main_mkoneszmap.o sub_clinput.o sub_fft.o sub_fft2image.o sub_fftsave.o \
		sub_mkbolo.o sub_mkcmb.o sub_mkdirtmap.o sub_setup.o sub_stdroutines.o

main_genfilter: main_genfilter.o sub_clinput.o sub_fft.o sub_fft2image.o sub_fftsave.o sub_setup.o sub_stdroutines.o sub_tauvalue.o

main_genfilterM: main_genfilterM.o sub_clinput.o sub_fft.o sub_fft2image.o sub_fftsave.o sub_setup.o sub_stdroutines.o sub_tauvalue.o

main_filtermap: main_filtermap.o sub_clinput.o sub_filter.o sub_fft2image.o sub_fft.o sub_fftsave.o sub_setup.o sub_stdroutines.o

main_filtermapM: main_filtermapM.o sub_clinput.o sub_filter.o sub_fft2image.o sub_fft.o sub_fftsave.o sub_setup.o sub_stdroutines.o

main_szpeakfinder: main_szpeakfinder.o sub_peak.o sub_stdroutines_cpp.o
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) $(LOADLIBES) -o $@

main_szpeakcmp: main_szpeakcmp.o sub_peak.o sub_cluster.o sub_stdroutines_cpp.o
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) $(LOADLIBES) -o $@

checkszsn: checkszsn.o sub_cluster.o sub_stdroutines_cpp.o
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) $(LOADLIBES) -o $@
#	cp checkszsn /home/moon/jiayiliu/2013Jan_sim/spt-lowzeta/mocks-2bands/


main_mock: main_mock.o 
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) $(LOADLIBES) -o $@

# auxillary files
formatcluster: aux_formatcluster.o sub_cluster.o
	$(CXX) $(CXXFLAGS) $^ -o $(RUN_DIR)/formatcluster

sub_stdroutines_cpp.o: sub_stdroutines.c
	$(CXX) $(CXXFLAGS) -c $^ -o $@


%.d: %.c
	$(CC) -MM $(CFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ :,g'< $@.$$$$ > $@; \
	rm -f $@.$$$$

%.d: %.cpp
	$(CXX) -MM $(CXXFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ :,g'< $@.$$$$ > $@; \
	rm -f $@.$$$$

sinclude $(csources:%.c=%.d)
sinclude $(cppsources:%.cpp=%.d)
