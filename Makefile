CXX=c++
CXXFLAGS=-std=c++11 -Wall -W -Wunused -O3 -I. 
#CXXFLAGS=-std=c++11 -Wall -W  -g -I.

all: rbam.o test.o rbamC.o bamreader.o bamrecord.o  uptr.o  reader.o
	$(CXX) $(CXXFLAGS) -o rbam rbam.o -lz
	$(CXX) $(CXXFLAGS) -o rbamC rbamC.o  bamreader.o bamrecord.o  -lz -lsequence -lhts
	$(CXX) $(CXXFLAGS) -o test test.o
	$(CXX) $(CXXFLAGS) -o uptr uptr.o
	$(CXX) $(CXXFLAGS) -o reader reader.o -lsequence

clean: 
	rm -f *.o test rbam rbamC

bamreader.o: bamreader.hpp bamrecord.hpp
bamrecord.o: bamrecord.hpp bamutil.hpp
rbamC.o: bamreader.hpp bamrecord.hpp bamutil.hpp
