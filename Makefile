CFLAGS = `root-config --cflags` -I/usr/local/digi
CXXFLAGS = `root-config --cflags` -I/usr/local/digi
LDFLAGS = `root-config --glibs` -L/usr/local/digi -lReadDigiData

pairbuilder : pairbuilder.cpp
	g++ -o $@ $^ `root-config --cflags --libs`

mupairbuilder : mupairbuilder.cpp
	g++ -o $@ $^ `root-config --cflags --libs`

cmbuilder : cmbuilder.cpp
	g++ -o $@ $^ `root-config --cflags --libs`
