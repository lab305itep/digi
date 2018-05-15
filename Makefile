CFLAGS = `root-config --cflags` -I/usr/local/digi -g
CXXFLAGS = `root-config --cflags` -I/usr/local/digi
LDFLAGS = `root-config --glibs` -L/usr/local/digi -lReadDigiData

pairbuilder : pairbuilder.cpp
	g++ -o $@ $^ `root-config --cflags --libs`

pairbuilder2 : pairbuilder2.cpp
	g++ -o $@ $^ `root-config --cflags --libs`

pairbuilder3 : pairbuilder3.cpp
	g++ -o $@ $^ `root-config --cflags --libs`

pairbuilder4 : pairbuilder4.cpp
	g++ -o $@ $^ `root-config --cflags --libs`

pairbuilder5 : pairbuilder5.cpp
	g++ -o $@ $^ `root-config --cflags --libs`

randombuilder : randombuilder.cpp
	g++ -o $@ $^ `root-config --cflags --libs`

randombuilder3 : randombuilder3.cpp
	g++ -o $@ $^ `root-config --cflags --libs`

randombuilder4 : randombuilder4.cpp
	g++ -o $@ $^ `root-config --cflags --libs`

mupairbuilder : mupairbuilder.cpp
	g++ -o $@ $^ `root-config --cflags --libs`

cmbuilder : cmbuilder.cpp
	g++ -o $@ $^ `root-config --cflags --libs`

mu12B : mu12B.cpp
	g++ -o $@ $^ `root-config --cflags --libs`
