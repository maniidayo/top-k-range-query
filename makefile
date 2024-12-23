CC      = g++
CFLAGS  = -O3 -mavx -std=c++17 -w
LDFLAGS =


SOURCES = utils.cpp containers/relation.cpp containers/offsets_templates.cpp containers/offsets.cpp indices/hierarchicalindex.cpp indices/hint.cpp indices/hint_m.cpp indices/hint_m_subs+sort.cpp indices/Interval.h
OBJECTS = $(SOURCES:.cpp=.o)

main: $(OBJECTS)
	$(CC) $(CFLAGS) $(LDFLAGS) utils.o containers/relation.o containers/offsets_templates.o containers/offsets.o indices/hierarchicalindex.o indices/hint_m.o indices/hint_m_subs+sort.o main.cpp -o TEST.exec $(LDADD)

main_weighted: $(OBJECTS)
	$(CC) $(CFLAGS) $(LDFLAGS) utils.o containers/relation.o containers/offsets_templates.o containers/offsets.o indices/hierarchicalindex.o indices/hint_m.o indices/hint_m_subs+sort.o main.cpp -o TEST.exec $(LDADD)


.cpp.o:
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf utils.o
	rm -rf containers/*.o
	rm -rf indices/*.o
	rm -rf TEST.exec
