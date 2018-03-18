CC=g++
CFLAGS="-Ofast"


cppCartoon: main.o cartoon.o cpdb/cpdb.o cpdb/utils.o
	$(CC) -o cppCartoon main.o cartoon.o cpdb/utils.o cpdb/cpdb.o 


# %.o: %.cpp
# 	$(CC) -c -o $@ $< $(CFLAGS)

cpdb/utils.o: cpdb/utils.cpp
	$(CC) -c -o cpdb/utils.o cpdb/utils.cpp $(CFLAGS)

cpdb/cpdb.o: cpdb/cpdb.cpp
	$(CC) -c -o cpdb/cpdb.o cpdb/cpdb.cpp $(CFLAGS)

main.o: main.cpp
	$(CC) -c -o main.o main.cpp $(CFLAGS)

cartoon.o: cartoon.cpp
	$(CC) -c -o cartoon.o cartoon.cpp $(CFLAGS)


clean:
	rm -f *.o cppCartoon
	rm -f cpdb/*.o
