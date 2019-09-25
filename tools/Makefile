run: tests.exe
	./tests.exe
	rm test_file.txt

tests.exe: tests.cpp catch-main.o tools*
	g++ -std=c++14 catch-main.o tests.cpp -o tests.exe -lgsl -lgslcblas -lm

catch-main.o: catch-main.cpp
	g++ -std=c++14 catch-main.cpp -c


