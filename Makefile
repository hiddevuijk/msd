
test.exe: main.cpp system.h vec2.h potential.h
	g++ main.cpp -std=c++11 -Wall -Werror  -o test.exe -O3



.PHONY : clean
clean:
	rm test.exe

