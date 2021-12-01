
test.exe: main.cpp system.h system2.h vec2.h potential.h msd.h ConfigFile.h\
          diffusion.h multi_diffusion.h
	g++ main.cpp -std=c++11 -Wall -Werror  -o test.exe -O3



.PHONY : clean
clean:
	rm test.exe

