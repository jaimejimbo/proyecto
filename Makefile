all:
	clear && clear
	g++ -g main.cc -o proyecto
        ./proyecto
	gnuplot graphs.plot
