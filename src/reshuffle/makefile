all:
	g++ -O3 -Wall -c -fmessage-length=0 -o iout_file.o iout_file.cpp
	g++ -O3 -Wall -c -fmessage-length=0 -o reshuffle_obj.o reshuffle.cpp
	g++ -O3 -Wall -c -fmessage-length=0 -o Parameters.o Parameters.cpp
	g++ -O3 -Wall -c -fmessage-length=0 -o main.o main.cpp
	g++ -o reshuffle reshuffle_obj.o main.o iout_file.o Parameters.o
	
clean:

	rm -f iout_file.o reshuffle_obj.o Parameters.o main.o reshuffle