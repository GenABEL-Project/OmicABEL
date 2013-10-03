all:
	g++ -O3 -Wall -c -fmessage-length=0 -o iout_file.o "D:\Work\SVNRepository\OmicABEL\src\reshuffle\iout_file.cpp" 
	g++ -O3 -Wall -c -fmessage-length=0 -o reshuffle.o "D:\Work\SVNRepository\OmicABEL\src\reshuffle\reshuffle.cpp"  
	g++ -O3 -Wall -c -fmessage-length=0 -o Parameters.o "D:\Work\SVNRepository\OmicABEL\src\reshuffle\Parameters.cpp" 
	g++ -O3 -Wall -c -fmessage-length=0 -o main.o "D:\Work\SVNRepository\OmicABEL\src\reshuffle\main.cpp" 
	g++ -o Reshuffle.exe "reshuffle.o" "main.o" "iout_file.o" "Parameters.o"

clean:

	del iout_file.o
	del reshuffle.o
	del Parameters.o
	del main.o
	del Reshuffle.exe