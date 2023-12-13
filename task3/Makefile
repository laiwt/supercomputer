compile_mpi_polus:
	g++ -O3 -std=c++11 -fopenmp main_mpi.cpp -o main_mpi -I/opt/ibm/spectrum_mpi/include -L/opt/ibm/spectrum_mpi/lib -lmpiprofilesupport -lmpi_ibm -DUSE_MPI

compile_mpi_local:
	g++ -O3 -std=c++11 main_mpi.cpp -o main_mpi -fopenmp -l msmpi -L "C:\Program Files (x86)\Microsoft SDKs\MPI\Lib\x86" -I "C:\Program Files (x86)\Microsoft SDKs\MPI\Include"

compile:
	g++ -O3 -std=c++11 -fopenmp main.cpp -o main

submit_mpi:
	./submit_mpi.sh

submit:
	./submit.sh

graph:
	python graph.py ./data/analytical.csv ./data/analytical.png
	python graph.py ./data/result.csv ./data/result.png
	python graph.py ./data/error.csv ./data/error.png

clear:
	rm ./*.txt