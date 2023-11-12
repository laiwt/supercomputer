compile:
	g++ -O3 -std=c++11 -fopenmp main.cpp -o main

submit:
	./submit.sh

graph:
	python graph.py ./data/analytical.csv ./data/analytical.png
	python graph.py ./data/result.csv ./data/result.png
	python graph.py ./data/error.csv ./data/error.png

clear:
	rm ./*.txt