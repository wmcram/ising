run: build
	./build/ising

build:
	mkdir -p build
	g++ main.cpp -o build/ising
