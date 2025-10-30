.PHONY: run
run: build
	./target/ising

.PHONY: build
build:
	mkdir -p target
	g++ main.cpp -o target/ising
