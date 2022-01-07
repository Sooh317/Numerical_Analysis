main: main.cpp distr/load.cpp distr/stencil5.cpp
	g++ -O2 main.cpp distr/load.cpp distr/stencil5.cpp

image:
	./a.out >| res.txt

clean:
	rm -f a.out res.txt res.png
