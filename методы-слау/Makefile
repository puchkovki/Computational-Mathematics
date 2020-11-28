a.out: MatrixCSR.o Norms.o Vector.o Jacobi.o Seidel.o Iterative.o ConjGradient.o main.o
	g++ -std=c++17 -O3 -Wall -Wextra -pedantic -fsanitize=address,undefined -o SLAE *.o; rm -rf *.o

MatrixCSR.o: MatrixCSR.cpp
	g++ -std=c++17 -O3 -Wall -Wextra -pedantic -fsanitize=address,undefined -c -o MatrixCSR.o MatrixCSR.cpp

Norms.o: Norms.cpp
	g++ -std=c++17 -O3 -Wall -Wextra -pedantic -fsanitize=address,undefined -c -o Norms.o Norms.cpp

Vector.o: Vector.cpp
	g++ -std=c++17 -O3 -Wall -Wextra -pedantic -fsanitize=address,undefined -c -o Vector.o Vector.cpp

Jacobi.o: Jacobi.cpp
	g++ -std=c++17 -O3 -Wall -Wextra -pedantic -fsanitize=address,undefined -c -o Jacobi.o Jacobi.cpp

Seidel.o: Seidel.cpp
	g++ -std=c++17 -O3 -Wall -Wextra -pedantic -fsanitize=address,undefined -c -o Seidel.o Seidel.cpp

Iterative.o: Iterative.cpp
	g++ -std=c++17 -O3 -Wall -Wextra -pedantic -fsanitize=address,undefined -c -o Iterative.o Iterative.cpp

ConjGradient.o: ConjGradient.cpp
	g++ -std=c++17 -O3 -Wall -Wextra -pedantic -fsanitize=address,undefined -c -o ConjGradient.o ConjGradient.cpp

main.o: main.cpp
	g++ -std=c++17 -O3 -Wall -Wextra -pedantic -fsanitize=address,undefined -c -o main.o main.cpp

clean:
	rm -rf *.o SLAE