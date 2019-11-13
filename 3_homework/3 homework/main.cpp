// Main interface
#include <fstream>              // File I/O operations
#include <sstream>              // String stream
#include <iomanip>              // Output formatting
#include <csignal>              // Signal handling
#include <unistd.h>             // Multiprocessing
#include <time.h>               // Timestamps

#include "General.hpp"          // General macroses
#include "MatrixCSR.hpp"        // MatrixCSR class - compressed sparce row matrices
#include "Vector.hpp"           // Vector class - standard multidimension linear algebra vector
#include "Norms.hpp"            // Second norm of vector
#include "Jacobi.hpp"           // Jacobi method for SLAE
#include "Seidel.hpp"           // Seidel method for SLAE
#include "Iterative.hpp"        // Simple single-parameter iterative method for SLAE
#include "ConjGradient.hpp"     // Conjugate gradient method for SLAE

using namespace std;

ofstream logger;

// Saves message to_print to logfile with timestamp
void log(ofstream &logger, string to_print) {

    time_t now = time(0);
    tm * local = localtime(&now);
    logger << put_time(local, "%F %T") << ": " << to_print << endl;
}

// Hendles the system interrupt signal
void signal_handler(int signal) {

    stringstream to_print;
    to_print << "Stopped by signal: " << signal << ".";
    log(logger, to_print.str());
    exit(1);
}

int main(void) {
    
    // Setting up signal handler
    signal(SIGINT, signal_handler);

    // Preparing logfile
    time_t now = time(0);
    tm * local = localtime(&now);

    stringstream logname;
    logname << "logs/SLAE_";
    logname << put_time(local, "%F") << "_";
    logname << put_time(local, "%T");
    logname << ".log";
    logger = ofstream(logname.str());

    logger << put_time(local, "%F %T") << ": Started;" << endl;

    log(logger, "Vectors b are ready");

    // Output vector
    Vector x;

    // Creating a process for graph creating
    pid_t pid = fork();
    if (pid == 0) {
        log(logger, "Graphs started");

        // Setting commas instead of dots in float numbers
        locale mylocale("");

        /* Creating output file streams for graphs */

        ofstream JacobiGraphOut;
        JacobiGraphOut.imbue(mylocale);
        JacobiGraphOut.precision(17);

        ofstream SeidelGraphOut;
        SeidelGraphOut.imbue(mylocale);
        SeidelGraphOut.precision(17);

        ofstream IterGraphOut;
        IterGraphOut.imbue(mylocale);
        IterGraphOut.precision(17);

        ofstream ConjGradOut;
        ConjGradOut.imbue(mylocale);
        ConjGradOut.precision(17);

        log(logger, "Graphs ready");

        /* n = 5 calculating */
        
        // Parsing matrix file
        MatrixCSR a_n5("matrix/matrix1_n5.mycsr");

        // Preparing right-side parts as x_i = 1 / (i + 1) for i = 0 to n - 1, where n is dimension
        Vector b_n5(a_n5.get_Dim(), 0.0);
        for (size_t i = 0; i < 5; ++i) {
            b_n5[i] = 1 / ((double) i + 1);
        }

        // Jacobi method
        JacobiGraphOut.open("tables/Jacobi_n5.csv");
        x = JacobiGraph(a_n5, b_n5, Norm2, EPS);
        for (size_t i = 0; i < x.size(); ++i) {
            JacobiGraphOut << i << ";" << x[i] << endl;
        }
        JacobiGraphOut.close();
        log(logger, "Graphs: Jacobi_n5 is done;");

        // Seidel method
        SeidelGraphOut.open("tables/Seidel_n5.csv");
        x = SeidelGraph(a_n5, b_n5, Norm2, EPS);
        for (size_t i = 0; i < x.size(); ++i) {
            SeidelGraphOut << i << ";" << x[i] << endl;
        }
        SeidelGraphOut.close();
        log(logger, "Graphs: Seidel_n5 is done;");

        // Simple single-parameter iterative method
        IterGraphOut.open("tables/Iterative_n5.csv");
        x = IterativeGraph(a_n5, b_n5, Norm2, EPS);
        for (size_t i = 0; i < x.size(); ++i) {
            IterGraphOut << i << ";" << x[i] << endl;
        }
        IterGraphOut.close();
        log(logger, "Graphs: Iterative_n5 is done;");

        // Conjugate gradient method
        ConjGradOut.open("tables/ConjugateGradient_n5.csv");
        x = ConjGradientGraph(a_n5, b_n5, Norm2, EPS);
        for (size_t i = 0; i < x.size(); ++i) {
            ConjGradOut << i << ";" << x[i] << endl;
        }
        ConjGradOut.close();
        log(logger, "Graphs: ConjugateGradient_n5 is done;");

        a_n5.~MatrixCSR();

        /* n = 100 calculating */
        
        // Parsing matrix file
        MatrixCSR a_n100("matrix/matrix1_n100.mycsr");

        // Preparing right-side parts as x_i = 1 / (i + 1) for i = 0 to n - 1, where n is dimension
        Vector b_n100(a_n100.get_Dim(), 0.0);
        for (size_t i = 0; i < 100; ++i) {
            b_n100[i] = 1 / ((double) i + 1);
        }

        // Jacobi method
        JacobiGraphOut.open("tables/Jacobi_n100.csv");
        x = JacobiGraph(a_n100, b_n100, Norm2, EPS);
        for (size_t i = 0; i < x.size(); ++i) {
            JacobiGraphOut << i << ";" << x[i] << endl;
        }
        JacobiGraphOut.close();
        log(logger, "Graphs: Jacobi_n100 is done;");

        // Seidel method
        SeidelGraphOut.open("tables/Seidel_n100.csv");
        x = SeidelGraph(a_n100, b_n100, Norm2, EPS);
        for (size_t i = 0; i < x.size(); ++i) {
            SeidelGraphOut << i << ";" << x[i] << endl;
        }
        SeidelGraphOut.close();
        log(logger, "Graphs: Seidel_n100 is done;");

        // Simple single-parameter iterative method
        IterGraphOut.open("tables/Iterative_n100.csv");
        x = IterativeGraph(a_n100, b_n100, Norm2, EPS);
        for (size_t i = 0; i < x.size(); ++i) {
            IterGraphOut << i << ";" << x[i] << endl;
        }
        IterGraphOut.close();
        log(logger, "Graphs: Iterative_n100 is done;");

        // Conjugate gradient method
        ConjGradOut.open("tables/ConjugateGradient_n100.csv");
        x = ConjGradientGraph(a_n100, b_n100, Norm2, EPS);
        for (size_t i = 0; i < x.size(); ++i) {
            ConjGradOut << i << ";" << x[i] << endl;
        }
        ConjGradOut.close();
        log(logger, "Graphs: ConjugateGradient_n100 is done;");

        a_n100.~MatrixCSR();

        /* n = 1000 calculating */
        
        // Parsing matrix file
        MatrixCSR a_n1000("matrix/matrix1_n1000.mycsr");

        // Preparing right-side parts as x_i = 1 / (i + 1) for i = 0 to n - 1, where n is dimension
        Vector b_n1000(a_n1000.get_Dim(), 0.0);
        for (size_t i = 0; i < 1000; ++i) {
            b_n1000[i] = 1 / ((double) i + 1);
        }

        // Jacobi method
        JacobiGraphOut.open("tables/Jacobi_n1000.csv");
        x = JacobiGraph(a_n1000, b_n1000, Norm2, EPS);
        for (size_t i = 0; i < x.size(); ++i) {
            JacobiGraphOut << i << ";" << x[i] << endl;
        }
        JacobiGraphOut.close();
        log(logger, "Graphs: Jacobi_n1000 is done;");

        // Seidel method
        SeidelGraphOut.open("tables/Seidel_n1000.csv");
        x = SeidelGraph(a_n1000, b_n1000, Norm2, EPS);
        for (size_t i = 0; i < x.size(); ++i) {
            SeidelGraphOut << i << ";" << x[i] << endl;
        }
        SeidelGraphOut.close();
        log(logger, "Graphs: Seidel_n1000 is done;");

        // Simple single-parameter iterative method
        IterGraphOut.open("tables/Iterative_n1000.csv");
        x = IterativeGraph(a_n1000, b_n1000, Norm2, EPS);
        for (size_t i = 0; i < x.size(); ++i) {
            IterGraphOut << i << ";" << x[i] << endl;
        }
        IterGraphOut.close();
        log(logger, "Graphs: Iterative_n1000 is done;");

        // Conjugate gradient method
        ConjGradOut.open("tables/ConjugateGradient_n1000.csv");
        x = ConjGradientGraph(a_n1000, b_n1000, Norm2, EPS);
        for (size_t i = 0; i < x.size(); ++i) {
            ConjGradOut << i << ";" << x[i] << endl;
        }
        ConjGradOut.close();
        log(logger, "Graphs: ConjugateGradient_n1000 is done;");

        a_n1000.~MatrixCSR();

        /* n = 10000 calculating */
        
        // Parsing matrix file
        MatrixCSR a_n10000("matrix/matrix1_n10000.mycsr");

        // Preparing right-side parts as x_i = 1 / (i + 1) for i = 0 to n - 1, where n is dimension
        Vector b_n10000(a_n10000.get_Dim(), 0.0);
        for (size_t i = 0; i < 10000; ++i) {
            b_n10000[i] = 1 / ((double) i + 1);
        }

        // Jacobi method
        JacobiGraphOut.open("tables/Jacobi_n10000.csv");
        x = JacobiGraph(a_n10000, b_n10000, Norm2, EPS);
        for (size_t i = 0; i < x.size(); ++i) {
            JacobiGraphOut << i << ";" << x[i] << endl;
        }
        JacobiGraphOut.close();
        log(logger, "Graphs: Jacobi_n10000 is done;");

        // Seidel method
        SeidelGraphOut.open("tables/Seidel_n10000.csv");
        x = SeidelGraph(a_n10000, b_n10000, Norm2, EPS);
        for (size_t i = 0; i < x.size(); ++i) {
            SeidelGraphOut << i << ";" << x[i] << endl;
        }
        SeidelGraphOut.close();
        log(logger, "Graphs: Seidel_n10000 is done;");

        // Simple single-parameter iterative method
        IterGraphOut.open("tables/Iterative_n10000.csv");
        x = IterativeGraph(a_n10000, b_n10000, Norm2, EPS);
        for (size_t i = 0; i < x.size(); ++i) {
            IterGraphOut << i << ";" << x[i] << endl;
        }
        IterGraphOut.close();
        log(logger, "Graphs: Iterative_n10000 is done;");

        // Conjugate gradient method
        ConjGradOut.open("tables/ConjugateGradient_n10000.csv");
        x = ConjGradientGraph(a_n10000, b_n10000, Norm2, EPS);
        for (size_t i = 0; i < x.size(); ++i) {
            ConjGradOut << i << ";" << x[i] << endl;
        }
        ConjGradOut.close();
        log(logger, "Graphs: ConjugateGradient_n10000 is done;");

        a_n10000.~MatrixCSR();
        return 0;
    }

    log(logger, "Solving started");

    // Creating outout file stream for calculating results
    ofstream Result;
    Result.precision(15);

    /* n = 5 calculating */
        
    // Parsing matrix file
    MatrixCSR a_n5("matrix/matrix1_n5.mycsr");

    // Preparing right-side parts as x_i = 1 / (i + 1) for i = 0 to n - 1, where n is dimension
    Vector b_n5(a_n5.get_Dim(), 0.0);
    for (size_t i = 0; i < 5; ++i) {
        b_n5[i] = 1 / ((double) i + 1);
    }

    log(logger, "Solving: n_5 started");
    
    Result.open("results/result_n5.txt");

    Result.width(20);
    Result << "b: ";
    Result << b_n5 << endl;

    // Jacobi method
    Result.width(20);
    Result << "Jacobi: ";
    x = Jacobi(a_n5, b_n5, Norm2, EPS);
    Result << x << endl;
    log(logger, "Solving: Jacobi is done;");

    // Seidel method
    Result.width(20);
    Result << "Seidel: ";
    x = Seidel(a_n5, b_n5, Norm2, EPS);
    Result << x << endl;
    log(logger, "Solving: Seidel is done;");

    // Simple single-parameter iterative method
    Result.width(20);
    Result << "Iterative: ";
    x = Iterative(a_n5, b_n5, Norm2, EPS);
    Result << x << endl;
    log(logger, "Solving: Iterative is done;");

    // Conjugate gradient method
    Result.width(20);
    Result << "Conjugate gradient: ";
    x = ConjGradient(a_n5, b_n5, Norm2, EPS);
    Result << x << endl;
    log(logger, "Solving: ConjugateGradient is done;");

    Result.close();

    a_n5.~MatrixCSR();

    /* n = 100 calculating */
        
    // Parsing matrix file
    MatrixCSR a_n100("matrix/matrix1_n100.mycsr");

    // Preparing right-side parts as x_i = 1 / (i + 1) for i = 0 to n - 1, where n is dimension
    Vector b_n100(a_n100.get_Dim(), 0.0);
    for (size_t i = 0; i < 100; ++i) {
        b_n100[i] = 1 / ((double) i + 1);
    }

    log(logger, "Solving: n_100 started");

    Result.open("results/result_n100.txt");

    Result.width(20);
    Result << "b: ";
    Result << b_n100 << endl;

    // Jacobi method
    Result.width(20);
    Result << "Jacobi: ";
    x = Jacobi(a_n100, b_n100, Norm2, EPS);
    Result << x << endl;
    log(logger, "Solving: Jacobi is done;");

    // Seidel method
    Result.width(20);
    Result << "Seidel: ";
    x = Seidel(a_n100, b_n100, Norm2, EPS);
    Result << x << endl;
    log(logger, "Solving: Seidel is done;");

    // Simple single-parameter iterative method
    Result.width(20);
    Result << "Iterative: ";
    x = Iterative(a_n100, b_n100, Norm2, EPS);
    Result << x << endl;
    log(logger, "Solving: Iterative is done;");

    // Conjugate gradient method
    Result.width(20);
    Result << "Conjugate gradient: ";
    x = ConjGradient(a_n100, b_n100, Norm2, EPS);
    Result << x << endl;
    log(logger, "Solving: ConjugateGradient is done;");

    Result.close();

    a_n100.~MatrixCSR();

    /* n = 1000 calculating */
        
    // Parsing matrix file
    MatrixCSR a_n1000("matrix/matrix1_n1000.mycsr");

    // Preparing right-side parts as x_i = 1 / (i + 1) for i = 0 to n - 1, where n is dimension
    Vector b_n1000(a_n1000.get_Dim(), 0.0);
    for (size_t i = 0; i < 1000; ++i) {
        b_n1000[i] = 1 / ((double) i + 1);
    }
    
    log(logger, "Solving: n_1000 started");

    Result.open("results/result_n1000.txt");

    Result.width(20);
    Result << "b: ";
    Result << b_n1000 << endl;

    // Jacobi method
    Result.width(20);
    Result << "Jacobi: ";
    x = Jacobi(a_n1000, b_n1000, Norm2, EPS);
    Result << x << endl;
    log(logger, "Solving: Jacobi is done;");

    // Seidel method
    Result.width(20);
    Result << "Seidel: ";
    x = Seidel(a_n1000, b_n1000, Norm2, EPS);
    Result << x << endl;
    log(logger, "Solving: Seidel is done;");

    // Simple single-parameter iterative method
    Result.width(20);
    Result << "Iterative: ";
    x = Iterative(a_n1000, b_n1000, Norm2, EPS);
    Result << x << endl;
    log(logger, "Solving: Iterative is done;");

    // Conjugate gradient method
    Result.width(20);
    Result << "Conjugate gradient: ";
    x = ConjGradient(a_n1000, b_n1000, Norm2, EPS);
    Result << x << endl;
    log(logger, "Solving: ConjugateGradient is done;");

    Result.close();

    a_n1000.~MatrixCSR();

    /* n = 10000 calculating */
        
    // Parsing matrix file
    MatrixCSR a_n10000("matrix/matrix1_n10000.mycsr");

    // Preparing right-side parts as x_i = 1 / (i + 1) for i = 0 to n - 1, where n is dimension
    Vector b_n10000(a_n10000.get_Dim(), 0.0);
    for (size_t i = 0; i < 10000; ++i) {
        b_n10000[i] = 1 / ((double) i + 1);
    }

    log(logger, "Solving: n_10000 started");

    Result.open("results/result_n10000.txt");

    Result.width(20);
    Result << "b: ";
    Result << b_n10000 << endl;

    // Jacobi method
    Result.width(20);
    Result << "Jacobi: ";
    x = Jacobi(a_n10000, b_n10000, Norm2, EPS);
    Result << x << endl;
    log(logger, "Solving: Jacobi is done;");

    // Seidel method
    Result.width(20);
    Result << "Seidel: ";
    x = Seidel(a_n10000, b_n10000, Norm2, EPS);
    Result << x << endl;
    log(logger, "Solving: Seidel is done;");

    // Simple single-parameter iterative method
    Result.width(20);
    Result << "Iterative: ";
    x = Iterative(a_n10000, b_n10000, Norm2, EPS);
    Result << x << endl;
    log(logger, "Solving: Iterative is done;");

    // Conjugate gradient method
    Result.width(20);
    Result << "Conjugate gradient: ";
    x = ConjGradient(a_n10000, b_n10000, Norm2, EPS);
    Result << x << endl;
    log(logger, "Solving: ConjugateGradient is done;");

    Result.close();

    a_n10000.~MatrixCSR();

    // Waiting for child process
    int res;
    wait(&res);

    log(logger, "Done.");
    
    return 0;
}
