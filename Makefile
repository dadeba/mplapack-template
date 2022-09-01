CXX=g++-12
CXXFLAGS=-O3 -fopenmp -I.

programs=Raxpy_bench__Float128 Raxpy_bench_double Raxpy_bench_gmp Raxpy_bench__Float16 \
Cgemm_bench__Float128 Cgemm_bench_double Cgemm_bench_gmp Cgemm_bench__Float16 \
Rgemm_bench__Float128 Rgemm_bench_double Rgemm_bench_gmp Rgemm_bench__Float16

all: $(programs)

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $<

Raxpy_bench__Float16: Raxpy_bench__Float16.o
	$(CXX) -o Raxpy_bench__Float16 Raxpy_bench__Float16.o

Raxpy_bench__Float128: Raxpy_bench__Float128.o
	$(CXX) -o Raxpy_bench__Float128 Raxpy_bench__Float128.o

Raxpy_bench_double: Raxpy_bench_double.o
	$(CXX) -o Raxpy_bench_double Raxpy_bench_double.o

Raxpy_bench_gmp: Raxpy_bench_gmp.o
	$(CXX) -o Raxpy_bench_gmp Raxpy_bench_gmp.o -lgmpxx -lgmp

Rgemm_bench__Float128: Rgemm_bench__Float128.o
	$(CXX) -o Rgemm_bench__Float128 Rgemm_bench__Float128.o

Rgemm_bench_double: Rgemm_bench_double.o
	$(CXX) -o Rgemm_bench_double Rgemm_bench_double.o

Rgemm_bench_gmp: Rgemm_bench_gmp.o
	$(CXX) -o Rgemm_bench_gmp Rgemm_bench_gmp.o -lgmpxx -lgmp

Rgemm_bench__Float16: Rgemm_bench__Float16.o
	$(CXX) -o Rgemm_bench__Float16 Rgemm_bench__Float16.o

Cgemm_bench__Float128: Cgemm_bench__Float128.o
	$(CXX) -o Cgemm_bench__Float128 Cgemm_bench__Float128.o

Cgemm_bench_double: Cgemm_bench_double.o
	$(CXX) -o Cgemm_bench_double Cgemm_bench_double.o

Cgemm_bench_gmp: Cgemm_bench_gmp.o
	$(CXX) -o Cgemm_bench_gmp Cgemm_bench_gmp.o -lgmpxx -lgmp

Cgemm_bench__Float16: Cgemm_bench__Float16.o
	$(CXX) -o Cgemm_bench__Float16 Cgemm_bench__Float16.o

clean:
	rm -rf *.o *~ $(programs) *bak
