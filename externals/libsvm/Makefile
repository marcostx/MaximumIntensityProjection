all: obj/svm.o
	ar csr lib/libsvm.a obj/svm.o
	@echo
	@echo "libsvm.a built..."
	@echo "DONE."

obj/svm.o: src/svm.cpp
	$(CC) -I include/ -O3 -Wall  -fPIC -fopenmp -pthread -pedantic -c src/svm.cpp -o obj/svm.o

clean:
	rm -f obj/*.o; rm -f lib/*.a;
