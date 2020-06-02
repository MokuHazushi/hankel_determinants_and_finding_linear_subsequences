LIB_DIR=/usr/local/lib
NTL_INC=/usr/local/include/NTL
GMPL_INC=/usr/local/include

lib:
	g++ -Wall -Wpedantic -c -g -O5 -Iinclude src/detfct.cpp -std=c++11 
	ar rvs libdetfct.a detfct.o
	rm -rf detfct.o

test: lib
	g++ -Wall -Wpedantic -g -O5 -Iinclude -I${NTL_INC} -I${GMPL_INC} test/detfct_test.cpp -o detfct_test -L${LIB_DIR} -L. -ldetfct -lntl -lgmp -lm -lpthread -std=c++11
	./detfct_test


clean:
	rm -rf libdetfct.a detfct_test
