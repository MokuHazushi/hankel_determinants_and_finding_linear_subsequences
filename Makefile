LIB_DIR=/usr/local/lib
NTL_INC=/usr/local/include/NTL
GMPL_INC=/usr/local/include

hanmat:
	g++ -Wall -Wpedantic -g -O5 -I${NTL_INC} -I${GMPL_INC} hanmat.cpp -o hanmat -std=c++11 -L${LIB_DIR} -lntl -lgmp -lm -lpthread

hanmat_sampling:
	g++ -Wall -Wpedantic -g -O5 -I${NTL_INC} -I${GMPL_INC} hanmat_sampling.cpp -o hanmat_sampling -std=c++11 -L${LIB_DIR} -lntl -lgmp -lm -lpthread

clean:
	rm -rf hanmat
	rm -rf hanmat_sampling
