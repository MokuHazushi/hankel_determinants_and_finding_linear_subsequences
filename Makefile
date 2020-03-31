LIB_DIR=/usr/local/lib
NTL_INC=/usr/local/include/NTL
GMPL_INC=/usr/local/include

all: clean
	g++ -Wall -Wpedantic -g -O5 -I${NTL_INC} -I${GMPL_INC} hanmat_mt.cpp -o hanmat_mt -std=c++11 -L${LIB_DIR} -lntl -lgmp -lm -lpthread

hanmat: clean
	g++ -Wall -Wpedantic -g -O5 -I${NTL_INC} -I${GMPL_INC} hanmat.cpp -o hanmat -std=c++11 -L${LIB_DIR} -lntl -lgmp -lm -lpthread

hanmat_sampling: clean
	g++ -Wall -Wpedantic -g -O5 -I${NTL_INC} -I${GMPL_INC} hanmat_sampling.cpp -o hanmat_sampling -std=c++11 -L${LIB_DIR} -lntl -lgmp -lm -lpthread

debug: clean
	g++ -Wall -Wpedantic -g -I${NTL_INC} -I${GMPL_INC} hanmat_mt.cpp -o hanmat_mt -std=c++11 -L${LIB_DIR} -lntl -lgmp -lm -lpthread
	gdb hanmat_mt

clean:
	rm -rf hanmat hanmat_sampling hanmat_mt
