# INSTALLATION

## Requirements
The application requires the following libraries to be installed :
- GNU Multiple Precision Arithmetic Library, GMP 6.2.0 ([link](https://gmplib.org/))
- Number Theory Library, NTL 11.4.3 ([link](https://www.shoup.net/ntl/))

## Installation of GMP library
IMPORTANT, GMP should be installed before NTL. The library has to be configured to be thread-safe, this is done by adding the flag `--enable-alloca=malloc-notreentrant` during configuration.
Requires root access:
```
wget https://gmplib.org/download/gmp/gmp-6.2.0.tar.lz
tar -xf gmp-6.2.0
cd gmp-6.2.0
./configure --enable-alloca=malloc-notreentrant
make
make check
sudo make install
```

## Installation of NTL library
NTL also needs to be configured to be thread-safe, this is done by adding the following flags `NTL_THREADS=on NTL_GMP_LIP=on` during configuration.
Requires root access:
```
wget https://www.shoup.net/ntl/ntl-11.4.3.tar.gz
tar -xf ntl-11.4.3.tar.gz
cd ntl-11.4.3.tar.gz
./configure NTL_THREADS=on NTL_GMP_LIP=on
make
make check
sudo make install
```

# USAGE
The application can be built using the provided makefile:
```
make hanmat
```

The executable `hanmat` can be run as usual C-compiled application:
```
./hanmat
```
