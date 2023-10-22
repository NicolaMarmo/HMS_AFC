# Eigen

## Installazione
Eigen is a set of open-source C++ libraries for linear algebra operations. 

The latest stable release can be downloaded from the official website (http://eigen.tuxfamily.org/). If the tar.bz2 archive is downloaded, the files can be extracted as

    $ sudo tar jxvf eigen-X.Y.Z.tar.bz2

To install Eigen to a default location, move into the new folder that contains all the extracted files by running

    $ cd eigen-X.Y.Z

then run the following commands

    $ sudo mkdir build
    $ cd build
    $ sudo cmake ..
    $ sudo make install

## Comdandi utili
Accesso agli elementi di un vettore `x.segment(i,n)  -->  x(i:i+n)`

    // x = [1, 2, 3, 4, 5, 6]
    x.segment(0,3)      // ans = {1, 2, 3}



Linspace-equivalent:    `VectorXd::LinSpaced(size, low, high)`    



    

