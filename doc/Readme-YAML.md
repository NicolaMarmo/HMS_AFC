## Note su yaml-cpp

YAML-cpp is an open-source parser and emitter of .yaml files. These are generally used as configuration files. While they are not used in the core library, they are used in most recent examples, so it is highly recommended to install this library.

Download the latest tar.gz release from https://github.com/jbeder/yaml-cpp/releases and extract it as

    $ sudo tar zxvf yaml-cpp-yaml-cpp-X.Y.Z.tar.gz

Then, enter the new folder and run

    $ sudo mkdir build
    $ cd build
    $ sudo cmake ..
    $ sudo make install

Both the archive and the folder can be safely deleted.

