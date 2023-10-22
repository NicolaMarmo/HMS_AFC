# Installazione
Install python
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    chmod u+x Miniconda3-latest-Linux-x86_64.sh
    ./Miniconda3-latest-Linux-x86_64.sh 
    // scegliere "lanciare conda init?" > YES!"


## Creazione di un ambiente "matematica"
    conda create -n "conda-math" python=3.8.2

    conda activate conda-math
    conda install numpy scipy
    conda install matplotlib
    
per usarlo
    conda activate conda-math


## Cose utili per Python VScode

per ignoreare errori quando si usa come linketer `flake8`.
Add to the file settings.json the following string

```json
"python.linting.flake8Args": [
    "--max-line-length=120",
    "--ignore=E402,F841,F401,E302,E305",
],
```

## Note su Miniconda
Install python
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    chmod u+x Miniconda3-latest-Linux-x86_64.sh
    ./Miniconda3-latest-Linux-x86_64.sh 
    // scegliere "lanciare conda init?" > YES!"

    conda create -n "conda-math" python=3.8.2

    conda activate conda-math
    conda install numpy scipy
    conda install matplotlib