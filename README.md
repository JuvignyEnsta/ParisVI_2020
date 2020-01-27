# ParisVI_2020
TP calcul intensif parallele Paris 6, Premier trimestre civil 2020

Cours, TP et projets pour la promotion 2021

## Installation des outils nécessaires aux TDs

### Linux/debian 

    sudo apt install  build-essential make g++ gdb libopenmpi-dev python3-numpy python3-mpi4py

### Windows 10

Installation d'ubuntu sous Windows 10 :

Allez sur le lien suivant  https://www.numerama.com/tech/158150-le-shell-bash-sous-windows-10-ce-quil-faut-savoir.html et suivez les instructions

Taper bash dans la bare de questionnement.
Une fois sous Linux :

    sudo apt update
    sudo apt install  build-essential make g++ gdb libopenmpi-dev python3-numpy python3-mpi4py

### Windows 7 ou 8

Installer msys 2 via le site :

https://www.msys2.org/


Lancer le shell via mingw64

    pacman -Ss gcc  # donne le nom du package contenant gcc
    pacman -S mingw64/mingw-w64-x86_64-gcc   # install le package gcc avec le nom copié/collé précédent
    pacman -Syu # Màj du système

Installer gcc, make

Pour MPI, installer MS-MPI : https://www.microsoft.com/en-us/download/details.aspx?id=57467

Dans .bashrc :

    export PATH=$PATH:"/c/Program Files/Microsoft MPI/Bin/"


Exemple de commande de compilation :

    g++ -o HelloWorld.exe HelloWorld.cpp -L/c/Program\ Files\ \(x86\)/Microsoft\ SDKs/MPI/Lib/x64 -I/c/Program\ Files\ \(x86\)/Microsoft\ SDKs/MPI/Include -lmsmpi

Pour mpi4py, télécharger le fichier d'installation. Demander à Xavier Juvigny pour changer le setup.py pour l'adapter à windows...

### Mac

Installer homebrew => https://brew.sh

Suivre les instructions

    brew install gcc open-mpi numpy scipy matplotlib
    brew update # màj

Puis 
    pip3 install mpi4py

A mettre dans le .bashrc:

    export OMPI_CC=gcc-9
    export OMPI_CXX=g++-9

