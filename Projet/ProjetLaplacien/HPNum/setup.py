#!/usr/bin/env python

from distutils.core import setup, Extension
import numpy as np

modules = [Extension('fem', ['CPyFem.cpp'],include_dirs=[np.__path__[0]+"/core/include"]),
           Extension('laplacian', ['CPyLaplacian.cpp'],include_dirs=[np.__path__[0]+"/core/include"]),
           Extension('mesh', ['CPyMesh.cpp'],include_dirs=[np.__path__[0]+"/core/include"]),
           Extension('splitter', sources=['splitter.cpp'], libraries=[],include_dirs=[np.__path__[0]+"/core/include"], extra_compile_args=['-std=c++11'])]

setup(name="HPNum",
      version="0.1", # Chaine de caractere pour ne pas avoir le probleme a la fin pour l'install du egginfo
      description='Python high performance numerical algorithms',
      author='Xavier Juvigny',
      author_email="juvigny@onera.fr",
      ext_modules = modules)

