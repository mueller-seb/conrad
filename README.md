conrad [![Build Status](https://travis-ci.org/bungun/conrad.svg?branch=master)](https://travis-ci.org/bungun/conrad) [![Documentation Status](https://readthedocs.org/projects/conrad/badge/?version=latest)](http://conrad.readthedocs.io/en/latest/?badge=latest)
===

**con**vex **rad**iation treatment planning

# HowTo
0. (install gcc)
1. create python 3.7 environment, e.g. by (ana)conda:

        conda create -n conrad python=3.7

2.      git clone https://github.com/mueller-seb/conrad/
3. select the correct branch

        git checkout cvxpy_1.0

4. downgrade setuptools<58 to avoid use_2to3 error

        pip install setuptools==57.5.0

5.      pip install -r requirements.txt
6.      pip install -e .

In case of problems compare library versions with requirements_py37.txt

## Build Documentation
7.      cd /docs
8.      pip install -r requirements.txt
9.      pip install mock
10.     install make
11.     make html


## DEPRECATED BRANCH (DO NOT USE)
Don't use the master branch of conrad, it is deprecated as it depends on CVXPY 0.x
However if you do so, you have to use gcc<=10 to compile old CVXPY libraries.

- multiprocess==0.70.5 if error "Python before 3.6 not supported" occurs
- CVXcanon==0.1.1 if error "canonInterface.py" occurs in hello_world.py
- scs==1.2.7 if SCS solving error occurs in CVX
- install matplotlib if 'NoneType' object has no attribute dvh_plot

see /oldrequirements for a full list of libraries
