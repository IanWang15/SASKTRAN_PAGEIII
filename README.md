# SASKTRAN_SAGEIII
SAGEIII project by SASKTRAN simulation

# Installation

## Build a python virtual environment on AWS

$conda create -n py36rtm python=3.6

(base) $ conda activate py36rtm

(py36rtm) $ conda install -c anaconda numpy

(py36rtm) $ conda install -c anaconda mpi4py

(py36rtm) $ echo ${PYTHONPATH}

(py36rtm) $ pip install sasktran -f https://arg.usask.ca/wheels/![image](https://user-images.githubusercontent.com/52504365/114101263-167b7d80-988b-11eb-83e2-23e7949413b7.png)

## update SASKTRAN to the develop branch:

$pip uninstall sasktran

$pip install sasktran -f https://arg.usask.ca/wheels/dev

the documentation here https://arg.usask.ca/docs/sasktran/dev.![image](https://user-images.githubusercontent.com/52504365/114101471-635f5400-988b-11eb-9aff-1bf97c0fff73.png)



