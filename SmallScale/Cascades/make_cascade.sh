#!/bin/bash


energy=1e18
theta=45
Ndivs=10
prefix="E${energy}_T${theta}_N${Ndivs}" #just an example naming coefficient for the cascade


echo "python make_cascade.py -N $Ndivs -E $energy --theta $theta -O $prefix"
python3 make_cascade.py -N $Ndivs -E $energy -T $theta -O $prefix


