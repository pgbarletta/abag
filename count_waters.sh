#!/bin/bash

for cada in `cat full_pdb.list`
do
    grep -c HOH structures/raw/$cada.pdb >> wat_count
done
