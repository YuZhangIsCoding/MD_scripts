#!/bin/bash

gen_itp_graphene_freestanding.py -sc $1 -c zigzag -o GPO.itp -n GPO -s 28 30
gen_itp_graphene_freestanding.py -sc -$1 -c zigzag -o GNE.itp -n GNE -s 28 30
