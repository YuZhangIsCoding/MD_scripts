#!/bin/bash

assign_itp_mxn_o.py -i base.gro -sc $1 -n MPO -o MPO.itp
assign_itp_mxn_o.py -i base.gro -sc -$1 -n MNE -o MNE.itp
