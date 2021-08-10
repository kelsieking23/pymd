#!/bin/bash

source /usr/local/gromacs/bin/GMXRC
gmx dump -f $1 > $2