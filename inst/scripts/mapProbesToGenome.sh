#!/bin/sh
cd SGD-0508
## cat chr01.fsa chr02.fsa chr03.fsa chr04.fsa chr05.fsa chr06.fsa chr07.fsa chr08.fsa chr09.fsa chr10.fsa chr11.fsa chr12.fsa chr13.fsa chr14.fsa chr15.fsa chr16.fsa  chrmt.fsa > scAll.fsa
/ebi/research/huber/users/huber/MUMmer3.18/mummer  -maxmatch -F -L -n -b -c -l 23  scAll.fsa ../Scerevisiaetilingprobe.fsa  > scAll.out
cd -
