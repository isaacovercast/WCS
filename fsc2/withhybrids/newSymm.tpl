//Number of population samples (demes)
2
//Population effective sizes (number of genes)
NCL
NLN
//Sample sizes
75 0
75 0
//Growth rates	: negative growth implies population expansion
EXB
EXN
//Number of migration matrices : 0 implies no migration between demes
2
//Migration matrix 0
0.0 M12
M12 0.0
//Migration matrix 1
0.0 0.0
0.0 0.0
//historical event: time, source, sink, migrants, new size, new growth rate, migr. matrix 
3  historical event
TD1 0 0 1 1  0   1
TD1 1 1 1 1  0   1
TD2 1 0 1 RS 0   1
//Number of independent loci [chromosome]
1 0
//Per chromosome: Number of linkage blocks
1
//per Block: data type, num loci, rec. rate and mut rate + optional parameters
FREQ 1 0 2.5e-8
