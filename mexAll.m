% minFunc
fprintf('Compiling minFunc files...\n');
mex -outdir minFunc minFunc/mcholC.c
mex -outdir minFunc minFunc/lbfgsC.c

