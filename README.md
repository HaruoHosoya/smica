This is code implementing an algorithm for overcomplete ICA based on the score matching method proposed by Aapo Hyvarinen[1].

The code relies on Schmidt's implementation of unconstrained opitmization, included in this distribution [2].

To use this code, first run initpath.m.  The entry function is smica.m.

You might need to recompile the optimizer, in which case please consult [2].  (You likely need only to run mexAll.m.)

[1] A. Hyvarinen. Estimation of non-normalized statistical models using score matching. Journal of Machine Learning Research, 6:695--709, 2005. 

[2] https://www.cs.ubc.ca/~schmidtm/Software/minFunc.html
