.. py:function:: prove(formula[,values=None,variables=None,niter=200,linprogiter=10000,translation={}]):

Try to prove that given formula is nonnegative assuming all variables are positive and formula is well defined for all positive values of variables.

Parameters:
formula - a rational function with integer coefficients.
values - vector of n positive rational values where n is a number of variables. If not specified, values=[1]*n.
variables - variables corresponding to values. If not specified, variables are sorted list of free symbols in the formula.
niter - maximum number of attempts to find a proof.
linprogiter - maximum number of iterations of internal function.
translation - a dictionary of phrases used in the function. See shiroTranslation for more details.

prove doesn't return any value - it only prints a proof if it finds one. It also shows status of the first attempt:

0 - found a proof with real coefficients,
1 - need bigger linprogiter parameter,
2 - function didn't find a proof,
3 - this shouldn't happen, but if it does, it probably means the same as 4,
4 - coefficients and/or powers were too big for a program.


Hint:
try to give such values, for which the formula is a small number. In particular, if there exist points for which the formula is equal to 0, one of those points should be given as this parameter.
