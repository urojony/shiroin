# Shiro Inequality Prover (ShiroIn)

A Python library for proving polynomial inequalities. It uses SciPy (more specifically: scipy.optimize.linprog) for computation and SymPy for reading and writing formulas.



## Prerequisites

* Python >=3.5

* pip

* git (if you install from GitHub)





## Installation

Copy one of these commands in terminal/cmd.

* Install from PyPI repository.

```pip install shiroin```

* Install from GitHub.

```pip install git+https://github.com/urojony/shiroin```



If you can't install libraries (for example you don't have administrator rights), try to add `--user` to command.

Note that ShiroIn works only for python3, if your default python is python2, then you should probably write `pip3` instead of `pip`.



## ShiroIn tutorial

For start, type ```python3``` in terminal/cmd.



```python

>>> from shiroin import *

>>> shiroSeed=1

```

First line obviously loads this package. Second one sets a seed for proving functions. If you don't write it, you can get slightly different proof each time you run the function. 



Now make some proofs. We will use problems from https://www.imomath.com/index.php?options=593&lmm=0.

The first inequality to prove is `a^2+b^2+c^2>= ab+ac+bc` for all real numbers. Function `prove` tries to prove that given formula is nonnegative, **assuming all variables are nonnegative**. In this case it's not a problem, since all powers on the left side are even, so if `|a|^2+|b|^2+|c|^2>= |ab|+|ac|+|bc|`, then `a^2+b^2+c^2= |a|^2+|b|^2+|c|^2>= |ab|+|ac|+|bc|>= ab+ac+bc`. There is a method we can use if this trick doesn't work, but it complicates the formula, so there is no need to use it here.

```python

>>> prove('a^2+b^2+c^2- (a*b+a*c+b*c)')

```

```

numerator: $$a^2-ab-ac+b^2-bc+c^2$$

denominator: $$1$$

status: 0

From weighted AM-GM inequality:



Program couldn't find a solution with integer coefficients. Try to multiple the formula by some integer and run this function again.

$$  ab+ac+bc  \le 

a^2+b^2+c^2  $$

```

Function prove prints several things. First two gives us a formula after expanding it. To proceed, a **numerator** has to be a **polynomial with integer coefficients**. Next one is status, which is return status of first use of ```scipy.optimize.linprog```. Possible outputs and explanations are

* 0 - found a proof with real coefficients,

* 1 - need more time, 

* 2 - function didn't find a proof,

* 3 - probably float loss of precision (which may happen if it has to work with numbers in the order of 2^46 or bigger).



Then we've got a hint. So let's use it!

```python

>>> prove('(a^2+b^2+c^2- (a*b+a*c+b*c))*2')

```

```

numerator: $$2a^2-2ab-2ac+2b^2-2bc+2c^2$$

denominator: $$1$$

status: 0

From weighted AM-GM inequality:

$$2ab \le a^2+b^2$$

$$2ac \le a^2+c^2$$

$$2bc \le b^2+c^2$$



$$  0  \le 

0  $$

Sum of all inequalities gives us a desired proof.

```

Now we have a proof! If you don't know weighted AM-GM inequality, check [wikipedia](https://en.wikipedia.org/wiki/Inequality_of_arithmetic_and_geometric_means#Weighted_AM%E2%80%93GM_inequality).

 

Now let's move to the next problem. Now we have to find all real numbers such that `a^2+b^2+c^2+d^2=a(b+c+d)`. At first glance it doesn't look like an inequality problem, but actually it is one. If you try to calculate both sides for different values, you can see that left side of the equation is never less than the right one. So let's try

```python

>>> prove('a^2+b^2+c^2+d^2-a*(b+c+d)')

```

```

numerator: $$a^2-ab-ac-ad+b^2+c^2+d^2$$

denominator: $$1$$

status: 2



Program couldn't find any proof.

```

This time `prove` didn't found the proof. But it doesn't mean that the inequality is not true! `prove` uses a list of values for which the formula should be small. There is no strict rule here, but the smaller that value is, the higher are chances to find a proof. List of values should correspond to the list of variables in alphabetical order. So let's try a=2 and b=c=d=1.



```python

>>> prove('a^2+b^2+c^2+d^2-a*(b+c+d)','2,1,1,1')

```

```

Substitute $ a \to 2a $

numerator: $$4a^2-2ab-2ac-2ad+b^2+c^2+d^2$$

denominator: $$1$$

status: 0

From weighted AM-GM inequality:

$$2ab \le a^2+b^2$$

$$2ac \le a^2+c^2$$

$$2ad \le a^2+d^2$$



$$  0  \le 

a^2  $$

Sum of all inequalities gives us a desired proof.

```

Function makes a substitution `a->2a` (which should be understood as `a=2a'`, but the author of the code is too lazy to rename variables) and try to prove new inequality. This time it succeeded. Moreover, if starting formula is equal to 0, then all these inequalities have to be equalities, so `a'^2=0` and eventually `a=0`. We can also try a little bit lower value for `a`.



```python

>>> prove('a^2+b^2+c^2+d^2-a*(b+c+d)','7/4,1,1,1')

```

```

Substitute $ a \to 7a/4 $

numerator: $$49a^2-28ab-28ac-28ad+16b^2+16c^2+16d^2$$

denominator: $$16$$

status: 0

From weighted AM-GM inequality:

$$28ab \le 14a^2+14b^2$$

$$28ac \le 14a^2+14c^2$$

$$28ad \le 14a^2+14d^2$$



$$  0  \le 

7a^2+2b^2+2c^2+2d^2  $$

Sum of all inequalities gives us a desired proof.

```

Now we can see that if `a^2+b^2+c^2+d^2-a(b+c+d)=0`, then `7a'^2+2b^2+2c^2+2d^2=0` and eventually `a=b=c=d=0`. Note that it looks only for nonnegative solutions. But using similar argumentation to the one in previous problem, if `(a,b,c,d)=(x,y,z,t)` is the solution of `a^2+b^2+c^2+d^2-a(b+c+d)=0`, then `(a,b,c,d)=(|x|,|y|,|z|,|t|)` is a solution, too. Since the only nonnegative solution is (0,0,0,0), it means that it is the only solution.



Let's skip the problem 3 and look solve the problem 4 instead. Now the inequality to prove is

`1/(1-x^2)+1/(1-y^2)>=2/(1-xy)` for `0



```python

>>> prove('1/(1-x^2)+1/(1-y^2)-2/(1-x*y)')

```

```

numerator: $$-x^3y+2x^2y^2-x^2-xy^3+2xy-y^2$$

denominator: $$x^3y^3-x^3y-x^2y^2+x^2-xy^3+xy+y^2-1$$

status: 2



Program couldn't find any proof.

```

`prove` assumes that formula is well-defined if all variables are positive, so it doesn't have to analyze the denominator (except of choosing the right sign). In this case it is not true, since if `x=1`, then `1-x^2=0`. Denominator is equal to `(x^2-1)(y^2-1)(xy-1)` which is negative for `0

So we will use a function makesubs to generate this substitutions. It has three basic parameters: `formula`, `intervals` and `values`. `intervals` are current limitations of variables, `values` are values of variables for which `formula` is small. Values should be inside intervals. This argument is optional but it's better to use it.

Let's go back to our problem. If `x=y`, then `1/(1-x^2)+1/(1-y^2)-2/(1-x*y)=0`, so it's the minimum value. So let values be `1/2,1/2` (**warning: do not use decimal point**, for example '0.5,0.5').

```python

>>> newformula,values=makesubs('1/(1-x^2)+1/(1-y^2)-2/(1-x*y)','[0,1],[0,1]','1/2,1/2')

```

```

Substitute $ x \to 1-1/(x+1) $

Substitute $ y \to 1-1/(y+1) $

```

```python

>>> prove(newformula*6,values)

```

```

numerator: $$12x^3y+6x^3-24x^2y^2-6x^2y+6x^2+12xy^3-6xy^2-12xy+6y^3+6y^2$$

denominator: $$4x^2y+2x^2+4xy^2+8xy+3x+2y^2+3y+1$$

status: 0

From weighted AM-GM inequality:

$$24x^2y^2 \le 12x^3y+12xy^3$$

$$6x^2y \le 4x^3+2y^3$$

$$6xy^2 \le 2x^3+4y^3$$

$$12xy \le 6x^2+6y^2$$



$$  0  \le 

0  $$

Sum of all inequalities gives us a desired proof.

```



Now let's get back to problem 3. Now our goal is to find the minimum of `a^2*b^2/c^2+b^2*c^2/a^2+c^2*a^2/b^2` assuming `a^2+b^2+c^2=1` and `a,b,c>0`. It is equivalent problem to finding minimum of `xy/z+yz/x+zx/y` assuming `x+y+z=1` and `x,y,z>0`. The first idea is to suppose that the minimum is reached when `x=y=z`. In that case, `x=y=z=1/3` and formula is equal to 1. Now we can substitute `z-> 1-x-y`. Constraints for variables are `x>0`, `y>0`, `x+y<1`. We can rewrite it as `x in (0,1-y)`, `y in (0,1)`. It's very similar to defining integration limits. After you define interval of x, you cannot use x in any next interval. **Warning:** at this moment `makesubs` **doesn't warn you if your intervals doesn't follow this rule!**

```python

>>> formula=Sm('xy/z+yz/x+zx/y-1').subs('z',S('1-x-y'))

>>> newformula,values=makesubs(formula,'[0,1-y],[0,1]','1/3,1/3')

```

```

Substitute $ x \to -y+1+(y-1)/(x+1) $

Substitute $ y \to 1-1/(y+1) $

```

```python

prove(newformula,values)

```

```

Substitute $ y \to y/2 $

numerator: $$x^4y^2+x^3y^2-2x^3y-4x^2y+4x^2+xy^2-2xy+y^2$$

denominator: $$x^3y^2+2x^3y+2x^2y^2+4x^2y+xy^2+2xy$$

status: 0

From weighted AM-GM inequality:

$$2x^3y \le x^4y^2+x^2$$

$$4x^2y \le x^3y^2+2x^2+xy^2$$

$$2xy \le x^2+y^2$$



$$  0  \le 

0  $$

Sum of all inequalities gives us a desired proof.

```

Proof is found, so the assumption that 1 is the minimum of `xy/z+yz/x+zx/y` was good.

Functions `S` and `Sm` creates a SymPy object from a string. The only difference is that `Sm` assumes that there are no multi-letter variables and adds a multiplication sign between every two terms which has no operator sign, so object `Sm(xy/z+yz/x+zx/y)` has 3 variables `x,y,z` and `S('xy/z+yz/x+zx/y')` has 6 variables `x,y,z,xy,yz,zx`.




