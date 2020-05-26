from __future__ import (absolute_import, division,
						print_function, unicode_literals)
import warnings,operator
warnings.filterwarnings("ignore")	
def vargen(n): 
	"""default function generating names for variables"""
	q=len(shiro.alphabet)
	x=shiro.alphabet[n%q]
	if n>=q:
		x+=str(n//q)
	return x
class Vars:
	pass
shiro=Vars()
shiro.display=shiro.warning=print
shiro.seed=None 
"""Seed is needed to select the weights in linprog function.
None means that the seed is random"""
shiro.translation={}
shiro.varind=0
shiro.varset=set() 
"""set of used symbols in the proof"""
shiro.vargen=vargen 
"""function generating names for variables"""
shiro.alphabet='abcdefghijklmnopqrstuvwxyz' 
"""list of names for variables"""
translationList=['numerator:','denominator:','status:',
'Substitute',"Formula after substitution:",
"Numerator after substitutions:","From weighted AM-GM inequality:",
'The sum of all inequalities gives us a proof of the inequality.',
"Program couldn't find a solution with integer coefficients. Try "+
"to multiple the formula by some integer and run this function again.",
"Program couldn't find any proof.",
"Try to set higher linprogiter parameter.",
"It looks like the formula is symmetric. You can assume without loss of"+
" generality that ","Try", 'From Jensen inequality:',
'Warning: intervals contain backwards dependencies. Consider changing order of variables and intervals.'
]
#Initialize english-english dictionary.
for phrase in translationList:
	shiro.translation[phrase]=phrase
from scipy.optimize import linprog,fmin
import random
from sympy import S,cancel,fraction,Pow,expand,solve,latex,oo,Poly,lambdify,srepr,gcd,Symbol
from sympy.parsing.sympy_parser import parse_expr, standard_transformations,\
implicit_multiplication_application, convert_xor
from collections import Counter
import re
def addsymbols(formula):
	formula=S(formula)
	funcsymbols=[x[10:-2] for x in re.findall(r"Function\(\'.*?\'\)",srepr(formula))]
	shiro.varset|=set(funcsymbols)|set(map(str,formula.free_symbols))
def newvar():
	while 1:
		x=shiro.vargen(shiro.varind)
		shiro.varind+=1
		if x not in shiro.varset:
			return S(x)
def newproof():
	shiro.varset=set()
	shiro.varind=0
def _remzero(coef,fun): 
	"""coef, fun represents an expression.
	For example, if expression=5f(2,3)+0f(4,6)+8f(1,4)
	then coef=[5,0,8], fun=[[2,3],[4,6],[1,4]]
	_remzero removes addends with coefficient equal to zero.
	In this example ncoef=[5,8], nfun=[[2,3],[1,4]]"""
	ncoef=[]
	nfun=[]
	for c,f in zip(coef,fun):
		if c>0:
			ncoef+=[c]
			nfun+=[f]
	return ncoef,nfun
# ~ def slatex(formula): #fancy function which makes latex code more readable, but still correct
	# ~ formula=re.sub(r'\^{(.)}',r'^\1',latex(formula,fold_short_frac=True).replace(' ','').replace('\\left(','(').replace('\\right)',')'))
	# ~ return re.sub(r'\{(\(.+?\))\}',r'\1',formula)
def _writ2(coef,fun,variables):
	return latex(Poly({fun:coef},gens=variables).as_expr())
def _writ(coef,fun,nullvar):
	return str(coef)+'f('+str(fun)[1:-1-(len(fun)==1)]+')'
def _check(coef,fun,res,rfun): 
	"""checks if rounding and all the floating point stuff works"""
	res2=[int(round(x)) for x in res]
	b1=[coef*x for x in fun]
	b2=[[x*y for y in rfuni] for x,rfuni in zip(res2,rfun)]
	return b1==[sum(x) for x in zip(*b2)] and coef==sum(res2)
def _powr(formula):
	if formula.func==Pow:
		return formula.args
	else:
		return [formula,S('1')]
def fractioncancel(formula): 
	"""workaround for buggy cancel function"""
	num,den=fraction(cancel(formula/S('tmp')))
	den=den.subs(S('tmp'),S('1'))
	return num,den
def ssolve(formula,variables): 
	"""workaround for inconsistent solve function"""
	result=solve(formula,variables)
	if type(result)==dict:
		result=[[result[var] for var in variables]]
	return result
#def sstr(formula):
#	return str(formula).replace('**','^').replace('*','').replace(' ','')
def assumeall(formula,**kwargs):
	"""Adds assumptions to all free symbols in formula.
	>>> assumeall('sqrt(x*y)-sqrt(x)*sqrt(y)',positive=True)
	0
	"""
	formula=S(formula)
	fs=formula.free_symbols
	for x in fs:
		y=Symbol(str(x),**kwargs)
		formula=formula.subs(x,y)
	return formula
def reducegens(formula):
	"""Reduces size of the generator of the polynomial
	>>> Poly('x+sqrt(x)')
	Poly(x + (sqrt(x)), x, sqrt(x), domain='ZZ')
	>>> reducegens('x+sqrt(x)')
	Poly((sqrt(x))**2 + (sqrt(x)), sqrt(x), domain='ZZ') """
	pol=Poly(formula)
	newgens={}
	ind={}
	for gen in pol.gens:
		base,pw=_powr(gen)
		coef,_=pw.as_coeff_mul()
		ml=pw/coef
		if base**ml in newgens:
			newgens[base**ml]=gcd(newgens[base**ml],coef)
		else:
			newgens[base**ml]=coef
			ind[base**ml]=S('tmp'+str(len(ind)))
	for gen in pol.gens:
		base,pw=_powr(gen)
		coef,_=pw.as_coeff_mul()
		ml=pw/coef
		pol=pol.replace(gen,ind[base**ml]**(coef/newgens[base**ml]))
	newpol=Poly(pol.as_expr())
	for gen in newgens:
		newpol=newpol.replace(ind[gen],gen**newgens[gen])
	return newpol
def Sm(formula):
	"""Adds multiplication signs and sympifies a formula.
	>>> Sm('(2x+y)(7+5xz)') 
	(2*x + y)*(5*x*z + 7)"""
	if type(formula)!=str:
		if _isiterable(formula): return type(formula)(map(Sm,formula))
		else: return S(formula)
	formula=formula.replace('{',' ').replace('}',' ')
	transformations = (standard_transformations +(implicit_multiplication_application,convert_xor))
	return parse_expr(formula, transformations=transformations)
# ~ def Sm(formula):
# ~ #Adds multiplication signs and sympifies a formula.
# ~ #For example, Sm('(2x+y)(7+5xz)') -> S('(2*x+y)*(7+5*x*z)')
	# ~ if type(formula)==str:
		# ~ formula=formula.replace(' ','')
		# ~ for i in range(2):
			# ~ formula=re.sub(r'([0-9a-zA-Z)])([(a-zA-Z])',r'\1*\2',formula)
	# ~ formula=S(formula)
	# ~ return formula
def _input2fraction(formula,variables,values):
	"""makes some substitutions and converts formula to a fraction
	with expanded numerator and denominator"""
	formula=S(formula)
	subst=[]
	for x,y in zip(variables,values):
		if y!=1:
			z=newvar()
			shiro.display(shiro.translation['Substitute']+' $'+latex(x)+'\\to '+latex(S(y)*S(z))+'$')
			subst+=[(x,z*y)]
	formula=formula.subs(subst)
	numerator,denominator=fractioncancel(formula)
	shiro.display(shiro.translation['numerator:']+' $'+latex(numerator)+'$')
	shiro.display(shiro.translation['denominator:']+' $'+latex(denominator)+'$')
	return (numerator,denominator)
def _formula2list(formula):
	"""Splits a polynomial to a difference of two polynomials with positive
	coefficients and extracts coefficients and powers of both polynomials.
	'variables' is used to set order of powers
	For example, If formula=5x^2-4xy+8y^3, variables=[x,y], then
	the program tries to prove that 
	0<=5x^2-4xy+8y^3
	4xy<=5x^2+8y^3
	returns [4],[(1,1)], [5,8],[(2,0),(0,3)], (x,y)"""
	formula=reducegens(assumeall(formula,positive=True))
	neg=(formula.abs()-formula)*S('1/2')
	pos=(formula.abs()+formula)*S('1/2')
	return neg.coeffs(),neg.monoms(),pos.coeffs(),pos.monoms(),Poly(formula).gens
def _list2proof(lcoef,lfun,rcoef,rfun,variables,itermax,linprogiter,_writ2=_writ2,theorem="From weighted AM-GM inequality:"):
	"""Now the formula is splitted on two polynomials with positive coefficients.
	we will call them LHS and RHS and our inequality to prove would
	be LHS<=RHS (instead of 0<=RHS-LHS).
		
	Suppose we are trying to prove that 
	30x^2y^2+60xy^4<=48x^3+56y^6 (assuming x,y>0).
	Program will try to find some a,b,c,d such that
	30x^2y^2<=ax^3+by^6  
	60xy^4<=cx^3+dy^6  
	where a+c<=48 and b+d<=56 (assumption 1).
	We need some additional equalities to meet assumptions
	of the weighted AM-GM inequality.
	a+b=30 and c+d=60 (assumption 2)
	3a+0b=30*2, 0a+6b=30*2, 3c+0d=60*1, 0c+6d=60*4 (assumption 3)

	The sketch of the algorithm.
	 for i in range(itermax):
		1. Create a vector of random numbers (weights).
		2. Try to find real solution of the problem (with linprog).
		3. If there is no solution (status: 2)
			3a. If the solution was never found, break.
			3b. Else, step back (to the bigger inequality)
		4. If the soltuion was found (status: 0)
			Check out which of variables (in example: a,b,c,d) looks like integer.
			If there are some inequalities with all integer coefficients, subtract 
			them from the original one.
			If LHS is empty, then break.""" 
	localseed=shiro.seed
	bufer=[]
	lcoef,lfun=_remzero(lcoef,lfun)
	rcoef,rfun=_remzero(rcoef,rfun)
	itern=0
	if len(lcoef)==0: #if LHS is empty
		shiro.display(shiro.translation['status:']+' 0')
		status=0
	elif len(rcoef)==0: 
	#if RHS is empty, but LHS is not
		shiro.display(shiro.translation['status:']+' 2')
		status=2
		itermax=0
	foundreal=0
	while len(lcoef)>0 and itern<itermax:
		itern+=1
		m=len(lcoef)
		n=len(rcoef)
		#lfunt=transposed matrix lfun (in fact, it's 
		#a list of lists)
		lfunt=list(map(list, zip(*lfun)))
		rfunt=list(map(list, zip(*rfun)))
		#A,b, - set of linear equalities 
		#A_ub,b_ub - set of linear inequalities
		#from linear program
		A=[]
		b=[]
		A_ub=[]
		b_ub=[]
		for i in range(m):
			for j in range(len(rfunt)): #assumption 3
				A+=[[0]*(m*n)] 
				A[-1][i*n:i*n+n]=rfunt[j]
				b+=[lfun[i][j]*lcoef[i]]
			A+=[[0]*(m*n)] #assumption 2
			A[-1][i*n:i*n+n]=[1]*n
			b+=[lcoef[i]]
		for j in range(n): #assumption 1
			A_ub+=[[0]*(m*n)] 
			A_ub[-1][j::n]=[1]*m
			b_ub+=[rcoef[j]]
		random.seed(localseed)
		vecc=[random.random() for i in range(m*n)]
		localseed=random.randint(1,1000000000)
		res=linprog(vecc,A_eq=A,b_eq=b,A_ub=A_ub,b_ub=b_ub,options={'maxiter':linprogiter})
		status=res.status
		if itern==1:
			shiro.display(shiro.translation['status:']+' '+str(status))
			if status==0:
				shiro.display(shiro.translation[theorem])
		if status==2: #if real solution of current inequality doesn't exist
			if foundreal==0: #if this is the first inequality, then break 
				break
			else:
				#step back 
				lcoef,lfun=oldlcoef,oldlfun
				rcoef,rfun=oldrcoef,oldrfun
				bufer=[]
				continue
		if status==0:#if found a solution with real coefficients
			for ineq in bufer:
				shiro.display(ineq)
			foundreal=1
			bufer=[]
			oldlfun,oldrfun=lfun,rfun
			oldlcoef,oldrcoef=lcoef[:],rcoef[:]
			for i in range(m):
				c=0
				for j in res.x[i*n:i*n+n]:#check if all coefficients
					#in an equality looks like integers
					if(abs(round(j)-j)>0.0001): 
						break
				else:
					#checks if rounding all coefficients doesn't make
					#inequality false
					isok=_check(lcoef[i],lfun[i],res.x[i*n:i*n+n],rfun)
					if not isok:
						continue
					bufer+=['']
					bufer[-1]+='$$'+_writ2(lcoef[i],lfun[i],variables)+' \\le '
					lcoef[i]=0
					for j in range(n):
						rcoef[j]-=int(round(res.x[i*n+j]))
					for j,k in zip(res.x[i*n:i*n+n],rfun):
						if j<0.0001:
							continue
						if(c):bufer[-1]+='+'
						else:c=1
						bufer[-1]+=_writ2(int(round(j)),k,variables)
					bufer[-1]+='$$'
			lcoef,lfun=_remzero(lcoef,lfun)
			rcoef,rfun=_remzero(rcoef,rfun)
	for ineq in bufer:
		shiro.display(ineq)
	lhs='+'.join([_writ2(c,f,variables) for c,f in zip(lcoef,lfun)])
	if lhs=='':
		lhs='0'
	elif status==0:
		shiro.display(shiro.translation[
		"Program couldn't find a solution with integer coefficients. Try "+
		"to multiple the formula by some integer and run this function again."])
		status=-1
	elif(status==2):
		shiro.display(shiro.translation["Program couldn't find any proof."])
		#return res.status
	elif status==1:
		shiro.display(shiro.translation["Try to set higher linprogiter parameter."])
	rhs='+'.join([_writ2(c,f,variables) for c,f in zip(rcoef,rfun)])
	if rhs=='':
		rhs='0'
	shiro.display('$$ '+latex(lhs)+' \\le '+latex(rhs)+' $$')
	if lhs=='0':
		shiro.display(shiro.translation['The sum of all inequalities gives us a proof of the inequality.'])
	return status
def _isiterable(obj):
	try:
		_ = (e for e in obj)
		return True
	except TypeError:
		return False
def _smakeiterable(x):
	if x=='':
		return []
	x=S(x)
	if _isiterable(x):
		return x
	return (x,)
def _smakeiterable2(x):
	if x=='':
		return []
	x=S(x)
	if len(x)==0:
		return []
	if _isiterable(x[0]):
		return x
	return (x,)
def prove(formula,values=None,variables=None,niter=200,linprogiter=10000):
	"""Tries to prove that formula>=0 assuming all variables are positive.
	
	Arguments:
	formula - sympy object or a string;
	values - list of numbers,
	by default values=[1]*len(variables);
	variables - list of variables corresponding to values,
	by default variables=sorted(formula.free_symbols,key=str);
	niter - number of iterations of the main loop;
	linprogiter - number of iterations in the 
	underlying scipy.optimize.linprog function.
	
	Return codes:
	-1 - found a proof with real coefficients, but not with the integer ones,
	0 - found a proof with integer coefficients,
	1 - linprogiter parameter was too small,
	2 - no proof was found,
	3/4 - the coefficients or powers of the formula were too big (may depend
	on version of SciPy).
	
	>>> prove('x^2-2*x+1')
	0
	>>> prove('x^2-4*x+4')
	2
	>>> prove('x^2-4*x+4',values='2')
	0
	
	In addition, prove displays a proof (if it has found one). The form of
	display depend on shiro.display.
	
	"""
	formula=S(formula)
	addsymbols(formula)
	if variables: variables=_smakeiterable(variables)
	else: variables=sorted(formula.free_symbols,key=str)
	if values: values=_smakeiterable(values)
	else: values=[1]*len(variables)
	num,den=_input2fraction(formula,variables,values)
	st=_list2proof(*(_formula2list(num)+(niter,linprogiter)))
	if st==2 and issymetric(num):
		fs=sorted([str(x) for x in num.free_symbols])
		shiro.display(shiro.translation["It looks like the formula is symmetric. "+
		"You can assume without loss of generality that "]+
		' >= '.join([str(x) for x in fs])+'. '+shiro.translation['Try'])
		shiro.display('prove(makesubs(S("'+str(num)+'"),'+
		str([(str(x),'oo') for x in variables[1:]])+')')
	return st
def powerprove(formula,values=None,variables=None,niter=200,linprogiter=10000):
	"""This is a bruteforce and ineffective function for proving inequalities.
	It can be used as the last resort.
	
	Arguments:
	formula - sympy object or a string;
	values - list of numbers,
	by default values=[1]*len(variables);
	variables - list of variables corresponding to values,
	by default variables=sorted(formula.free_symbols,key=str);
	niter - number of iterations of the main loop;
	linprogiter - number of iterations in the 
	underlying scipy.optimize.linprog function.
	
	Returns a counter of codes returned by 'prove' function.
	
	>>> powerprove('(x^2-1)^2')
	Counter({0: 2})"""
	formula=S(formula)
	addsymbols(formula)
	if variables: variables=_smakeiterable(variables)
	else: variables=sorted(formula.free_symbols,key=str)
	if values: values=_smakeiterable(values)
	else: values=[1]*len(variables)
	num,den=_input2fraction(formula,variables,values)
	statusses=[]

	for i in range(1<<len(variables)): #tricky substitutions to improve speed
		shiro.display(r'_______________________')
		subst1=[]
		subst2=[]
		substout=[]
		for j in range(len(variables)):
			z=newvar()
			x=variables[j]
			subst2+=[(x,1+z)]
			if i&(1<<j):
				subst1+=[(x,1/x)]
				substout+=[latex(x)+'\\to 1/(1+'+latex(z)+')']
			else:
				substout+=[latex(x)+'\\to 1+'+latex(z)]
		shiro.display(shiro.translation['Substitute']+ ' $'+','.join(substout)+'$')
		num1=fractioncancel(num.subs(subst1))[0]
		num2=expand(num1.subs(subst2))
		shiro.display(shiro.translation["Numerator after substitutions:"]+' $'+latex(num2)+'$')
		statusses+=[_list2proof(*(_formula2list(num2)+(niter,linprogiter)))]
	return Counter(statusses)
def makesubs(formula,intervals,values=None,variables=None,numden=False):
	"""Generates a new formula which satisfies this condition:
	for all positive variables new formula is nonnegative iff
	for all variables in corresponding intervals old formula is nonnegative.
	>>> newproof()
	>>> makesubs('1-x^2','[0,1]')
	Substitute $x\to 1 - \frac{1}{a + 1}$ #depend on shiro.display
	(2*a + 1)/(a**2 + 2*a + 1)
	>>> makesubs('1-x^2','[0,1]',values='1/2')
	Substitute $x\to 1 - \frac{1}{b + 1}$ #depend on shiro.display
	((2*b + 1)/(b**2 + 2*b + 1), [1])
	>>> makesubs('1-x^2','[0,1]',values='1/2',numden=True)
	Substitute $x\to 1 - \frac{1}{c + 1}$ #depend on shiro.display
	(2*c + 1, c**2 + 2*c + 1, [1])
	"""
	formula=S(formula)
	addsymbols(formula)
	intervals=_smakeiterable2(intervals)
	if variables: variables=_smakeiterable(variables)
	else: variables=sorted(formula.free_symbols,key=str)
	if values!=None:
		values=_smakeiterable(values)
		equations=[var-value for var,value in zip(variables,values)]
	else:
		equations=[]
	newvars=[]
	warn=0
	usedvars=set()
	for var,interval in zip(variables,intervals):
		end1,end2=interval
		z=newvar()
		newvars+=[z]
		usedvars|={var}
		if (end1.free_symbols|end2.free_symbols)&usedvars:
			warn=1
		if end1 in {S('-oo'),S('oo')}:
			end1,end2=end2,end1
		if {end1,end2}=={S('-oo'),S('oo')}:
			sub1=sub2=(z-1/z)
		elif end2==S('oo'):
			sub1=sub2=(end1+z)
		elif end2==S('-oo'):
			sub1=sub2=end1-z
		else:
			sub1=end2+(end1-end2)/z
			sub2=end2+(end1-end2)/(1+z)
		formula=formula.subs(var,sub1)
		shiro.display(shiro.translation['Substitute']+" $"+latex(var)+'\\to '+latex(sub2)+'$')
		equations=[equation.subs(var,sub1) for equation in equations]
	num,den=fractioncancel(formula)
	for var,interval in zip(newvars,intervals):
		if {interval[0],interval[1]} & {S('oo'),S('-oo')}==set():
			num=num.subs(var,var+1)
			den=den.subs(var,var+1)
			equations=[equation.subs(var,var+1) for equation in equations]
	if values:
		values=ssolve(equations,newvars)
		if len(values):
			values=values[0]
	num,den=expand(num),expand(den)
	#shiro.display(shiro.translation["Formula after substitution:"],"$$",latex(num/den),'$$')
	if warn:
		shiro.warning(shiro.translation[
		'Warning: intervals contain backwards dependencies. Consider changing order of variables and intervals.'])
	if values and numden:
		return num,den,values
	elif values:
		return num/den,values
	elif numden:
		return num,den
	else:
		return num/den
def _formula2listf(formula):
	"""Splits a polynomial to a difference of two formulas with positive
	coefficients and extracts coefficients and function
	arguments of both formulas."""
	lfun=[]
	lcoef=[]
	rfun=[]
	rcoef=[]
	for addend in formula.as_ordered_terms():
		coef,facts=addend.as_coeff_mul()
		if(coef<0):
			lcoef+=[-coef]
			lfun+=[facts[0].args]
		else:
			rcoef+=[coef]
			rfun+=[facts[0].args]
	return(lcoef,lfun,rcoef,rfun,None)
def provef(formula,niter=200,linprogiter=10000):
	"""This function is similar to prove, formula is a linear combination of
	values of f:R^k->R instead of a polynomial. provef checks if a formula
	is nonnegative for any nonnegative and convex function f. If so, it
	displays a proof of nonnegativity. 
	
	Return codes:
	-1 - found a proof with real coefficients, but not with the integer ones,
	0 - found a proof with integer coefficients,
	1 - linprogiter parameter was too small,
	2 - no proof was found,
	3/4 - the coefficients or arguments in the formula were too big (may depend
	on version of SciPy).
	
	>>> provef('f(5)-f(2)+f(1)')
	-1
	"""
	formula=S(formula)
	addsymbols(formula)
	num,den=_input2fraction(formula,[],[])
	return _list2proof(*(_formula2listf(num)+(niter,linprogiter,_writ,'From Jensen inequality:')))
def issymetric(formula): 
	"""checks if formula is symmetric
	and has at least two variables"""
	formula=S(formula)
	addsymbols(formula)
	if len(formula.free_symbols)<2:
		return False
	ls=list(formula.free_symbols)
	a=ls[0]
	for b in ls[1:]:
		if expand(formula-formula.subs({a:b, b:a}, simultaneous=True))!=S(0):
			return False
	return True
def cyclize(formula,oper=operator.add,variables=None,init=None):
	""">>>cyclize('a^2*b')
	a**2*b + a*b**2
	>>> cyclize('a^2*b',variables='a,b,c')
	a**2*b + a*c**2 + b**2*c
	"""
	formula=S(formula)
	addsymbols(formula)
	if variables==None:
		variables=sorted(formula.free_symbols,key=str)
	else:
		variables=S(variables)
	if len(variables)==0:
		return init
	variables=list(variables) #if variables is a tuple, change it to a list
	variables+=[variables[0]]
	subst=list(zip(variables[:-1],variables[1:]))
	if init==None:
		init=formula
	else:
		init=oper(init,formula)
	for _ in variables[2:]:
		formula=formula.subs(subst,simultaneous=True)
		init=oper(init,formula)
	return init
def symmetrize(formula,oper=operator.add,variables=None,init=None):
	""">>> symmetrize('a^2*b')
	a^2*b+b^2*a
	>>> symmetrize('a^2*b',variables='a,b,c')=
	a**2*b + a**2*c + a*b**2 + a*c**2 + b**2*c + b*c**2"""
	formula=S(formula)
	addsymbols(formula)
	if variables==None:
		variables=sorted(formula.free_symbols,key=str)
	else:
		variables=S(variables)
	for i in range(1,len(variables)):
		formula=cyclize(formula,oper,variables[:i+1])
	return formula
def findvalues(formula,values=None,variables=None,**kwargs):
	"""finds a candidate for parameter "values" in "prove" function
	Arguments:
	formula - sympy object or a string;
	values - list of numbers, starting point for fmin function,
	by default values=[1]*len(variables);
	variables - list of variables corresponding to values,
	by default variables=sorted(formula.free_symbols,key=str);
	**kwargs - list of parameters sent to fmin function
	>>> findvalues('a^2+b^2+c^2+d^2-a*(b+c+d)')
	Optimization terminated successfully.
			 Current function value: 1.154701
			 Iterations: 68
			 Function evaluations: 127
	(1.4339109663193974, 0.8278441585048405, 0.8279027492686561, 0.8278930696996669)
	>>> findvalues('a^2+b^2+c^2+d^2-a*(b+c+d)',disp=0)
	(1.4339109663193974, 0.8278441585048405, 0.8279027492686561, 0.8278930696996669)
	"""
	formula=S(formula)
	addsymbols(formula)
	num,den=fractioncancel(formula)
	if variables==None:
		variables=sorted(num.free_symbols,key=str)
	num=num.subs(zip(variables,list(map(lambda x:x**2,variables))))
	num=Poly(num)
	newformula=S((num.abs()+num)/(num.abs()-num))
	f=lambdify(variables,newformula)
	f2=lambda x:f(*x)
	if values==None:
		values=[1.0]*len(variables)
	else:
		values=S(values)
	tup=tuple(fmin(f2,values,**kwargs))
	return tuple([x*x for x in tup])
