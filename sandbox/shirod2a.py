from __future__ import print_function
import warnings,operator
warnings.filterwarnings("ignore")    
#Seed is needed to select the weights in linprog function.
#None means that the seed is random.
shiroSeed=None
shiroTranslation={}
translationList=['numerator:','denominator:','status:',
'Substitute',"Formula after substitution:",
"Numerator after substitutions:","From weighted AM-GM inequality:",
'Sum of all inequalities gives us a desired proof.',
"Program couldn't find a solution with integer coefficients. Try "+
"to multiple the formula by some integer and run this function again.",
"Program couldn't find any proof.",
"Try to set higher linprogiter parameter.",
"It looks like the formula is symmetric. You can assume without loss of"+
" generality that ","Try"
]
#Initialize english-english dictionary.
for phrase in translationList:
	shiroTranslation[phrase]=phrase
from scipy.optimize import linprog
import random
from sympy import S,cancel,fraction,Pow,expand,solve,latex
import re
def _remzero(coef,fun): 
#coef, fun represents an expression.
#For example, if expression=5f(2,3)+0f(4,6)+8f(1,4)
#then coef=[5,0,8], fun=[[2,3],[4,6],[1,4]]
#_remzero removes addends with coefficient equal to zero.
#In this example ncoef=[5,8], nfun=[[2,3],[1,4]]
	ncoef=[]
	nfun=[]
	for c,f in zip(coef,fun):
		if c>0:
			ncoef+=[c]
			nfun+=[f]
	return ncoef,nfun
def slatex(formula): #fancy function which makes latex code more readable, but still correct
	formula=re.sub('\^{(.)}',r'^\1',latex(formula,fold_short_frac=True).replace(' ','').replace('\\left(','(').replace('\\right)',')'))
	return re.sub('\{(\(.+?\))\}',r'\1',formula)
def _transcomplete(translation):
	trans=shiroTranslation.copy()
	trans.update(translation)
	return trans
def _writ2(coef,fun,variables):
	return slatex(S((str(coef)+'*'+'*'.join([str(x)+'^'+str(y) for x,y in zip(variables,fun)]))))
def _writ(coef,fun,nullvar):
	return str(coef)+'f('+str(fun)[1:-1-(len(fun)==1)]+')'
def _check(coef,fun,res,rfun): 
#checks if rounding and all the floating point stuff works
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
#workaround for buggy cancel function
	num,den=fraction(cancel(formula/S('tmp')))
	den=den.subs('tmp','1')
	return num,den
def ssolve(formula,variables): 
#workaround for inconsistent solve function
	result=solve(formula,variables)
	if type(result)==dict:
		result=[[result[var] for var in variables]]
	return result
def sstr(formula):
	return str(formula).replace('**','^').replace('*','').replace(' ','')
def Sm(formula):
#Adds multiplication signs and sympifies a formula.
#For example, Sm('(2x+y)(7+5xz)') -> S('(2*x+y)*(7+5*x*z)')
	if type(formula)==str:
		formula.replace(' ','')
		for i in range(2):
			formula=re.sub(r'([0-9a-zA-Z)])([(a-zA-Z])',r'\1*\2',formula)
	formula=S(formula)
	return formula
def _input2fraction(formula,variables,values,translation):
#makes some substitutions and converts formula to a fraction
#with expanded numerator and denominator
	formula=S(formula)
	subst=[]
	for x,y in zip(variables,values):
		if y!=1:
			print(translation['Substitute'],'$',x,'\\to',slatex(S(y)*S(x)),'$')
			subst+=[(x,x*y)]
	formula=formula.subs(subst)
	numerator,denominator=fractioncancel(formula)
	print(translation['numerator:'],'$$'+slatex(numerator)+'$$')
	print(translation['denominator:'],'$$'+slatex(denominator)+'$$')
	return (numerator,denominator)
def _formula2list(formula,variables):
#Splits a polynomial to a difference of two polynomials with positive
#coefficients and extracts coefficients and powers of both polynomials.
#'variables' is used to set order of powers
#For example, If formula=5x^2-4xy+8y^3, variables=[x,y], then
#the program tries to prove that 
#0<=5x^2-4xy+8y^3
#4xy<=5x^2+8y^3
#lcoef=[4]
#lfun=[[1,1]]
#rcoef=[5,8]
#rfun=[[2,0],[0,3]]
	lfun=[]
	lcoef=[]
	rfun=[]
	rcoef=[]
	varorder=dict(zip(variables,range(len(variables))))
	for addend in formula.as_ordered_terms():
		coef,facts=addend.as_coeff_mul()
		powers=[0]*len(variables)
		for var in variables:
			powers[varorder[var]]=0
		for fact in facts:
			var,pw=_powr(fact)
			powers[varorder[var]]=int(pw)
		if(coef<0):
			lcoef+=[-coef]
			lfun+=[powers]
		else:
			rcoef+=[coef]
			rfun+=[powers]
	return(lcoef,lfun,rcoef,rfun)
def _list2proof(lcoef,lfun,rcoef,rfun,variables,itermax,linprogiter,translation,_writ2=_writ2):
#Now the formula is splitted on two polynomials with positive coefficients.
#we will call them LHS and RHS and our inequality to prove would
#be LHS<=RHS (instead of 0<=RHS-LHS).

#suppose we are trying to prove that 
#30x^2y^2+60xy^4<=48x^3+56y^6 (assuming x,y>0)
#program will try to find some a,b,c,d such that
#30x^2y^2<=ax^3+by^6  
#60xy^4<=cx^3+dy^6  
#where a+c<=48 and b+d<=56 (assumption 1)
#We need some additional equalities to meet assumptions
#of the weighted AM-GM inequality.
#a+b=30 and c+d=60 (assumption 2)
#3a+0b=30*2, 0a+6b=30*2, 3c+0d=60*1, 0c+6d=60*4 (assumption 3)

#The sketch of the algorithm.
# for i in range(itermax):
	#1. Create a vector of random numbers (weights).
	#2. Try to find real solution of the problem (with linprog).
	#3. If there is no solution (status: 2)
		#3a. If the solution was never found, break.
		#3b. Else, step back (to the bigger inequality)
	#4. If the soltuion was found (status: 0)
		#Check out which of variables (in example: a,b,c,d) looks like integer.
		#If there are some inequalities with all integer coefficients, subtract 
		#them from the original one.
		#If LHS is empty, then break. 
	localseed=shiroSeed
	bufer=''
	itern=0
	if len(lcoef)==0: #if LHS is empty
		print(translation['status:'], 0)
		status=0
	elif len(rcoef)==0: 
	#if RHS is empty, but LHS is not
		print(translation['status:'], 2)
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
			print(translation['status:'],status)
			if status==0:
				print(translation['From weighted AM-GM inequality:'])
		if status==2: #if real solution of current inequality doesn't exist
			if foundreal==0: #if this is the first inequality, then break 
				break
			else:
				#step back 
				lcoef,lfun=oldlcoef,oldlfun
				rcoef,rfun=oldrcoef,oldrfun
				bufer=''
				continue
		if status==0:#if found a solution with real coefficients
			print(bufer,end='')
			foundreal=1
			bufer=''
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
					bufer+='$$'+_writ2(lcoef[i],lfun[i],variables)+' \\le '
					lcoef[i]=0
					for j in range(n):
						rcoef[j]-=int(round(res.x[i*n+j]))
					for j,k in zip(res.x[i*n:i*n+n],rfun):
						if j<0.0001:
							continue
						if(c):bufer+='+'
						else:c=1
						bufer+=_writ2(int(round(j)),k,variables)
					bufer+='$$\n'
			lcoef,lfun=_remzero(lcoef,lfun)
			rcoef,rfun=_remzero(rcoef,rfun)
	print(bufer)
	lhs='+'.join([_writ2(c,f,variables) for c,f in zip(lcoef,lfun)])
	if lhs=='':
		lhs='0'
	elif status==0:
		print(translation[
		"Program couldn't find a solution with integer coefficients. Try "+
		"to multiple the formula by some integer and run this function again."])
	elif(status==2):
		print(translation["Program couldn't find any proof."])
		#return res.status
	elif status==1:
		print(translation["Try to set higher linprogiter parameter."])
	print('$$ ',slatex(lhs),' \\le ')
	rhs='+'.join([_writ2(c,f,variables) for c,f in zip(rcoef,rfun)])
	if rhs=='':
		rhs='0'
	print(slatex(rhs),' $$')
	if lhs=='0':
		print(translation['Sum of all inequalities gives us a desired proof.'])
	return status
def _isiterable(obj):
	try:
		_ = (e for e in obj)
		return True
	except TypeError:
		return False
def _smakeiterable(x):
	x=S(x)
	if _isiterable(x):
		return x
	return (x,)
def _smakeiterable2(x):
	x=S(x)
	if _isiterable(x[0]):
		return x
	return (x,)
def prove(formula,values=None,variables=None,niter=200,linprogiter=10000,translation={}):
#tries to prove that formula>=0 assuming all variables are positive
	formula=S(formula)
	if variables: variables=_smakeiterable(variables)
	else: variables=sorted(formula.free_symbols,key=str)
	if values: values=_smakeiterable(values)
	else: values=[1]*len(variables)
	translation=_transcomplete(translation)
	num,den=_input2fraction(formula,variables,values,translation)
	st=_list2proof(*(_formula2list(num,variables)+(variables,niter,linprogiter,translation)))
	if st==2 and issymetric(num):
		fs=sorted([str(x) for x in num.free_symbols])
		print(translation["It looks like the formula is symmetric. "+
		"You can assume without loss of generality that "],
		' >= '.join([str(x) for x in fs]),translation['Try'])
		print('prove(makesubs(S("',num,'"),',
		[(str(x),'inf') for x in variables[1:]],')')
	return st
def powerprove(formula,values=None,variables=None,niter=200,linprogiter=10000,translation={}):
#This is a bruteforce and ineffective function for proving inequalities.
#It can be used as the last resort.
	formula=S(formula)
	if variables: variables=_smakeiterable(variables)
	else: variables=sorted(formula.free_symbols,key=str)
	if values: values=_smakeiterable(values)
	else: values=[1]*len(variables)
	translation=_transcomplete(translation)
	num,den=_input2fraction(formula,variables,values,translation)
	subst2=[]
	for j in range(len(variables)):
		subst2+=[(variables[j],1+variables[j])]
	for i in range(1<<len(variables)): #tricky substitutions to improve speed
		print('\n\\hline\n')
		subst1=[]
		substout=[]
		for j in range(len(variables)):
			if i&(1<<j):
				subst1+=[(variables[j],1/variables[j])]
				substout+=[str(variables[j])+'\\to 1/(1+'+str(variables[j])+')']
			else:
				substout+=[str(variables[j])+'\\to 1+'+str(variables[j])]
		print(translation['Substitute'], '$'+','.join(substout)+'$')
		num1=fractioncancel(num.subs(subst1))[0]
		num2=expand(num1.subs(subst2))
		print(translation["Numerator after substitutions:"],slatex(num2))
		_list2proof(*(_formula2list(num2,variables)+(variables,niter,linprogiter,translation)))
def makesubs(formula,intervals,values=None,variables=None,numden=False,translation={}):
#This function generates a new formula which satisfies this condition:
#for all positive variables new formula is nonnegative iff
#for all variables in corresponding intervals old formula is nonnegative
	formula=S(formula)
	translation=_transcomplete(translation)
	intervals=_smakeiterable2(intervals)
	if variables: variables=_smakeiterable(variables)
	else: variables=sorted(formula.free_symbols,key=str)
	if values!=None:
		values=_smakeiterable(values)
		equations=[var-value for var,value in zip(variables,values)]
	else:
		equations=[]
	for var,interval in zip(variables,intervals):
		end1,end2=interval
		if end1 in {S('-inf'),S('inf')}:
			end1,end2=end2,end1
		if {end1,end2}=={S('-inf'),S('inf')}:
			formula=formula.subs(var,var-1/var)
			equations=[equation.subs(var,var-1/var) for equation in equations]
			print(translation['Substitute'],'$',var,'\\to',var-1/var,'$')
		elif end2==S('inf'):
			formula=formula.subs(var,end1+var)
			equations=[equation.subs(var,end1+var) for equation in equations]
			print(translation['Substitute'],'$',var,'\\to',sstr(end1+var),'$')
		elif end2==S('-inf'):
			formula=formula.subs(var,end1-var)
			equations=[equation.subs(var,end1-var) for equation in equations]
			print(translation['Substitute'], "$",var,'\\to',sstr(end1-var),'$')
		else:
			formula=formula.subs(var,end2+(end1-end2)/var)
			print(translation['Substitute'], "$",var,'\\to',sstr(end2+(end1-end2)/(1+var)),'$')
			equations=[equation.subs(var,end2+(end1-end2)/var) for equation in equations]
	num,den=fractioncancel(formula)
	for var,interval in zip(variables,intervals):
		if {interval[0],interval[1]} & {S('inf'),S('-inf')}==set():
			num=num.subs(var,var+1)
			den=den.subs(var,var+1)
			equations=[equation.subs(var,var+1) for equation in equations]
	if values:
		values=ssolve(equations,variables)
		if len(values):
			values=values[0]
	num,den=expand(num),expand(den)
	#print(translation["Formula after substitution:"],"$$",slatex(num/den),'$$')
	if values and numden:
		return num,den,values
	elif values:
		return num/den,values
	elif numden:
		return num,den
	else:
		return num/den
def _formula2listf(formula):
#Splits a polynomial to a difference of two formulas with positive
#coefficients and extracts coefficients and function
#arguments of both formulas.
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
	return(lcoef,lfun,rcoef,rfun)
def provef(formula,niter=200,linprogiter=10000,translation={}):
#this function is similar to prove, formula is a linear combination of
#values of f:R^k->R instead of a polynomial. provef checks if a formula
#is nonnegative for any nonnegative and convex function f. If so, it
#provides a proof of nonnegativity.
	formula=S(formula)
	translation=_transcomplete(translation)
	num,den=_input2fraction(formula,[],[],translation)
	_list2proof(*(_formula2listf(num)+(None,niter,linprogiter,translation,_writ)))
def issymetric(formula): #checks if formula is symmetric
#and has at least two variables	
	if len(formula.free_symbols)<2:
		return False
	ls=list(formula.free_symbols)
	a=ls[0]
	for b in ls[1:]:
		if expand(formula-formula.subs({a:b, b:a}, simultaneous=True))!=S(0):
			return False
	return True
def cyclize(formula,oper=operator.add,variables=None,init=None):
#cyclize('a^2*b')=S('a^2*b+b^2*a')
#cyclize('a^2*b',variables='a,b,c')=S('a^2*b+b^2*c+c^2*a')
	formula=S(formula)
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
#symmetrize('a^2*b')=S('a^2*b+b^2*a')
#symmetrize('a^2*b',variables='a,b,c')=
#=S('a^2*b+a^2*c+b^2*a+b^2*c+c^2*a+c^2*b')
	formula=S(formula)
	if variables==None:
		variables=sorted(formula.free_symbols,key=str)
	else:
		variables=S(variables)
	for i in range(1,len(variables)):
		formula=cyclize(formula,oper,variables[:i+1])
	return formula
