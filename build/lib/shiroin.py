import warnings
warnings.filterwarnings("ignore")
shiroSeed=None
shiroTranslation={}
translationList=['numerator:','denominator:','status:',
'Substitute',"Formula after substitution:",
"Numerator after substitutions:","From weighted AM-GM inequality:",
'Sum of all inequalities gives us a desired proof.',
"Program couldn't find a solution with integer coefficients. Try "+
"to multiple the formula by some integer and run this function again.",
"Program couldn't find any proof.",
"Try to set higher linprogiter parameter."
]
for phrase in translationList:
	shiroTranslation[phrase]=phrase
from scipy.optimize import linprog
import random
from sympy import S,cancel,fraction,Pow,expand,solve,latex
import re
from collections.abc import Iterable
def _remzero(coef,fun): #removes addends with coefficient equal to zero
	ncoef=[]
	nfun=[]
	for c,f in zip(coef,fun):
		if c>0:
			ncoef+=[c]
			nfun+=[f]
	return ncoef,nfun
def slatex(formula):
	return re.sub('{(.)}',r'\1',latex(formula,fold_short_frac=True).replace(' ',''))
def _transcomplete(translation):
	trans=shiroTranslation.copy()
	trans.update(translation)
	return trans
def _writ2(coef,fun,variables):
	return slatex(S((str(coef)+'*'+'*'.join([str(x)+'^'+str(y) for x,y in zip(variables,fun)]))))
def _check(coef,fun,res,rfun): #checks if rounding and all floating point stuff works
	res2=[int(round(x)) for x in res]
	b1=[coef*x for x in fun]
	b2=[[x*y for y in rfuni] for x,rfuni in zip(res2,rfun)]
	return b1==[sum(x) for x in zip(*b2)] and coef==sum(res2)
def _powr(formula):
	if formula.func==Pow:
		return formula.args
	else:
		return [formula,S('1')]
def fractioncancel(formula): #workaround for buggy cancel function
	num,den=fraction(cancel(formula/S('tmp')))
	den=den.subs('tmp','1')
	return num,den
def ssolve(formula,variables): #workaround for inconsistent solve function
	result=solve(formula,variables)
	if type(result)==dict:
		result=[[result[var] for var in variables]]
	return result
def sstr(formula):
	return str(formula).replace('**','^').replace('*','').replace(' ','')
def Sm(formula):
	if type(formula)==str:
		formula.replace(' ','')
		for i in range(2):
			formula=re.sub(r'([0-9a-zA-Z)])([(a-zA-Z])',r'\1*\2',formula)
	formula=S(formula)
	return formula
def _input2fraction(formula,variables,values,translation):
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
	return numerator,denominator
def _formula2list(formula,variables):
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
def _list2proof(lcoef,lfun,rcoef,rfun,variables,itermax,linprogiter,translation):
	localseed=shiroSeed
	itern=0
	if len(lcoef)==0:
		print(translation['status:'], 0)
	elif len(rcoef)==0:
		print(translation['status:'], 2)
		itermax=0
	while len(lcoef)>0 and itern<itermax:
		itern+=1
		m=len(lcoef)
		n=len(rcoef)
		lfunt=list(map(list, zip(*lfun)))
		rfunt=list(map(list, zip(*rfun)))
		A=[]
		b=[]
		A_ub=[]
		b_ub=[]
		for i in range(m):
			for j in range(len(rfunt)):
				A+=[[0]*(m*n)]
				A[-1][i*n:i*n+n]=rfunt[j]
				b+=[lfun[i][j]*lcoef[i]]
			A+=[[0]*(m*n)]
			A[-1][i*n:i*n+n]=[1]*n
			b+=[lcoef[i]]
		for j in range(n):
			A_ub+=[[0]*(m*n)]
			A_ub[-1][j::n]=[1]*m
			b_ub+=[rcoef[j]]
		random.seed(localseed)
		vecc=[random.random() for i in range(m*n)]
		localseed=random.randint(1,1000000000)
		res=linprog(vecc,A_eq=A,b_eq=b,A_ub=A_ub,b_ub=b_ub,options={'maxiter':linprogiter})
		if itern==1:
			print(translation['status:'],res.status)
			if res.status==0:
				print(translation['From weighted AM-GM inequality:'])
		if res.status>1:
			break
		for i in range(m):
			c=0
			for j in res.x[i*n:i*n+n]:
				if(abs(round(j)-j)>0.0001):
					break
			else:
				isok=_check(lcoef[i],lfun[i],res.x[i*n:i*n+n],rfun)
				if not isok:
					continue
				print('$$',_writ2(lcoef[i],lfun[i],variables),' \\le ',sep='',end='')
				lcoef[i]=0
				for j in range(n):
					rcoef[j]-=int(round(res.x[i*n+j]))
				for j,k in zip(res.x[i*n:i*n+n],rfun):
					if j<0.0001:
						continue
					if(c):print('+',end='')
					else:c=1
					print(_writ2(int(round(j)),k,variables),end='')
				print('$$')
		lcoef,lfun=_remzero(lcoef,lfun)
		rcoef,rfun=_remzero(rcoef,rfun)
	print()
	lhs='+'.join([_writ2(c,f,variables) for c,f in zip(lcoef,lfun)])
	if lhs=='':
		lhs='0'
	elif res.status==0:
		print(translation[
		"Program couldn't find a solution with integer coefficients. Try "+
		"to multiple the formula by some integer and run this function again."])
	elif(res.status==2):
		print(translation["Program couldn't find any proof."])
		return
	elif res.status==1:
		print(translation["Try to set higher linprogiter parameter."])
	print('$$ ',slatex(lhs),' \\le ')
	rhs='+'.join([_writ2(c,f,variables) for c,f in zip(rcoef,rfun)])
	if rhs=='':
		rhs='0'
	print(slatex(rhs),' $$')
	if lhs=='0':
		print(translation['Sum of all inequalities gives us a desired proof.'])
def _smakeiterable(x):
	x=S(x)
	if isinstance(x,Iterable):
		return x
	return (x,)
def _smakeiterable2(x):
	x=S(x)
	if isinstance(x[0],Iterable):
		return x
	return (x,)
def prove(formula,values=None,variables=None,niter=200,linprogiter=10000,translation={}):
	formula=S(formula)
	if variables: variables=_smakeiterable(variables)
	else: variables=sorted(formula.free_symbols,key=str)
	if values: values=_smakeiterable(values)
	else: values=[1]*len(variables)
	translation=_transcomplete(translation)
	num,den=_input2fraction(formula,variables,values,translation)
	_list2proof(*(_formula2list(num,variables)),variables,niter,linprogiter,translation)
def powerprove(formula,values=None,variables=None,niter=200,linprogiter=10000,translation={}):
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
		_list2proof(*(_formula2list(num2,variables)),variables,niter,linprogiter,translation)
def makesubs(formula,intervals,values=None,variables=None,numden=False,translation={}):
	formula=S(formula)
	translation=_transcomplete(translation)
	intervals=_smakeiterable2(intervals)
	if variables: variables=_smakeiterable(variables)
	else: variables=sorted(formula.free_symbols,key=str)
	if values:
		values=_smakeiterable(values)
		equations=[var-value for var,value in zip(variables,values)]
	else:
		equations=[]
	for var,interval in zip(variables,intervals):
		closed,opn=interval
		if {opn,closed}=={S('-inf'),S('inf')}:
			formula=formula.subs(var,var-1/var)
			equations=[equation.subs(var,var-1/var) for equation in equations]
			print(translation['Substitute'],'$',var,'\\to',var-1/var,'$')
		elif opn==S('inf'):
			formula=formula.subs(var,closed+var)
			equations=[equation.subs(var,closed+var) for equation in equations]
			print(translation['Substitute'],'$',var,'\\to',sstr(closed+var),'$')
		elif opn==S('-inf'):
			formula=formula.subs(var,closed-var)
			equations=[equation.subs(var,closed-var) for equation in equations]
			print(translation['Substitute'], "$",var,'\\to',sstr(closed-var),'$')
		else:
			formula=formula.subs(var,opn+(closed-opn)/var)
			print(translation['Substitute'], "$",var,'\\to',sstr(opn+(closed-opn)/(1+var)),'$')
			equations=[equation.subs(var,opn+(closed-opn)/var) for equation in equations]
	num,den=fractioncancel(formula)
	for var,interval in zip(variables,intervals):
		if {interval[0],interval[1]} & {S('inf'),S('-inf')}==set():
			num=num.subs(var,var+1)
			den=den.subs(var,var+1)
			equations=[equation.subs(var,var+1) for equation in equations]
	#print(equations)
	values=ssolve(equations,variables)[0]
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
