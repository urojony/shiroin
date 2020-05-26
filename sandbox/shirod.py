import warnings,operator
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
"Try to set higher linprogiter parameter.",
"It looks like the formula is symmetric. You can assume without loss of"+
" generality that ","Try"
]
class proof:
	pass
for phrase in translationList:
	shiroTranslation[phrase]=phrase
from scipy.optimize import linprog
import random
from sympy import S,cancel,fraction,Pow,expand,solve,latex
import re
from collections.abc import Iterable
def blank(*args,**kwargs):
	return
def merge(ob1, ob2):
    ob1.__dict__.update(ob2.__dict__)
    return ob1
def _remzero(coef,fun): #removes addends with coefficient equal to zero
	ncoef=[]
	nfun=[]
	for c,f in zip(coef,fun):
		if c>0:
			ncoef+=[c]
			nfun+=[f]
	return ncoef,nfun
def slatex(formula):
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
def _input2fraction(formula,variables,values,translation,vprint):
	formula=S(formula)
	subst=[]
	for x,y in zip(variables,values):
		if y!=1:
			vprint(translation['Substitute']+' $ '+x+' \\to '+slatex(S(y)*S(x))+' $')
			subst+=[(x,x*y)]
	formula=formula.subs(subst)
	numerator,denominator=fractioncancel(formula)
	vprint(translation['numerator:']+' $$'+slatex(numerator)+'$$')
	vprint(translation['denominator:']+' $$'+slatex(denominator)+'$$')
	x=proof()
	x.subs=subst
	x.numerator=numerator
	x.denominator=denominator
	return x
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
def _list2proof(lcoef,lfun,rcoef,rfun,variables,itermax,linprogiter,translation,vprint,_writ2=_writ2):
	localseed=shiroSeed
	bufer=''
	itern=0
	if len(lcoef)==0:
		vprint(translation['status:'], 0)
		status=0
	elif len(rcoef)==0:
		vprint(translation['status:'], 2)
		status=2
		itermax=0
	ok=0
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
		status=res.status
		#print(lcoef)
		if itern==1:
			vprint(translation['status:']+' '+status)
			if status==0:
				print(translation['From weighted AM-GM inequality:'])
		if status==2:
			if ok==0:
				break
			else:
				lcoef,lfun=oldlcoef,oldlfun
				rcoef,rfun=oldrcoef,oldrfun
				bufer=''
				continue
		if status==0:
			vprint(bufer,end='')
			ok=1
			bufer=''
		oldlfun,oldrfun=lfun,rfun
		oldlcoef,oldrcoef=lcoef[:],rcoef[:]
		for i in range(m):
			c=0
			for j in res.x[i*n:i*n+n]:
				if(abs(round(j)-j)>0.0001):
					break
			else:
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
	vprint(bufer)
	x.inequalities+=[bufer]
	lhs='+'.join([_writ2(c,f,variables) for c,f in zip(lcoef,lfun)])
	if lhs=='':
		lhs='0'
	elif status==0:
		vprint(translation[
		"Program couldn't find a solution with integer coefficients. Try "+
		"to multiple the formula by some integer and run this function again."])
	elif(status==2):
		vprint(translation["Program couldn't find any proof."])
		#return res.status
	elif status==1:
		vprint(translation["Try to set higher linprogiter parameter."])
	bufer='$$ '+slatex(lhs)+' \\le '
	vprint(bufer)
	rhs='+'.join([_writ2(c,f,variables) for c,f in zip(rcoef,rfun)])
	if rhs=='':
		rhs='0'
	bufer2=slatex(rhs)+' $$'
	vprint(bufer2)
	x.inequalities+=[bufer+bufer2]
	if lhs=='0':
		vprint(translation['Sum of all inequalities gives us a desired proof.'])
	x.status=status
	return x
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
def prove(formula,values=None,variables=None,niter=200,linprogiter=10000,translation={},verbose=1):
	if verbose:vprint=print
	else:vprint=blank
	formula=S(formula)
	if variables: variables=_smakeiterable(variables)
	else: variables=sorted(formula.free_symbols,key=str)
	if values: values=_smakeiterable(values)
	else: values=[1]*len(variables)
	translation=_transcomplete(translation)
	x=_input2fraction(formula,variables,values,translation,vprint)
	st=_list2proof(*(_formula2list(x.numerator,variables)),variables,niter,linprogiter,translation,vprint)
	x=merge(x,st)
	if st.status==2 and issymetric(num):
		fs=sorted([str(x) for x in num.free_symbols])
		vprint(translation["It looks like the formula is symmetric. "+
		"You can assume without loss of generality that "],
		' >= '.join([str(x) for x in fs]),translation['Try'])
		vprint('prove(makesubs(S("',num,'"),',
		[(str(x),'inf') for x in variables[1:]],')')
	if verbose==1:
		return
	return x
def powerprove(formula,values=None,variables=None,niter=200,linprogiter=10000,translation={},verbose=1):
	if verbose:vprint=print
	else:vprint=blank
	formula=S(formula)
	if variables: variables=_smakeiterable(variables)
	else: variables=sorted(formula.free_symbols,key=str)
	if values: values=_smakeiterable(values)
	else: values=[1]*len(variables)
	translation=_transcomplete(translation)
	x=_input2fraction(formula,variables,values,translation)
	subst2=[]
	for j in range(len(variables)):
		subst2+=[(variables[j],1+variables[j])]
	for i in range(1<<len(variables)): #tricky substitutions to improve speed
		vprint('\n\\hline\n')
		subst1=[]
		substout=[]
		for j in range(len(variables)):
			if i&(1<<j):
				subst1+=[(variables[j],1/variables[j])]
				substout+=[str(variables[j])+'\\to 1/(1+'+str(variables[j])+')']
			else:
				substout+=[str(variables[j])+'\\to 1+'+str(variables[j])]
		vprint(translation['Substitute']+' $'+','.join(substout)+'$')
		num1=fractioncancel(x.num.subs(subst1))[0]
		num2=expand(num1.subs(subst2))
		vprint(translation["Numerator after substitutions:"]+' '+slatex(num2))
		y=_list2proof(*(_formula2list(num2,variables)),variables,niter,linprogiter,translation,vprint)
		prof+=[
def makesubs(formula,intervals,values=None,variables=None,numden=False,translation={},verbose=1):
	if verbose:vprint=print
	else:vprint=blank
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
		closed,opn=interval
		if {opn,closed}=={S('-inf'),S('inf')}:
			formula=formula.subs(var,var-1/var)
			equations=[equation.subs(var,var-1/var) for equation in equations]
			vprint(translation['Substitute'],'$',var,'\\to',var-1/var,'$')
		elif opn==S('inf'):
			formula=formula.subs(var,closed+var)
			equations=[equation.subs(var,closed+var) for equation in equations]
			vprint(translation['Substitute'],'$',var,'\\to',sstr(closed+var),'$')
		elif opn==S('-inf'):
			formula=formula.subs(var,closed-var)
			equations=[equation.subs(var,closed-var) for equation in equations]
			vprint(translation['Substitute'], "$",var,'\\to',sstr(closed-var),'$')
		else:
			formula=formula.subs(var,opn+(closed-opn)/var)
			vprint(translation['Substitute'], "$",var,'\\to',sstr(opn+(closed-opn)/(1+var)),'$')
			equations=[equation.subs(var,opn+(closed-opn)/var) for equation in equations]
	num,den=fractioncancel(formula)
	for var,interval in zip(variables,intervals):
		if {interval[0],interval[1]} & {S('inf'),S('-inf')}==set():
			num=num.subs(var,var+1)
			den=den.subs(var,var+1)
			equations=[equation.subs(var,var+1) for equation in equations]
	#print(equations)
	if values:
		values=ssolve(equations,variables)
		if len(values):
			vaues=values[0]
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
def provef(formula,niter=200,linprogiter=10000,translation={},verbose=True):
	if verbose:vprint=print
	else:vprint=blank
	formula=S(formula)
	translation=_transcomplete(translation)
	num,den=_input2fraction(formula,[],[],translation)
	_list2proof(*(_formula2listf(num)),None,niter,linprogiter,translation,vprint,_writ)
def issymetric(formula):
	if len(formula.free_symbols)<2:
		return True
	ls=list(formula.free_symbols)
	a=ls[0]
	for b in ls[1:]:
		if expand(formula-formula.subs({a:b, b:a}, simultaneous=True))!=S(0):
			return False
	return True
def cyclize(formula,oper=operator.add,variables=None,init=None):
	formula=S(formula)
	if variables==None:
		variables=sorted(formula.free_symbols,key=str)
	if len(variables)==0:
		return init
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
	if variables==None:
		variables=sorted(formula.free_symbols,key=str)
	for i in range(1,len(variables)):
		formula=cyclize(formula,oper,variables[:i+1])
	return formula
