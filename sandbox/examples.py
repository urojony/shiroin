from shiroindev import *
from sympy import powdenest
shiroSeed=1
###https://www.imomath.com/index.php?options=593&lmm=0
#Problem 1
prove('(a^2+b^2+c^2-a*b-a*c-b*c)*2')
#Problem 2
formula,values=Sm('(a^2+b^2+c^2+d^2-a(b+c+d))'),'2,1,1,1'
prove(formula,values)
#Problem 3
formula=Sm('(ab/c+bc/a+ca/b)-1')
formula=formula.subs('c',S('1-a-b'))
values,intervals='1/3,1/3','[0,1-b],[0,1]'
formula,values=makesubs(formula,intervals,values)
prove(formula,values)
#Problem 4
formula,values=Sm('1/(1-x^2)+1/(1-y^2)-2/(1-xy)'),'1/2,1/2'
formula*=3
intervals='[1,0],[1,0]' #intervals [1,0] and [0,1] are equivalent
formula,values=makesubs(formula,intervals,values)
prove(formula,values)
#Problem 5
formula=Sm('a^3+b^3-a^2b-ab^2')
formula*=3
prove(formula)
#Problem 7
formula=Sm('a/(b+c)+b/(c+a)+c/(a+b)-3/2')
formula*=3
prove(formula)
#Problem 10
formula=Sm('a/(b+c)+b/(c+d)+c/(d+a)+d/(a+b)-2')
#without loss of generality, a>=c and b>=d
prove(makesubs(formula*2,'[c,inf],[d,inf]')) 
#Problem 11
formula=Sm('a^3/(a^2+ab+b^2)+b^3/(b^2+bc+c^2)+c^3/(c^2+ca+a^2)-(a+b+c)/3')
formula*=12
prove(formula)
#Problem 13
formula=Sm('(a^2+c^2+e^2)(b^2+d^2+f^2)-(ab+cd+ef)^2')
prove(formula)
#Problem 14
formula=Sm('(5a^3-ab^2)/(a+b)+(5b^3-bc^2)/(b+c)+(5c^3-ca^2)/(c+a)-2(a^2+b^2+c^2)')
prove(formula)


# ~ ###https://www.imomath.com/index.php?options=608&lmm=0
#Problem 1
sub=S('[[a,x^2+1],[b,y^2+1],[c,z^2+1]]')
left=S('(a-1)^(1/2)+(b-1)^(1/2)+(c-1)^(1/2)').subs(sub)
right=S('(c*a*b+c)^(1/2)').subs(sub)
formula=powdenest(right**2-left**2,force=True) #getting rid of roots in formula
prove(makesubs(formula,'[0,x/(x^2+1)],[0,1/x]',variables='z,y,x'))
print('\n\n')
prove(makesubs(formula,'[x/(x^2+1),inf],[0,1/x]',variables='z,y,x'))
print('\n\n')
prove(makesubs(formula,'[0,x/(x^2+1)],[1/x,inf]',variables='z,y,x'))
print('\n\n')
prove(makesubs(formula,'[x/(x^2+1),inf],[1/x,inf]',variables='z,y,x'))
#Problem 3
formula=Sm('(1/a-1)(1/b-1)(1/c-1)-8')
formula=formula.subs('c',S('1-a-b'))
num,den,values=makesubs(formula,'[1-b,0],[1,0]','1/3,1/3',numden=True)
#print(den)
prove(num,values)
#Problem 5
formula=Sm('s^2-ay-bz-cx').subs([['x','s-a'],['y','s-b'],['z','s-c']])
prove(makesubs(formula,'[0,s],[0,s],[0,s]'))
#Problem 6
left=S('u^2+v^2+w^2+x^2+y^2')
right=S('2/sqrt(3)*(u*v+v*w+w*x+x*y)')
formula=left**2-right**2
prove(formula,'1,sqrt(3),2,sqrt(3),1')
#Problem 17
formula=Sm('(a+b)/(2b+c)+(b+c)/(2c+a)+(c+a)/(2a+b)-2')
prove(formula*3)

#https://math.stackexchange.com/questions/1775572/olympiad-inequality-sum-limits-cyc-fracx48x35y3-geqslant-fracxy?fbclid=IwAR3Nr2-aWTHBTXGMBeR7ww8Z4OtUH4m1fpZZAi1RsajsTU4pXtlf1Gcba3A
# ~ formula=Sm('x^4/(8x^3+5y^3)+y^4/(8y^3+5z^3)+z^4/(8z^3+5x^3)-(x+y+z)/13')
# ~ prove(makesubs(formula,'[z,inf],[y,inf]',variables='x,z'))
# ~ print('\n\n')
# ~ prove(makesubs(formula,'[y,inf],[z,inf]'))

###other problems
formula=Sm('-(3a + 2b + c)(2a^3 + 3b^2 + 6c + 1) + (4a + 4b + 4c)(a^4 + b^3 + c^2 + 3)')
powerprove(formula)

formula=Sm('-80a^2b^3c^2+3a^3b^4+20b^3c^4+20b^4c^3+26a^5b^2+28a^5c^2+15b^2c^5+1b^5c^2+4a^5bc')
prove(formula)

prove(Sm('a^2-2a'),'1','a')

formula=Sm('a^3/(a^2+ab+b^2)+b^3/(b^2+bc+c^2)+c^3/(c^2+ca+a^2)-(a+b+c)/3')
formula*=6
prove(formula,'1,1,1','a,b,c',linprogiter=3)

provef('3*f(6,0)-3*f(4,2)-3*f(2,4)+4*f(0,6)')
