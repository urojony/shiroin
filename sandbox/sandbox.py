from shirod import *
t=[0]*5
for shiroSeed in range(1,6):
	print(shiroSeed)
	#formula0=Sm('x^4/(8x^3+5y^3)+y^4/(8y^3+5z^3)+z^4/(8z^3+5x^3)-(x+y+z)/13')
	formula=Sm('x^4/(ax^3+by^3)+y^4/(ay^3+bz^3)+z^4/(az^3+bx^3)-(x+y+z)/(a+b)')
	formula2=formula.subs([['a',8],['b',5]])
	t[prove(makesubs(formula2,'[z,inf],[y,inf]',variables='x,z'))]+=1
	print('\n\n')
# ~ print(t)
#prove(makesubs(formula2,'[y,inf],[z,inf]',variables='x,y'))
# ~ #$$1630xy^6z^3 \le 326x^2y^8+978xy^5z^4+326y^7z^3$$
# ~ left=Sm('13890x^3y^5z^2+27780x^3y^4z^3+15100x^3y^3z^4+12000x^2y^6z^2+32910x^2y^5z^3+29490x^2y^4z^4+2205x^2y^3z^5')
# ~ right1=Sm('325x^7y^3+375x^7y^2z+375x^7yz^2+125x^7z^3+1755x^6y^4+3180x^6y^3z+3450x^6y^2z^2+2500x^6yz^3+675x^6z^4+3705x^5y^5+8085x^5y^4z+9330x^5y^3z^2')
# ~ right2=Sm('+9450x^5y^2z^3+5925x^5yz^4+1425x^5z^5+3770x^4y^6+7635x^4y^5z+3810x^4y^4z^2+3895x^4y^3z^3+8625x^4y^2z^4+6450x^4yz^5+1575x^4z^6+1950x^3y^7')
# ~ right3=Sm('+1590x^3y^6z+2955x^3y^2z^5+5305x^3yz^6+1375x^3z^7+454x^2y^8+570x^2y^7z+12180x^2y^2z^6+7440x^2yz^7+1425x^2z^8+780xy^8z+960xy^7z^2')
# ~ right4=Sm('+1872xy^5z^4+20145xy^4z^5+30465xy^3z^6+21420xy^2z^7+7515xyz^8+1075xz^9+780y^8z^2+3964y^7z^3+11960y^6z^4+21645y^5z^5+26130y^4z^6+20475y^3z^7+9945y^2z^8+2730yz^9+325z^10')
# ~ prove(right1+right2+right3+right4-left)
#prove(Sm('x^4/(8x^3+5y^3)-28x/169+15y/169'))

# ~ prove(Sm('-(27780x^3y^4z^3+32910x^2y^5z^3)+325x^7y^3+3180x^6y^3z+3160x^6y^2z^2+2500x^6yz^3+675x^6z^4+2566x^5y^5+4323x^5y^4z+1756x^5y^3z^2+6205x^5y^2z^3+5925x^5yz^4+1425x^5z^5+2675x^4y^6+3810x^4y^4z^2+8460x^4y^2z^4+5285x^4yz^5+1575x^4z^6+1950x^3y^7+1590x^3y^6z+2955x^3y^2z^5+5305x^3yz^6+1375x^3z^7+12055x^2y^2z^6+7440x^2yz^7+1425x^2z^8+2850xy^5z^4+1660xy^4z^5+30465xy^3z^6+21420xy^2z^7+7515xyz^8+1075xz^9+4290y^7z^3+21645y^5z^5+18850y^4z^6+20475y^3z^7+9945y^2z^8+2730yz^9+325z^10'))
class proof:
	m=0
from itertools import permutations
print(list(permutations(['x','y','z']))[1:3])
from operator import add,mul
from functools import reduce
print(add(5,8))
print(list(range(8,5)))
print(cyclize(Sm('x+2y+3z')),symmetrize(Sm('x+2y+3z')))
x=proof()
x.m=5
x.u=8
x.v=10
print(x.m,x.u,x,x.v)
