.. py:function:: Sm(formula)

   Convert the formula to a SymPy object assuming that there are no functions nor multiletter variables (add multiplication sign between any digit/letter/")" and letter/"(", then run sympy.sympify)
  
Example:

	>>> Sm('(x^2+5)(yz^2+5)')
	(x**2 + 5)*(y*z**2 + 5)
