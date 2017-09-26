"""
This script tests fmin_slsqp using Example 14.4 from Numerical Methods for
Engineers by Steven Chapra and Raymond Canale.  This example maximizes the
function f(x) = 2*x0*x1 + 2*x0 - x0**2 - 2*x1**2, which has a maximum
at x0=2, x1=1.
"""

from scipy.optimize import fmin_slsqp
from numpy import array

def testfunc(x, *args):
    """
    Parameters
    ----------
    x : list
        A list of two elements, where x[0] represents x and
        x[1] represents y in the following equation.
    args : tuple
        First element of args is a multiplier for f.
        Since the objective function should be maximized, and the scipy
        optimizers can only minimize functions, it is nessessary to
        multiply the objective function by -1 to achieve the desired
        solution.
    Returns
    -------
    res : float
        The result, equal to ``2*x*y + 2*x - x**2 - 2*y**2``.

    """
    try:
        sign = args[0]
    except:
        sign = 1.0
    return sign*(2*x[0]*x[1] + 2*x[0] - x[0]**2 - 2*x[1]**2)

def testfunc_deriv(x,*args):
    """ This is the derivative of testfunc, returning a numpy array
    representing df/dx and df/dy """
    try:
        sign = args[0]
    except:
        sign = 1.0
    dfdx0 = sign*(-2*x[0] + 2*x[1] + 2)
    dfdx1 = sign*(2*x[0] - 4*x[1])
    return array([ dfdx0, dfdx1 ])

def test_eqcons(x,*args):
    """ Lefthandside of the equality constraint """
    return array([ x[0]**3-x[1] ])

def test_ieqcons(x,*args):
    """ Lefthandside of inequality constraint """
    return array([ x[1]-1 ])

def test_fprime_eqcons(x,*args):
    """ First derivative of equality constraint """
    return array([ 3.0*(x[0]**2.0), -1.0 ])

def test_fprime_ieqcons(x,*args):
    """ First derivative of inequality constraint """
    return array([ 0.0, 1.0 ])

from time import time

print "Unbounded optimization."
print "Derivatives of objective function approximated."
t0 = time()
result = fmin_slsqp(testfunc, [-1.0,1.0], args=(-1.0,), iprint=2, full_output=1)
print "Elapsed time:", 1000*(time()-t0), "ms"
print "Results", result, "\n\n"

print "Unbounded optimization."
print "Derivatives of objective function provided."
t0 = time()
result = fmin_slsqp(testfunc, [-1.0,1.0], fprime=testfunc_deriv, args=(-1.0,),
               iprint=2, full_output=1)
print "Elapsed time:", 1000*(time()-t0), "ms"
print "Results", result, "\n\n"

print "Bound optimization (equality constraints)."
print "Constraints implemented via lambda function."
print "Derivatives of objective function approximated."
print "Derivatives of constraints approximated."
t0 = time()
result = fmin_slsqp(testfunc, [-1.0,1.0], args=(-1.0,),
               eqcons=[lambda x, args: x[0]-x[1] ], iprint=2, full_output=1)
print "Elapsed time:", 1000*(time()-t0), "ms"
print "Results", result, "\n\n"

print "Bound optimization (equality constraints)."
print "Constraints implemented via lambda."
print "Derivatives of objective function provided."
print "Derivatives of constraints approximated."
t0 = time()
result = fmin_slsqp(testfunc, [-1.0,1.0], fprime=testfunc_deriv, args=(-1.0,),
               eqcons=[lambda x, args: x[0]-x[1] ], iprint=2, full_output=1)
print "Elapsed time:", 1000*(time()-t0), "ms"
print "Results", result, "\n\n"

print "Bound optimization (equality and inequality constraints)."
print "Constraints implemented via lambda."
print "Derivatives of objective function provided."
print "Derivatives of constraints approximated."
t0 = time()
result = fmin_slsqp(testfunc,[-1.0,1.0], fprime=testfunc_deriv, args=(-1.0,),
               eqcons=[lambda x, args: x[0]-x[1] ],
               ieqcons=[lambda x, args: x[0]-.5], iprint=2, full_output=1)
print "Elapsed time:", 1000*(time()-t0), "ms"
print "Results", result, "\n\n"

print "Bound optimization (equality and inequality constraints)."
print "Constraints implemented via function."
print "Derivatives of objective function provided."
print "Derivatives of constraints approximated."
t0 = time()
result = fmin_slsqp(testfunc, [-1.0,1.0], fprime=testfunc_deriv, args=(-1.0,),
               f_eqcons=test_eqcons, f_ieqcons=test_ieqcons,
               iprint=2, full_output=1)
print "Elapsed time:", 1000*(time()-t0), "ms"
print "Results", result, "\n\n"

print "Bound optimization (equality and inequality constraints)."
print "Constraints implemented via function."
print "All derivatives provided."
t0 = time()
result = fmin_slsqp(testfunc,[-1.0,1.0], fprime=testfunc_deriv, args=(-1.0,),
               f_eqcons=test_eqcons, fprime_eqcons=test_fprime_eqcons,
               f_ieqcons=test_ieqcons, fprime_ieqcons=test_fprime_ieqcons,
               iprint=2, full_output=1)
print "Elapsed time:", 1000*(time()-t0), "ms"
print "Results", result, "\n\n"
