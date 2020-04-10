"""

Tool for general non-linear fitting based on the Marquardt algorithm. To setup a fit you need to::

   1. define parameters and assign them reasonable initial values
   2. define a fitting function in terms of the parameters
   3. define the fit and carry it out

Depending on the outcome, one might need to change the initial values of the parameters or
decide to change which parameters are fixed and which are to be fitted.

1. Setting up the parameters::

    >>> import LT_Fit.parameters as P  # get the parameter module
    >>> import LT_fit.gen_fit as G     # load the genfit module
    >>> a = P.Parameter(1., 'amp')     # parameter a, called 'amp', initialized to 1. 
    >>> b = P.Parameter(3., 'omega')   # parameter b, called 'omega', intialized to 3.
    >>> c = P.Parameter(1.5, 'phase')  # parameter c, called 'phase', initialize to 1.5

2. Define the fitting function::

    >>> def f(x):
    ....    return a()*sin(b() * x + c() )

3. Now do the fit::

    >>> F = G.gen_slsqp( f, [a, b, c], x = xval, y = y_exp) #  the list [a,b,c] defines which parameters are fitted
    >>> F = G.gen_slsqp( f, [a], x = xval, y = y_exp, yerr = sigy) #  only a is fitted, but take the exp. error into account
    >>> F = G.gen_slsq0( f, [a, c], x = xval, y = y_exp, yerr = sigy) #  a and c are fitted b is kept constant
    
To change the initial values and repeat the fit::

    >>> a.Set(1.5)
    >>> b.Set(3.2)
    >>> c.Set(2.0)
    >>> F = G.genfit( f, [a, b, c], x = xval, y = y_exp)

Finally plot the data and the fit::

    >>> import LT.box as B
    >>> B.plot_exp(xval, y_exp, sigy)   # plot the data with errorbars
    >>> B.plot_line(F.xpl, F.ypl)       # plot the fitted function as a line
    >>> show()                          # show the new plot


----------------------------------------------------------


"""
#
from scipy import optimize
from LT_Fit.parameters import *
from numpy import *
import pdb


class  gen_slsqp:
    """

    general non-linear fit  based on the Marquardt algorithm

    Important keywords:
    
    ==========   ===================================================================
    Keyword      Meaning
    ==========   ===================================================================
    y            (:func:`numpy.array`) array of experimental values (mandatory)
    x            (:func:`numpy.array`) array of independent variables 
    y_err        (:func:`numpy.array`) array of errors
    nplot        number of points to be used for plotting the fit
    ftol         minimal change in chi square to determine if the fit has converged
    ==========   ===================================================================

    Additional keyword arguments are passed on to :func:`scipy.optimize.leastsq`

    """
    def __init__(self, function, parameters, \
                 x = None, \
                 y = None, \
                 y_err = None, \
                 nplot = 100, \
                 full_output = 1, \
                 ftol = 0.001, \
                 print_results = True, \
                 **kwargs):
        self.print_results = print_results
        self.y = y
        if x is None:
            if y is not None:
                self.x = arange(y.shape[0])
            else:
                self.x = x
        else:
            self.x = x
        self.y_err = y_err
        self.nplot = nplot
        self.parameters = parameters # the array stores the addresses of the parameter objects
        self.func = function
        # carry out the fit
        if y is None:
            print 'No values to fit, use set_yval to set them before fitting !'
            return
        return self.fit(full_output = full_output, \
                 ftol = ftol, \
                 **kwargs)

    def f(self, params):
        # define minimzation function for least square        i = 0
        i = 0
        for p in self.parameters:
            # store the current values in the array p back into the parameter
            # classes to be used by the user defined function
            # check if there is only 1 parameter
            ps = params.shape
            if len(ps) > 0:
                p.set(params[i])
            else:
                p.set(params)
            i += 1
        # now calculate the difference between data and fit function
        # and return it to the optimization routine
        if self.y_err == None:
            return self.y - self.func(self.x)
        else:
            return (self.y - self.func(self.x))/self.y_err
    # end of the minimization function

    def fit(self, full_output = 1, ftol = 0.001,  **kwargs):
        # this is the minimzation routine
        # store the current parameter values in a list to be passed to the fitting function 
        # this is an implicit loop, it is equivalent to
        # p =[]
        # for param in parameters:
        #     p.append(param())
        # clear the parameter errors:
        for p in self.parameters:
            p.err = 0.
        p = [param() for param in self.parameters]
        self.fit_result = optimize.fmin_slsqp(self.f,\
                                              p, \
                                              iprint=2, \
                                              full_output=1,\
                                              **kwargs)
        #
        #self.fit_result = optimize.leastsq( self.f, p,\
        #                                    full_output = full_output, \
        #                                    ftol = ftol, \
        #                                    **kwargs)
        # now calculate the fitted values
        fit_func = self.func(self.x)
        # self.covar = self.fit_result[1]
        # final total chi square
        p_fin = self.fit_result[0]
        # number of degrees of freedom
        ps = p_fin.shape
        if len(ps) > 0:
            # if there was only 1 data point
            self.n_dof = len(self.y) - len(p_fin)
        else:
            self.n_dof = len(self.y)
        self.chi2 = sum( power( self.f( self.fit_result[0]) , 2) )
        self.chi2_red = self.chi2/self.n_dof
        self.xpl = []
        self.ypl = []
        self.stat = {'fitted values':fit_func, \
                     'parameters':self.fit_result[0], \
                     'leastsq output':self.fit_result}
        if (self.nplot > 0):
            self.xpl = linspace(self.x.min(), self.x.max(), self.nplot+1)
            self.ypl = self.func(self.xpl)
        # set the parameter errors
        try:
            for i,p in enumerate(self.parameters):
                p.err = sqrt( self.covar[i,i] )
        except:
            print "gen_fit : problem with fit, parameter errors,  check initial parameters !"
            print 'covariance matrix : ', self.covar
        if self.print_results:
            print '----------------------------------------------------------------------'
            print 'fit results : '
            print '----------------------------------------------------------------------'
            print 'chisquare = ',self.chi2
            print 'red. chisquare = ',self.chi2_red
            print 'parameters: '
            self.show_parameters()

    def set_yval(self, y, y_err = None):
        """
        
        set the array of values to be fitted

        """
        self.y = y
        self.y_err = y_err

    def set_xval(self, x):
        """
        
        set the array of values x

        """
        self.x = x

    def show_parameters(self):
        """

        show the fitted parameters


        """
        for i,p in enumerate(self.parameters):
            print 'parameter ', i, ' : ', p




