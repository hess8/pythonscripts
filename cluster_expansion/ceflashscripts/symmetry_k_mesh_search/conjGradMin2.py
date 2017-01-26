import warnings
import sys
# import numpy
from scipy._lib.six import callable
from numpy import (atleast_1d, eye, mgrid, argmin, zeros, shape, squeeze,
                   vectorize, asarray, sqrt, Inf, asfarray, isinf, intc,
                   amax, amin, dot)
from numpy.linalg import det, inv, eig
from numpy.linalg import norm as npnorm
from scipy._lib._util import getargspec_no_self as _getargspec

from dcsrch import dcsrch

_status_message = {'success': 'Optimization terminated successfully.',
                   'maxfev': 'Maximum number of function evaluations has '
                              'been exceeded.',
                   'maxiter': 'Maximum number of iterations has been '
                              'exceeded.',
                   'pr_loss': 'Desired error not necessarily achieved due '
                              'to precision loss.'}
# _epsilon = sqrt(finfo(float).eps)
# _epsilon = rc

def fmin_cg(self, x0, epsilon, fprime=None, args=(), gtol=1e-5, norm=Inf,
            maxiter=None, full_output=0, disp=1, retall=0, callback=None):

    opts = {'gtol': gtol,
            'norm': norm,
            'disp': disp,
            'maxiter': maxiter,
            'return_all': retall}

    res = self.minimize_cg(x0, epsilon, args, callback=callback, **opts)

    if full_output:
        retlist = res['x'], res['fun'], res['nfev'], res['njev'], res['status']
        if retall:
            retlist += (res['allvecs'], )
        return retlist
    else:
        if retall:
            return res['x'], res['allvecs']
        else:
            return res['x']
        
def minimize_cg(self,x0, epsilon, args=(), jac=None, callback=None,
                 gtol=1e-5, norm=Inf, maxiter=None,
                 disp=False, return_all=False):
    """
    Minimization of scalar function of one or more variables using the
    conjugate gradient algorithm.
    Options
    -------
    disp : bool
        Set to True to print convergence messages.
    maxiter : int
        Maximum number of iterations to perform.
    gtol : float
        Gradient norm must be less than `gtol` before successful
        termination.
    norm : float
        Order of norm (Inf is max, -Inf is min).
    eps : float or ndarray
        If `jac` is approximated, use this value for the step size.
    """

    retall = return_all
    x0 = asarray(x0).flatten()
    if maxiter is None:
        maxiter = len(x0) * 200
#     func_calls, f = wrap_function(f, args)
#     grad_calls, myfprime = approx_fprime, (f, epsilon))
    gfk = self.approx_fprime(x0,epsilon)
    k = 0
    xk = x0

    # Sets the initial step guess to dx ~ 1
    old_fval = self.energy(xk)
    old_old_fval = old_fval + npnorm(gfk) / 2

    if retall:
        allvecs = [xk]
    warnflag = 0
    self.pk = -gfk
    gnorm = vecnorm(gfk, ord=norm)
    while (gnorm > gtol) and (k < maxiter):
        print 'k',k
        deltak = dot(gfk, gfk)

#         try:
        alpha_k, fc, gc, old_fval, old_old_fval, gfkp1 = \
                     self.line_search_wolfe1(xk, epsilon, gfk, old_fval,
                                          old_old_fval, c2=0.4, amin=1e-100, amax=1e100)
        sys.exit('stop')
        xk = xk + alpha_k * self.pk
        if retall:
            allvecs.append(xk)
        if gfkp1 is None:
            gfkp1 = myfprime(xk)
        yk = gfkp1 - gfk
        beta_k = max(0, dot(yk, gfkp1) / deltak)
        self.pk = -gfkp1 + beta_k * self.pk
        gfk = gfkp1
        gnorm = vecnorm(gfk, ord=norm)
        if callback is not None:
            callback(xk)
        k += 1

    fval = old_fval
    if warnflag == 2:
        msg = _status_message['pr_loss']
        if disp:
            print("Warning: " + msg)
            print("         Current function value: %f" % fval)
            print("         Iterations: %d" % k)
            print("         Function evaluations: %d" % func_calls[0])
            print("         Gradient evaluations: %d" % grad_calls[0])

    elif k >= maxiter:
        warnflag = 1
        msg = _status_message['maxiter']
        if disp:
            print("Warning: " + msg)
            print("         Current function value: %f" % fval)
            print("         Iterations: %d" % k)
            print("         Function evaluations: %d" % func_calls[0])
            print("         Gradient evaluations: %d" % grad_calls[0])
    else:
        msg = _status_message['success']
        if disp:
            print(msg)
            print("         Current function value: %f" % fval)
            print("         Iterations: %d" % k)
            print("         Function evaluations: %d" % func_calls[0])
            print("         Gradient evaluations: %d" % grad_calls[0])

    result = OptimizeResult(fun=fval, jac=gfk, nfev=func_calls[0],
                            njev=grad_calls[0], status=warnflag,
                            success=(warnflag == 0), message=msg, x=xk,
                            nit=k)
    if retall:
        result['allvecs'] = allvecs
    return result

def approx_fprime(self, xk, epsilon):
    """
    See ``approx_fprime``.  An optional initial function value arg is added.
    """
    f0 = self.energy(xk)
    grad = zeros((len(xk),), float)
    ei = zeros((len(xk),), float)
    for k in range(len(xk)):
        ei[k] = 1.0
        d = epsilon * ei #d is changing each time
        grad[k] = (self.energy(xk + d) - f0) / d[k]
        ei[k] = 0.0
    return grad


def vecnorm(x, ord=2):
    if ord == Inf:
        return amax(abs(x))
    elif ord == -Inf:
        return amin(abs(x))
    else:
        return sum(abs(x)**ord, axis=0)**(1.0 / ord)
    
def phi(self,xk,s):
    self.fc[0] += 1
    ener = self.energy(xk + s*self.pk)
    return ener

def derphi(self,xk,epsilon,s):
    self.gval[0] = self.approx_fprime(xk + s*self.pk,epsilon)
    self.gc[0] += 1
#     else:
#         fc[0] += len(xk) + 1
    return dot(self.gval[0],self.pk)
    
def line_search_wolfe1(self,xk, epsilon,gfk, old_fval, old_old_fval,
                       c1=1e-4, c2=0.9, amax=50, amin=1e-8,
                       xtol=1e-14):
    """
    As `scalar_search_wolfe1` but do a line search to direction `pk`
    Parameters
    ----------
    f : callable
        Function `f(x)`
    fprime : callable
        Gradient of `f`
    xk : array_like
        Current point
    pk : array_like
        Search direction
    gfk : array_like, optional
        Gradient of `f` at point `xk`
    old_fval : float, optional
        Value of `f` at point `xk`
    old_old_fval : float, optional
        Value of `f` at point preceding `xk`
    The rest of the parameters are the same as for `scalar_search_wolfe1`.
    Returns
    -------
    stp, f_count, g_count, fval, old_fval
        As in `line_search_wolfe1`
    gval : array
        Gradient of `f` at the final point
    """

#     if isinstance(fprime, tuple):
#         eps = fprime[1]
#         fprime = fprime[0]
#         newargs = (f, eps) + args
#         gradient = False
#     else:
#     newargs = args
    self.gval = [gfk]
    self.gc = [0]
    self.fc = [0]
    derphi0 = dot(gfk, self.pk)
    stp, fval, old_fval = self.scalar_search_wolfe1(xk,epsilon,old_fval, old_old_fval, derphi0)
    return stp, self.fc[0], self.gc[0], fval, old_fval, self.gval[0]



def scalar_search_wolfe1(self, xk,epsilon,phi0, old_phi0, derphi0,
                         c1=1e-4, c2=0.9,
                         amax=50, amin=1e-8, xtol=1e-14):
    """
    Scalar function search for alpha that satisfies strong Wolfe conditions
    alpha > 0 is assumed to be a descent direction.
    Parameters
    ----------
    phi : callable phi(alpha)
        Function at point `alpha`
    derphi : callable dphi(alpha)
        Derivative `d phi(alpha)/ds`. Returns a scalar.
    phi0 : float, optional
        Value of `f` at 0
    old_phi0 : float, optional
        Value of `f` at the previous point
    derphi0 : float, optional
        Value `derphi` at 0
    amax : float, optional
        Maximum step size
    c1, c2 : float, optional
        Wolfe parameters
    Returns
    -------
    alpha : float
        Step size, or None if no suitable step was found
    phi : float
        Value of `phi` at the new point `alpha`
    phi0 : float
        Value of `phi` at `alpha=0`
    Notes
    -----
    Uses routine DCSRCH from MINPACK.
    """       


    alpha1 = min(1.0, 1.01*2*(phi0 - old_phi0)/derphi0)
    if alpha1 < 0:
        alpha1 = 1.0
    phi1 = phi0
    derphi1 = derphi0
#     isave = zeros((2,), intc)
#     dsave = zeros((13,), float)
    isave = zeros((3,), intc) #bch increase these by 1 over original because dcsrch starts counting at 1
    dsave = zeros((14,), float)
    task = 'START'

    maxiter = 30
    for i in xrange(maxiter):
        stp, task, phi1, derphi1, isave, dsave  = dcsrch(alpha1, phi1, derphi1,
                                                   c1, c2, xtol, task,
                                                   amin, amax, isave, dsave)
#def dcsrch(stp, f, g, ftol, gtol, xtol, task, stpmin, stpmax, isave, dsave):
#return stp, task, isave, dsave
        
        if task[:2] == 'FG':
#             alpha1 = stp
            phi1 = self.phi(xk,stp)
            derphi1 = self.derphi(xk,epsilon,stp)
        else:
            break
    else:
        # maxiter reached, the line search did not converge
        stp = None

    if task[:5] == 'ERROR' or task[:4] == 'WARN':
        stp = None  # failed

    return stp, phi1, phi0


class OptimizeResult(dict):
    """ Represents the optimization result.
    Attributes
    ----------
    x : ndarray
        The solution of the optimization.
    success : bool
        Whether or not the optimizer exited successfully.
    status : int
        Termination status of the optimizer. Its value depends on the
        underlying solver. Refer to `message` for details.
    message : str
        Description of the cause of the termination.
    fun, jac, hess: ndarray
        Values of objective function, its Jacobian and its Hessian (if
        available). The Hessians may be approximations, see the documentation
        of the function in question.
    hess_inv : object
        Inverse of the objective function's Hessian; may be an approximation.
        Not available for all solvers. The type of this attribute may be
        either np.ndarray or scipy.sparse.linalg.LinearOperator.
    nfev, njev, nhev : int
        Number of evaluations of the objective functions and of its
        Jacobian and Hessian.
    nit : int
        Number of iterations performed by the optimizer.
    maxcv : float
        The maximum constraint violation.
    Notes
    -----
    There may be additional attributes not listed above depending of the
    specific solver. Since this class is essentially a subclass of dict
    with attribute accessors, one can see which attributes are available
    using the `keys()` method.
    """
    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError:
            raise AttributeError(name)

    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

    def __repr__(self):
        if self.keys():
            m = max(map(len, list(self.keys()))) + 1
            return '\n'.join([k.rjust(m) + ': ' + repr(v)
                              for k, v in sorted(self.items())])
        else:
            return self.__class__.__name__ + "()"

    def __dir__(self):
        return list(self.keys())
