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

#     if full_output:
#         retlist = res['x'], res['fun'], res['nfev'], res['njev'], res['status']
#         if retall:
#             retlist += (res['allvecs'], )
#         return retlist
#     else:
#         if retall:
#             return res['x'], res['allvecs']
#         else:
#             return res['x']
        
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
#     grad = self.approx_fprime(x0,epsilon)
    k = 0
    xk = x0

    # Sets the initial step guess to dx ~ 1
    old_fval,grad = self.enerGrad(xk)
    old_old_fval = old_fval + npnorm(grad) / 2

    if retall:
        allvecs = [xk]
    warnflag = 0
    self.pk = -grad
    gnorm = vecnorm(grad, ord=norm)
    methodMin = 'steepest'
    while (gnorm > gtol) and (k < maxiter):# and (abs(old_fval - old_old_fval)>0.01):
#         print 'k,gnorm, energy',k,gnorm,old_fval
        if methodMin == 'conjGrad':       
            deltak = dot(grad, grad)
            stp_k, old_fval, old_old_fval, gradp1 = \
                         self.line_search_wolfe1(xk, epsilon, grad, old_fval,
                                              old_old_fval, c2=0.4, amin=1e-100, amax=1e100)
            print 'step',stp_k

            xk = xk + stp_k * self.pk
            if retall:
                allvecs.append(xk)
            yk = gradp1 - grad
            beta_k = max(0, dot(yk, gradp1) / deltak)
            self.pk = -gradp1 + beta_k * self.pk
            grad = gradp1
        elif methodMin == 'steepest':
            stp_k = 0.02
#             print 'force step to be', stp_k
            xk = xk + stp_k * self.pk
            fnew,grad = self.enerGrad(xk)
            self.plotPos(self.points,self.npoints,'pos','{}_'.format(str(k)))
            print 'xk  ',xk
            print 'grad', grad  
            print
            self.pk = -grad         
        gnorm = vecnorm(grad, ord=norm)
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
#             print("         Function evaluations: %d" % func_calls[0])
#             print("         Gradient evaluations: %d" % grad_calls[0])
    

#     result = OptimizeResult(fun=fval, jac=grad, nfev=func_calls[0],
#                             njev=grad_calls[0], status=warnflag,
#                             success=(warnflag == 0), message=msg, x=xk,
#                             nit=k)
#     if retall:
#         result['allvecs'] = allvecs
#     return result

# def approx_fprime(self, xk, epsilon):
#     """
#     See ``approx_fprime``.  An optional initial function value arg is added.
#     """
#     f0 = self.enerGrad(xk)
#     grad = zeros((len(xk),), float)
#     ei = zeros((len(xk),), float)
#     for k in range(len(xk)):
#         ei[k] = 1.0
#         d = epsilon * ei #d epsilon only at k, 0 otherwise
#         grad[k] = (self.enerGrad(xk + d) - f0) / d[k]
#         ei[k] = 0.0
#     return grad


def vecnorm(x, ord=2):
    if ord == Inf:
        return amax(abs(x))
    elif ord == -Inf:
        return amin(abs(x))
    else:
        return sum(abs(x)**ord, axis=0)**(1.0 / ord)
    

    
def line_search_wolfe1(self,xk, epsilon, grad, old_fval, old_old_fval,
                       c1=1e-4, c2=0.9, amax=50, amin=1e-8,
                       xtol=1e-14):
    """
    See https://github.com/scipy/scipy/blob/v0.13.0/scipy/optimize/linesearch.py#L187
    ***Modified***

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
    grad : array_like, optional
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
    derf0 = dot(grad, self.pk)
    stp, fval, old_fval, grad  = self.scalar_search_wolfe1(xk,epsilon,old_fval, 
                                                    old_old_fval, derf0, grad)
    return stp, fval, old_fval, grad



def scalar_search_wolfe1(self, xk,epsilon,f0, old_f0, derf0, grad,
                         c1=1e-4, c2=0.9,
                         amax=50, amin=1e-8, xtol=1e-14):
    """
    See https://github.com/scipy/scipy/blob/v0.13.0/scipy/optimize/linesearch.py#L187
    ***Modified***
    
    Scalar function search for stp that satisfies strong Wolfe conditions
    stp > 0 is assumed to be a descent direction.
    Parameters
    ----------
    f : callable f(stp)
        Function at point `stp`
    derf : callable df(stp)
        Derivative `d f(stp)/ds`. Returns a scalar.
    f0 : float, optional
        Value of `f` at 0
    old_f0 : float, optional
        Value of `f` at the previous point
    derf0 : float, optional
        Value `derf` at 0
    amax : float, optional
        Maximum step size
    c1, c2 : float, optional
        Wolfe parameters
    Returns
    -------
    stp : float
        Step size, or None if no suitable step was found
    f : float
        Value of `f` at the new point `stp`
    f0 : float
        Value of `f` at `stp=0`
    Notes
    -----
    Uses routine DCSRCH from MINPACK.
    """       


    stp1 = min(1.0, 1.01*2*(f0 - old_f0)/derf0)
    if stp1 < 0:
        stp1 = 1.0
    f1 = f0
    derf1 = derf0
#     isave = zeros((2,), intc)
#     dsave = zeros((13,), float)
    isave = zeros((3,), intc) #bch increase these by 1 over original because dcsrch starts counting at 1
    dsave = zeros((14,), float)
    task = 'START'

    maxiter = 30
    for i in xrange(maxiter):
        stp, task, f1, derf1, isave, dsave  = dcsrch(stp1, f1, derf1,
                                                   c1, c2, xtol, task,
                                                   amin, amax, isave, dsave)
#def dcsrch(stp, f, g, ftol, gtol, xtol, task, stpmin, stpmax, isave, dsave):
#return stp, task, isave, dsave
#         stpA = stp
#         maxMove = self.rmax/10
#         move = npnorm(stp*self.pk)
#         if move > maxMove:
#             stp  = stp*maxMove/move
#             print'step too large',stpA, 'now',stp
#             print
        
        if task[:2] == 'FG':
#             stp1 = stp
            f1, grad = self.enerGrad(xk + stp*self.pk)
            derf1 = dot(grad,self.pk)
#             def f(self,xk,s):
#     self.fc[0] += 1
#     ener = self.enerGrad(xk + s*self.pk)
#     return ener
# 
# def derf(self,xk,epsilon,s):
#     self.gval[0] = self.approx_fprime(xk + s*self.pk,epsilon)
#     self.gc[0] += 1
# #     else:
# #         fc[0] += len(xk) + 1
#     return dot(self.gval[0],self.pk)
#             
            
            
        else:
            break
    else: #no break was found
        # maxiter reached, the line search did not converge
        stp = None

    if task[:5] == 'ERROR' or task[:4] == 'WARN':
        stp = None  # failed

    return stp, f1, f0, grad


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
