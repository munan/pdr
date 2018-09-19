import numpy as np
from scipy.optimize import brentq
from scipy.optimize import newton
import shielding
import const

def fH2_CR(nH, kcr=2.0e-16, Zd=1.):
    kgr=3.0e-17*Zd
    R = kgr * nH / (2. * kcr)
    b = - (1.65*1.5 + 2.*R)
    a = 1.65*0.7
    c = R
    x_H2 = (-b - np.sqrt(b*b - 4*a*c) )/(2.*a)
    return x_H2

def fH2_CR_FUV(nH, GH2, kcr=2.0e-16, Zd=1.):
    """Perhaps we will need G_H2 that's shelf and dust shielded in the
    implementation."""
    kgr=3.0e-17*Zd
    k_FUV = 5.7e-11 * GH2
    R = kgr * nH / (2. * kcr)
    b = - (1.65*1.5 + 2.*R + k_FUV/(2.*kcr))
    a = 1.65*0.7
    c = R
    x_H2 = (-b - np.sqrt(b*b - 4*a*c) )/(2.*a)
    return x_H2

def fHplus_e(x_e, xCplus, nH, T, G_PE, kcr=2.0e-16, Zd=1.):
    small_ = 1e-50
    x_H2 = fH2_CR(nH, kcr=kcr, Zd=Zd)
    x_H = max(1. - 2. * x_H2 - (x_e-xCplus), 0)
    k_Hplus_e = 2.753e-14 * pow( 315614.0 / T, 1.5) * pow( 
                 1.0 + pow( 115188.0 / T, 0.407) , -2.242 )
    cgr_Hplus = [12.25, 8.074e-6, 1.378, 5.087e2, 1.586e-2, 0.4723, 1.102e-5]
    psi_gr = 1.7 * G_PE * np.sqrt(T)/(nH * x_e + small_)
    k_Hplus_gr = 1e-14*cgr_Hplus[0]/(1 + cgr_Hplus[1]*psi_gr**cgr_Hplus[2]*(
             1 + cgr_Hplus[3]*T**cgr_Hplus[4] *
             psi_gr**(-cgr_Hplus[5]-cgr_Hplus[6]*np.log(T))))
    k_Hplus_gr *= Zd
    k_cr_H = kcr * (2.3*x_H2 + 1.5*x_H)
    c = k_cr_H * x_H / nH 
    x_Hplus = c/(k_Hplus_e *  x_e + k_Hplus_gr)
    return min(x_Hplus, 1.)

def fCplus_e(x_e, nH, T, G_PE, G_CI, kcr=2.0e-16, Zd=1., xCstd=1.6e-4):
    small_ = 1e-50
    x_H2 = fH2_CR(nH, kcr=kcr, Zd=Zd)
    k_C_cr = 3.85 * kcr
    k_C_photo = 3.5e-10*G_CI
    k_Cplus_e = Cplus_rec_rate(T)
    psi_gr = 1.7 * G_PE * np.sqrt(T)/(nH * x_e + small_)
    cgr_Cplus = [45.58, 6.089e-3, 1.128, 4.331e2, 4.845e-2, 0.8120, 1.333e-4]
    k_Cplus_gr = 1e-14*cgr_Cplus[0]/(1 + cgr_Cplus[1]*psi_gr**cgr_Cplus[2]*(
             1 + cgr_Cplus[3]*T**cgr_Cplus[4] *
             psi_gr**(-cgr_Cplus[5]-cgr_Cplus[6]*np.log(T))))
    k_Cplus_gr *= Zd 
    k_Cplus_H2 = 3.3e-13 * T**(-1.3) * np.exp(-23./T)
    c = (k_C_cr + k_C_photo) / nH
    al = k_Cplus_e*x_e + k_Cplus_gr + k_Cplus_H2*x_H2 + c
    ar = xCstd * Zd * c
    x_Cplus = ar / al
    return x_Cplus

def fe(nH, T, G_PE, G_CI, kcr=2.0e-16, Zd=1.):
    """Assuming x_Hplus + x_Cplus = x_e, only works for scalar arguments."""
    def fun(x, nH, T, G_PE, G_CI, kcr, Zd):
        x_Cplus = fCplus_e(x, nH, T, G_PE, G_CI, kcr=kcr, Zd=Zd)
        x_Hplus = fHplus_e(x, x_Cplus, nH, T, G_PE, kcr=kcr, Zd=Zd)
        eq = x_Hplus + x_Cplus - x
        return eq
    f1 = fun(1., nH, T, G_PE, G_CI, kcr, Zd)
    if f1 > 0.:
        xeq = 1.
    else:
        xeq = brentq(fun, 0., 1., args=(nH, T, G_PE, G_CI, kcr, Zd),
                     rtol=1e-2, xtol=1e-10)
    return xeq

def fe_iter(nH, T, G_PE, G_CI, kcr=2.0e-16, Zd=1., maxiter=20, x0=0.5):
    def fun(x, nH, T, G_PE, G_CI, kcr, Zd):
        x_Cplus = fCplus_e(x, nH, T, G_PE, G_CI, kcr=kcr, Zd=Zd)
        x_Hplus = fHplus_e(x, x_Cplus, nH, T, G_PE, kcr=kcr, Zd=Zd)
        return x_Hplus + x_Cplus
    xprev = x0
    niter = 0
    rtol = 1e-2
    small_x = 1e-6
    small_f = 1e-20
    a = 0.
    fa = 0.
    b = 1.
    flag = True
    fprev = fun(xprev, nH, T, G_PE, G_CI, kcr, Zd) - xprev
    while True:
        x = fun(xprev, nH, T, G_PE, G_CI, kcr, Zd)
        if abs((x-xprev)/(xprev+small_x)) < rtol or abs(fprev) < small_x:
            break
        if niter > maxiter:
            print "fe_iter: WARNING: niter > maxiter, x={:.2e}, xprev={:.2e}.".format(
                    x, xprev)
            break
        f = fun(x, nH, T, G_PE, G_CI, kcr, Zd) - x
        if abs(f - fprev) < small_f:
            xnew = (x + xprev)/2.
        else:
            xnew = x - f * (x-xprev)/(f-fprev)
            if xnew < a or xnew > b:
                xnew = (a + b)/2.
        if flag:
            fa = fun(a, nH, T, G_PE, G_CI, kcr, Zd)-a
        fnew = fun(xnew, nH, T, G_PE, G_CI, kcr, Zd)-xnew
        if fa*fnew < 0:
            b = xnew
            flag = False
        else:
            a = xnew
            flag = True
        xprev = xnew
        fprev = fnew
        niter = niter + 1
    return x, niter
def fe_dekker(nH, T, G_PE, G_CI, kcr=2.0e-16, Zd=1., maxiter=20, x0=0.5, a=0, b=1.):
    def fun(x):
        x_Cplus = fCplus_e(x, nH, T, G_PE, G_CI, kcr=kcr, Zd=Zd)
        x_Hplus = fHplus_e(x, x_Cplus, nH, T, G_PE, kcr=kcr, Zd=Zd)
        return x_Hplus + x_Cplus - x
    f_small = 1.0e-10
    x_small = 1.0e-10
    rtol = 1.0e-2
    delta = 1e-10
    fa = fun(a)
    fb = fun(b)
    if abs(fa) < f_small:
        return a, 0
    if abs(fb) < f_small:
        return b, 0
    if fa*fb > 0.:
        raise RuntimeError("Signs of f(a) and f(b) must be opposite.")
    s = b
    fs = fb
    niter = 0
    while True:
        if abs(fs) < f_small or abs((b-a)/(b+x_small)) < rtol:
            #print "return: niter={:d}, s={:.2e}".format(niter, s)
            return s, niter
        if niter > maxiter:
            #print "fe_iter: WARNING: niter > maxiter"
            break
        if abs(fa - fb) < f_small:
            s = (a + b) * 0.5
        else:
            s = b - fb * (b - a) / (fb - fa)
        fs = fun(s)
        if fa*fs < 0:
            b = s
            fb = fs
        else:
            a = s
            fa = fs
        niter = niter + 1
    return s, niter

def fe_brent(nH, T, G_PE, G_CI, kcr=2.0e-16, Zd=1., maxiter=20, x0=0.5, a=0., b=1.):
    def fun(x):
        x_Cplus = fCplus_e(x, nH, T, G_PE, G_CI, kcr=kcr, Zd=Zd)
        x_Hplus = fHplus_e(x, x_Cplus, nH, T, G_PE, kcr=kcr, Zd=Zd)
        return x_Hplus + x_Cplus - x
    def swap(x, y):
        return y, x
    f_small = 1.0e-20
    x_small = 1.0e-10
    rtol = 1.0e-2
    delta = 1e-10
    fa = fun(a)
    fb = fun(b)
    if abs(fa) < f_small:
        return a, 0
    if abs(fb) < f_small:
        return b, 0
    if fa*fb > 0.:
        raise RuntimeError("Signs of f(a) and f(b) must be opposite.")
    if abs(fa) < abs(fb):
        a, b = swap(a, b)
        fa, fb = swap(fa, fb)
    c = a
    fc = fa
    s = b
    d = 0.0
    mflag = True
    niter = 0
    fs = fb
    while True:
        if abs(fs) < f_small or abs((b-a)/(b+x_small)) < rtol:
            #print "return: niter={:d}, s={:.2e}".format(niter, s)
            return s, niter
        if niter > maxiter:
            #print "fe_iter: WARNING: niter > maxiter"
            break
        if abs(fa-fc) > f_small and abs(fb-fc) > f_small:
            # use inverse quadratic interopolation
            s =  (   ( a * fb * fc / ((fa - fb) * (fa - fc)) )
                + ( b * fa * fc / ((fb - fa) * (fb - fc)) )
                + ( c * fa * fb / ((fc - fa) * (fc - fb)) )   )
        else:
            # secant method
            s = b - fb * (b - a) / (fb - fa)
        if ( ((s < (3 * a + b) * 0.25) or (s > b)) or (
            mflag and (abs(s-b) >= (abs(b-c) * 0.5)) ) or (
            ~mflag and (abs(s-b) >= (abs(c-d) * 0.5)) ) or (
            mflag and (abs(b-c) < delta) ) or (
            ~mflag and (abs(c-d) < delta) ) ):
            # bisection method
            s = (a+b)*0.5
            mflag = True
        else:
            mflag = False
        fs = fun(s)
        d = c
        c = b
        fc = fb
        if fa * fs < 0:
            b = s
            fb = fs
        else:
            a = s
            fa = fs
        if abs(fa) < abs(fb):
            a, b = swap(a, b)
            fa, fb = swap(fa, fb)
        niter = niter + 1
    return s, niter

def Cplus_rec_rate(T):
    A = 2.995e-9
    B = 0.7849
    T0 =  6.670e-3
    T1 = 1.943e6
    C = 0.1597
    T2 = 4.955e4
    BN = B + C * np.exp(-T2/T);
    term1 = np.sqrt(T/T0);
    term2 = np.sqrt(T/T1);
    alpharr = A / ( term1*pow(1.0+term1, 1.0-BN) * pow(1.0+term2, 1.0+BN) )
    alphadr = pow( T, -3.0/2.0 ) * ( 6.346e-9 * np.exp(-1.217e1/T) +
        9.793e-09 * np.exp(-7.38e1/T) + 1.634e-06 * np.exp(-1.523e+04/T) )
    return (alpharr+alphadr)

def fCO(nH, x_Cplus, x_H2, GCO, kcr=2.0e-16, Zd=1., xCstd=1.6e-4):
    """Use the fitting function in GOW17. Only works for scalar arguments.
       Use x_Cplus x_H2 to limit the CO abundance."""
    kcr16 = kcr/1.0e-16
    term1 = 4.0e3*Zd*kcr16**(-2)
    term1 = np.fmax(term1, 1.)
    ncrit2 = 2.*( term1 )**(GCO**(1./3.)) * (50.*kcr16/Zd**1.4)
    x_CO = nH/ncrit2;
    if type(nH) == type(0.) or type(nH) == type(0):
        if nH >= ncrit2:
            x_CO = 1.
    else:
        indx = (nH >= ncrit2)
        x_CO[indx] = 1.
    x_CO *= xCstd*Zd
    x_CO_max1 = xCstd*Zd - x_Cplus
    x_CO_max2 = xCstd*Zd*x_H2*2.
    x_CO = np.fmin(x_CO, x_CO_max1)
    x_CO = np.fmin(x_CO, x_CO_max2)
    return x_CO

def fCO_exp(nH, GCO, kcr=2.0e-16, Zd=1., xCstd=1.6e-4):
    """Use the expomential smoothing. For solar neighborhood, the linear
    smoothing in f_CO seem to work better"""
    kcr16 = kcr/1.0e-16
    ncrit = ( 4.0e3*Zd*kcr16**(-2) )**(GCO**(1./3.)) * (50.*kcr16/Zd**1.4)
    x_CO = 1. - np.exp(-nH/ncrit)
    x_CO *= xCstd*Zd
    return x_CO

