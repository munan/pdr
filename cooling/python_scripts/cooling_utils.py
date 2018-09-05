import numpy as np
from scipy.optimize import brentq
from scipy.optimize import newton
import shielding

def fH2_CR(nH, kcr=2.0e-16, Zd=1.):
    kgr=3.0e-17*Zd
    R = kgr * nH / (2. * kcr)
    b = - (1.65*1.5 + 2.*R)
    a = 1.65*0.7
    c = R
    x_H2 = (-b - np.sqrt(b*b - 4*a*c) )/(2.*a)
    return x_H2

def fH2_CR_FUV(nH, NH, G0, kcr=2.0e-16, b=3.0e5, xH2_fs=1., Zd=1.):
    """Perhaps we will need G_H2 that's shelf and dust shielded in the
    implementation."""
    kgr=3.0e-17*Zd
    f_dust = shielding.fH2_dust(NH)
    f_self = shielding.fH2_self(NH/2. * xH2_fs, b)
    k_FUV = 5.7e-11 * G0 * f_dust * f_self
    R = kgr * nH / (2. * kcr)
    b = - (1.65*1.5 + 2.*R + k_FUV/(2.*kcr))
    a = 1.65*0.7
    c = R
    x_H2 = (-b - np.sqrt(b*b - 4*a*c) )/(2.*a)
    return x_H2

def fHplus_e(x_e, nH, T, G_PE, kcr=2.0e-16, Zd=1.):
    x_H2 = fH2_CR(nH, kcr=kcr, Zd=Zd)
    x_H = 1. - 2. * x_H2
    k_Hplus_e = 2.753e-14 * pow( 315614.0 / T, 1.5) * pow( 
                 1.0 + pow( 115188.0 / T, 0.407) , -2.242 )
    cgr_Hplus = [12.25, 8.074e-6, 1.378, 5.087e2, 1.586e-2, 0.4723, 1.102e-5]
    psi_gr = 1.7 * G_PE * np.sqrt(T)/(nH * x_e)
    k_Hplus_gr = 1e-14*cgr_Hplus[0]/(1 + cgr_Hplus[1]*psi_gr**cgr_Hplus[2]*(
             1 + cgr_Hplus[3]*T**cgr_Hplus[4] *
             psi_gr**(-cgr_Hplus[5]-cgr_Hplus[6]*np.log(T))))
    k_Hplus_gr *= Zd
    k_cr_H = kcr * (2.3*x_H2 + 1.5*x_H)
    c = k_cr_H * x_H / nH 
    x_Hplus = c/(k_Hplus_e *  x_e + k_Hplus_gr)
    return x_Hplus

def fCplus_e(x_e, nH, T, G_PE, G_CI, kcr=2.0e-16, Zd=1., xCstd=1.6e-4):
    x_H2 = fH2_CR(nH, kcr=kcr, Zd=Zd)
    k_C_cr = 3.85 * kcr
    k_C_photo = 3.5e-10*G_CI
    k_Cplus_e = Cplus_rec_rate(T)
    psi_gr = 1.7 * G_PE * np.sqrt(T)/(nH * x_e)
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
        x_Hplus = fHplus_e(x, nH, T, G_PE, kcr=kcr, Zd=Zd)
        x_Cplus = fCplus_e(x, nH, T, G_PE, G_CI, kcr=kcr, Zd=Zd)
        eq = x_Hplus + x_Cplus - x
        return eq
    f1 = fun(1., nH, T, G_PE, G_CI, kcr, Zd)
    if f1 > 0.:
        xeq = 1.
    else:
        xeq = brentq(fun, 0., 1., args=(nH, T, G_PE, G_CI, kcr, Zd),
                     rtol=1e-2, xtol=1e-10)
    return xeq

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

