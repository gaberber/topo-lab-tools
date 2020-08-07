k= 8.617e-5 #eV/K Boltzmann
T= 30e-3 #temp of 1 mK
c=1
import numpy as np
from scipy.optimize import curve_fit

def gap(V, delta, Z, broadness, offset=0):
    E = np.arange(1.5*min(V), 1.5*max(V), 1.5*k*T)

    def u2(E, delta): #u0 squared, BCS coefficient
        return .5 * (1+np.sqrt((E**2 - delta**2)/E**2))

    def v2(E, delta): # v0 squared, BCS coefficient
        return (1-u2(E, delta))

    def gamma2(E, delta, Z): # gamma squared
        return (u2(E, delta) + (Z**2)*(u2(E, delta) - v2(E, delta)))**2

    def A(E, delta, Z):  #Andreev reflection
        if E <= delta:
            return delta**2/(E**2 + (delta**2 - E**2)*(1+2*Z**2)**2)
        if E > delta:
            return (u2(E, delta)*v2(E, delta))/gamma2(E, delta, Z)

    def B(E, delta, Z): #Ordinary reflection
        if E < delta:
            return 1 - delta**2/(E**2 + (delta**2 - E**2)*(1+2*Z**2)**2)
        if E >= delta:
            return ((u2(E, delta)-v2(E, delta))**2*Z**2*(1+Z**2))/gamma2(E, delta, Z)

    def f0_deriv(e, v):
        A = (e-c*v)/(k*T)
        #print(A)
        #A = np.sign(A)* min(20, abs(A))
        #print(type(A))
        A[abs(A)>20] = 20
        _f0_deriv = (c/(2*k*T))/(1+np.cosh(A))
        return _f0_deriv

    def per_e(e, v, delta, Z, broadness):
        _f0_deriv = f0_deriv(e, v)
        return _f0_deriv*(1 + A(e+1j*broadness, delta, Z)-B(e+1j*broadness, delta, Z))

    def DIDV(V, delta, Z, broadness):
        return np.sum(np.array([abs(per_e(e, V, delta, Z, broadness)) for e in E])/len(E), axis=0)

    norm= 2/DIDV(np.array([[0]]), delta, 0, 0)
    return (norm * DIDV(abs(V-offset), delta, Z, broadness)).reshape(len(V))
    #return np.array([norm * DIDV(abs(v-offset), delta, Z, broadness) for v in V])

def double_gap(V, delta1=2e-4, delta2=2e-4, Z1=1.0, Z2=1.0, broadness1 =1e-6, broadness2=1e-6, w=.5, offset=0):
    E = np.arange(1.5*min(V), 1.5*max(V), 1.5*k*T)
    def u2(E, delta): #u0 squared, BCS coefficient
        return .5 * (1+np.sqrt((E**2 - delta**2)/E**2))

    def v2(E, delta): # v0 squared, BCS coefficient
        return (1-u2(E, delta))

    def gamma2(E, delta, Z): # gamma squared
        return (u2(E, delta) + (Z**2)*(u2(E, delta) - v2(E, delta)))**2

    def A(E, delta, Z):
        if E <= delta:
            return delta**2/(E**2 + (delta**2 - E**2)*(1+2*Z**2)**2)
        if E > delta:
            return (u2(E, delta)*v2(E, delta))/gamma2(E, delta, Z)

    def B(E, delta, Z):
        if E < delta:
            return 1 - delta**2/(E**2 + (delta**2 - E**2)*(1+2*Z**2)**2)
        if E >= delta:
            return ((u2(E, delta)-v2(E, delta))**2*Z**2*(1+Z**2))/gamma2(E, delta, Z)

    def f0_deriv(e, v):
        A = (e-c*v)/(k*T)
        A = np.sign(A)* min(20, abs(A))
        _f0_deriv = (c/(2*k*T))/(1+np.cosh(A))
        return _f0_deriv

    def per_e(e, v, delta, Z, broadness):
        _f0_deriv = f0_deriv(e, v)
        return _f0_deriv*(1 + A(e+1j*broadness, delta, Z)-B(e+1j*broadness, delta, Z))

    def DIDV(V, delta, Z, broadness):
        return np.sum(np.array([abs(per_e(e, V, delta, Z, broadness)) for e in E])/len(E))

    norm1, norm2 = 2/DIDV(np.array([[0]]), delta1, 0, 0), 2/DIDV(np.array([[0]]), delta2, 0, 0)
    return (w * norm1 *DIDV(abs(V-offset), delta1, Z1, broadness1) +
                     (1-w) * norm2 *DIDV(abs(V-offset), delta2, Z2, broadness2) ).reshape(len(V))

def fit_gap(Vn, 
           gn, 
        delta=2.25e-4, 
            Z=2, 
            broadness=1e-6, 
            offset=0,
           bounds = (
            (1e-5,0,0, -.1),
            (2.7e-3,4,1e-4, .1)
            )):
    """
    Fits a NS linecurve
    
    Arguments:
        Vn: Vbias range in volts
        gn: conductance data in (single) conductance quanta
        ansatz fitting parameters like: delta, Z, broadening, offset
        delta is the gap, Z is the barrier parameter, broadness the Dynes parameter 
        bounds: tuple of bounds for fitting parameters
    Returns:
        g_fit
        popt
        pcov
    """
    p0=delta, Z, broadness, offset
    popt, pcov = curve_fit(gap, Vn, gn, p0 = p0, bounds = bounds)
    delta, Z, broadness, offset = popt
    g_fit = gap(Vn, delta, Z, broadness, offset)
    return g_fit, popt, pcov

def fit_double_gap(Vn, 
                   gn, 
                   p0=[2e-4,2e-4, 1.0,1.0, 1e-6, 1e-6, .5, 0],
                   bounds = (
                    (1e-3, 1e-3 ,0 , 0, 0, 0, 0, -.1),
                    (2.7e-3, 2.7e-3, 4, 4, 1e-4, 1e-4, 1, .1)
                    )):
    """
    Fits a linecurve that shows two seperate gaps
    
    Arguments:
        Vn: Vbias range in volts
        gn: conductance data in (single) conductance quanta
        ansatz fitting parameters like: delta1, delta2, Z1, Z2, broadness1, broadness2, w, offset
        delta is the gap, Z is the barrier parameter, broadness the Dynes parameter, w weight factor of both gaps, offset is zero bias offset
        bounds: tuple of bounds for fitting parameters
    Returns:
        g_fit
        popt
        pcov
    """
    popt, pcov = curve_fit(double_gap, Vn, gn, p0 = p0, bounds = bounds)
    g_fit = two_gap(Vn, **popt)
    return g_fit, popt, pcov