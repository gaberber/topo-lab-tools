k= 8.617e-5 #eV/K Boltzmann
T= 55e-3 #temp of 1 mK
import numpy as np
from scipy.optimize import curve_fit

def gap(V, delta, Z, broadness, offset=0):
    def ingap(E, delta):
        return E <= delta

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

    def f0_deriv(expo):
        return ((c/(k*T))*expo)/(expo+1)**2

    def f0(expo):
        return 1/(expo+1)

    def per_e(e, v, delta, Z, broadness):
        expo_f = exp(v)
        expo_f_deriv = exp(e-c*v)
        return (f0_deriv(expo_f_deriv)-f0(expo_f))*(1 + A(e, delta, Z)-B(e, delta, Z))

    def exp(E):
        return np.exp(np.real(E)/(k*T))

    def DIDV(V, delta, Z, broadness):
        return np.sum(np.array([np.real(per_e(e+1j*broadness, V, delta, Z, broadness)) for e in E])  /len(E))

    #Rn = .5*norm(delta)
    #return scale*(1/Rn)*np.array([DIDV(abs(v-offset)) for v in V]) data not well normalized
    norm= 2/DIDV(0, delta, 0, 0)
    return np.array([norm * DIDV(abs(v-offset), delta, Z, broadness) for v in V])

def double_gap(V, delta1=2e-4, delta2=2e-4, Z1=1.0, Z2=1.0, broadness1 =1e-6, broadness2=1e-6, w=.5, offset=0):
    def ingap(E, delta):
        return E <= delta

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

    def f0_deriv(expo):
        return ((c/(k*T))*expo)/(expo+1)**2

    def f0(expo):
        return 1/(expo+1)

    def per_e(e, v, delta, Z, broadness):
        expo_f = exp(v)
        expo_f_deriv = exp(e-c*v)
        return (f0_deriv(expo_f_deriv)-f0(expo_f))*(1 + A(e, delta, Z)-B(e, delta, Z))

    def exp(E):
        return np.exp(np.real(E)/(k*T))

    def DIDV(V, delta, Z, broadness):
        return np.sum(np.array([np.real(per_e(e+1j*broadness, V, delta, Z, broadness)) for e in E])  /len(E))

    #Rn = .5*norm(delta)
    #return scale*(1/Rn)*np.array([DIDV(abs(v-offset)) for v in V]) data not well normalized
    norm1, norm2 = 2/DIDV(0, delta1, 0, 0), 2/DIDV(0, delta2, 0, 0)
    return np.array([w * norm1 *DIDV(abs(v-offset), delta1, Z1, broadness1) +
                     (1-w) * norm2 *DIDV(abs(v-offset), delta2, Z2, broadness2) for v in V])

def fit_gap(Vn, 
           gn, 
           p0=[2e-4,1.0,1e-6,0],
           bounds = (
            (2.4e-4,0,0, -.1),
            (2.7e-4,2,1e-4, .1)
            )):
    """
    Fits a NS linecurve
    
    Arguments:
        Vn: Vbias range in volts
        gn: conductance data in (single) conductance quanta
        p0: ansatz fitting parameters like: [delta1, delta2, Z1, Z2, broadness1, broadness2, w, offset]
        delta is the gap, Z is the barrier parameter, broadness the Dynes parameter, w weight factor of both gaps, offset is zero bias offset
        bounds: tuple of bounds for fitting parameters
    """
    popt, pcov = curve_fit(gap, Vn, gn, p0 = p0, bounds = bounds)
    g_fit = gap(Vn, **popt)
    return g_fit, popt, pcov

def fit_double_gap(Vn, 
                   gn, 
                   p0=[2e-4,2e-4, 1.0,1.0, 1e-6, 1e-6, .5, 0],
                   bounds = (
                    (2.4e-4,2e-4,0,0,0,0,0, -.1),
                    (2.7e-4,2.4e-4,2,2,1e-4,1e-4,1, .1)
                    )):
    """
    Fits a linecurve that shows two seperate gaps
    
    Arguments:
        Vn: Vbias range in volts
        gn: conductance data in (single) conductance quanta
        p0: ansatz fitting parameters like: [delta1, delta2, Z1, Z2, broadness1, broadness2, w, offset]
        delta is the gap, Z is the barrier parameter, broadness the Dynes parameter, w weight factor of both gaps, offset is zero bias offset
        bounds: tuple of bounds for fitting parameters
    """
    popt, pcov = curve_fit(double_gap, Vn, gn, p0 = p0, bounds = bounds)
    g_fit = two_gap(Vn, **popt)
    return g_fit, popt, pcov