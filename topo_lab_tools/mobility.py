"""Script for fitting pinch-off curves to extract _Mobility, series resistance and threshold voltage.
Currents are in amps. Vwf is the work function difference between gate metal and semiconductor.
Capacitance can be included from the lookup table or by using the analytical function "_Capacitance_cyl_plane", 
make sure to use the right correction factor in the function (e.g. 1.44 for SiO2 back-gate data)
Check material below (effective mass for InSb/InAs)!
Cut the data right below the threshold voltage to get an accurate estimate for the _Mobility.
For the fundamental constants the file "constants.json" is required.

Written by Sebastian Heedt"""

import numpy as np
import pandas as pd
import math
pi = math.pi
import scipy, scipy.stats, scipy.integrate
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import scipy.constants
import json
import codecs
mu_B = scipy.constants.physical_constants['Bohr magneton'][0]
k_B = scipy.constants.Boltzmann
m_e =  scipy.constants.m_e
e =  scipy.constants.e
h =  scipy.constants.h
hbar =  scipy.constants.hbar
R0 = scipy.constants.physical_constants['inverse of conductance quantum'][0]
Phi0 = scipy.constants.physical_constants['mag. flux quantum'][0]
eps_0 = scipy.constants.epsilon_0

# Fermi velocity:
def _vF(n, dn):
    # n electron concentration in cm^{-3}
    n = n * 1e6
    dn = dn * 1e6
    v_F = hbar / meff * (3 * pi ** 2 * n) ** (1 / 3)
    dv_F = dn * hbar / meff * pi ** 2 * (3 * pi ** 2 * n) ** (-2 / 3)
    return v_F, dv_F


# Fermi energy in meV:
def _E_F(n):
    # n electron concentration in cm^{-3}
    EF = float(format(0.5 * meff * _vF(n, dn=0)[0] ** 2 / e * 1e3, '.3f'))
    EF = EF / 1e3 * e
    return EF


# diffusion constant in 3D:
def _D3(n, mu):
    # n electron concentration in cm^{-3}
    # mu electron _Mobility in cm^2/Vs
    D = 1 / 3 * _vF(n, 0)[0] ** 2 * mu * 1e-4 * meff / e
    return D


# mean free path:
def _l_e(n, dn, mu, dmu):
    # n electron concentration in cm^{-3}
    # mu electron _Mobility in cm^2/Vs
    le = _vF(n, dn)[0] * mu * 1e-4 * meff / e
    dle = np.sqrt((_vF(n, dn)[1] * mu * 1e-4 * meff / e) ** 2 + (_vF(n, dn)[0] * dmu * 1e-4 * meff / e) ** 2)
    return le, dle


def _Capacitance_cyl_plane(LG, d_NW, t_ox, eps_r, hex):
    # LG is the effectie gate length
    # d_NW is the top-view nanowire width (i.e. vertex-to-vertex for the hexagonal cross-section)

    # corrections to the cylinder-on-plane model:
    #                TG (LLO)   BG     BG w/o LLO
    # undoped NWs     15.2       4.4    1.44 [2.35 (w/o SS)]
    # doped NWs       21.9       6.0    ~2.0
    # LaLuO3: eps_r0=26.9
    # Al2O3: eps_r0=8.4
    # SiO2: eps_r0=3.9
    correction = 1.44 / eps_r  # 15.2/26.9

    if hex == True:
        A = 3 * np.sqrt(3) * d_NW ** 2 / 8  # area of the hexagonal cross-section
        d = d_NW * np.sqrt(3 * np.sqrt(3) / 2 / pi)  # equivalent circular diameter
        C = float(format(2 * pi * eps_r * eps_0 * LG / np.arccosh((t_ox + d / 2) / (d / 2)) * correction * 1e15, '.4f'))
        return C / 1e15
    else:
        C = float(
            format(2 * pi * eps_r * eps_0 * LG / np.arccosh((t_ox + d_NW / 2) / (d_NW / 2)) * correction * 1e15, '.4f'))
        return C / 1e15


def _Capacitance_plate(LG, d_NW, t_ox, eps_r):
    # LG is the effectie gate length
    # d_NW is the top-view quantum wire width
    A = LG * d_NW
    C = float(format(eps_r * eps_0 * A / t_ox * 1e15, '.4f'))
    return C / 1e15


def _Mobility(g, L, LG, C, Vsd):
    # g is the transconductance dI/dVsd in S
    # LG is the effective gate length
    mu = float(format(g * L * LG / C / Vsd * 1e4, '.1f'))
    return mu * 1e-4


def _Mobility_wErrors(g, dg, L, dL, LG, dLG, C, Vsd):
    # g is the transconductance dI/dVsd in S
    # LG is the effective gate length
    dC = 0.1 * C
    mu = float(format(g * L * LG / C / Vsd * 1e4, '.1f'))
    dmu = float(format(np.sqrt((L * LG / (C * Vsd)) ** 2 * dg ** 2 + (g * LG / (C * Vsd)) ** 2 * dL ** 2 + (
    g * L / (C * Vsd)) ** 2 * dLG ** 2 + \
                               (g * L * LG / (C ** 2 * Vsd)) ** 2 * dC ** 2) * 1e4, '.1f'))
    return mu * 1e-4, dmu * 1e-4


def _Density(Vth, Vwf, LG, d_NW, C, hex):
    # Vth is the threshold voltage
    # Vwf is the gate metal work function minus the semiconductor electron affinity
    # d_NW is the top-view nanowire width (i.e. vertex-to-vertex for the hexagonal cross-section)

    if hex == True:
        A = 3 * np.sqrt(3) * d_NW ** 2 / 8  # area of the hexagonal cross-section
        d = d_NW * np.sqrt(3 * np.sqrt(3) / 2 / pi)  # equivalent circular diameter
        n = float(format(C * (Vwf - Vth) / e / LG / A / 1e17 / 1e6, '.2f'))
        return n * 1e17 * 1e6
    else:
        A = pi * (d_NW / 2) ** 2
        n = float(format(C * (Vwf - Vth) / e / LG / A / 1e17 / 1e6, '.2f'))
        return n * 1e17 * 1e6


def _Density_wErrors(Vth, dVth, Vwf, LG, dLG, d_NW, dd_NW, C, hex):
    # Vth is the threshold voltage
    # Vwf is the gate metal work function minus the semiconductor electron affinity
    # d_NW is the top-view nanowire width (i.e. vertex-to-vertex for the hexagonal cross-section)
    dC = 0.2 * C

    if hex == True:
        A = 3 * np.sqrt(3) * d_NW ** 2 / 8  # area of the hexagonal cross-section
        d = d_NW * np.sqrt(3 * np.sqrt(3) / 2 / pi)  # equivalent circular diameter
        dd = dd_NW * np.sqrt(3 * np.sqrt(3) / 2 / pi)
        n = float(format(C * (Vwf - Vth) / e / LG / A / 1e17 / 1e6, '.2f'))
        dn = float(format(
            np.sqrt((C / (e * LG * A)) ** 2 * dVth ** 2 + (C * abs(Vwf - Vth) / (e * LG ** 2 * A)) ** 2 * dLG ** 2 + \
                    (3 / 4 * np.sqrt(3) * d * C * abs(Vwf - Vth) / (e * LG * A ** 2)) ** 2 * dd ** 2 + \
                    (abs(Vwf - Vth) / (e * LG * A)) ** 2 * dC ** 2) / 1e17 / 1e6, '.2f'))
        return n * 1e17 * 1e6, dn * 1e17 * 1e6
    else:
        A = pi * (d_NW / 2) ** 2
        d = d_NW
        dd = dd_NW
        n = float(format(C * (Vwf - Vth) / e / LG / A / 1e17 / 1e6, '.2f'))
        dn = float(format(
            np.sqrt((C / (e * LG * A)) ** 2 * dVth ** 2 + (C * abs(Vwf - Vth) / (e * LG ** 2 * A)) ** 2 * dLG ** 2 + \
                    (3 / 4 * np.sqrt(3) * d * C * abs(Vwf - Vth) / (e * LG * A ** 2)) ** 2 * dd ** 2 + \
                    (abs(Vwf - Vth) / (e * LG * A)) ** 2 * dC ** 2) / 1e17 / 1e6, '.2f'))
        return n * 1e17 * 1e6, dn * 1e17 * 1e6


def _Density_drift(mu, rho):
    # Vwf is the gate metal work function minus the semiconductor electron affinity
    # d_NW is the top-view nanowire width (i.e. vertex-to-vertex for the hexagonal cross-section)

    n = float(format(1 / rho / e / (mu * 1e-4) / 1e17 / 1e6, '.2f'))
    return n * 1e17 * 1e6


def _IV_FET(Vgate, Vsd, L, LG, C, mu, Vth, Rs):
    # print(mu, Vth, Rs)
    # I = Vsd/(L*LG/(mu*C*(Vg-Vth)+Rs)
    return Vsd / (L * LG / (mu * C * abs(Vgate - Vth)) + Rs)

def fit_mobility(Vgate, I, Vbias, L, d_NW, NW_mat, ox_mat, t_ox=285e-9, LG=1.9e-6, mu=3e-1, Vth=-.2, Rs=7272):
    """Fits the mobility of a single nanowire IV curve
    Required Arguments:
        Vgate: the back gate range (volts)
        I: the measured current range (amperes)
        Vbias: the source-drain bias voltage in (volts)
        L: the nanowire length (meters)
        d_NW: the nanowire diameter (meters)
        NW_mat: material of nanowire
        ox_mat: oxide material choose from ("Al2O3", "HfO2", "Si3N4", "SiO2", "ZrO2", "air")  
        
    Optional Arguments:
        mu: guess for the mobility incm^2/(v*s) 
        Vth: guess for theshold voltage, i.e. where the device opens up
        Rs: contact resistance of lead with nanowire
        t_ox: oxide thickness
    
    returns:
        popt: the found fit parameters like [mu, Vth, Rs]
        I: the current according to the fit"""
    global meff
    eps_r_dc = {"Al2O3": 9.0, "HfO2": 25.0, "Si3N4": 7.0, "SiO2": 3.9, "ZrO2": 25.0, "air": 1.0}
    meff_dc = {"InSb": 0.014 * m_e, "InAs": 0.026 * m_e}
    
    meff = meff_dc[NW_mat]
    eps_r = eps_r_dc[ox_mat] 
    Vwf = 0.0
    guess = [mu, Vth, Rs]
    C = _Capacitance_cyl_plane(LG, d_NW, t_ox, eps_r, hex)
    popt, pcov = curve_fit(lambda Vgate, mu, Vth, Rs: _IV_FET(Vgate=Vgate, Vsd=Vbias, L=L, LG=LG, \
                                             C=C, mu=mu, Vth=Vth, Rs=Rs), Vgate, I, \
                                             p0=guess, xtol = 0.00000001, maxfev = 200)
    mu, Vth, Rs = popt
    I = _IV_FET(Vgate=Vgate, Vsd=Vbias, L=L, LG=LG, C=C, mu=mu, Vth=Vth, Rs=Rs)
    return popt, I

def IV_FET(Vgate, Vbias, L, d_NW, NW_mat, ox_mat, t_ox=285e-9, LG=1.9e-6, mu=3e-1, Vth=-.2, Rs=7272):
    """Fits the mobility of a single nanowire IV curve
    Required Arguments:
        Vgate: the back gate range (volts)
        Vbias: the source-drain bias voltage in (volts)
        L: the nanowire length (meters)
        d_NW: the nanowire diameter (meters)
        NW_mat: material of nanowire
        ox_mat: oxide material choose from ("Al2O3", "HfO2", "Si3N4", "SiO2", "ZrO2", "air")
        
        
        
    Optional Arguments:
        mu: guess for the mobility (m^2/(v*s) )
        Vth: guess for theshold voltage, i.e. where the device opens up
        Rs: contact resistance of lead with nanowire
        t_ox: oxide thickness
    returns:
        I: the current corresponding to the IV_FET"""
    global meff
    eps_r_dc = {"Al2O3": 9.0, "HfO2": 25.0, "Si3N4": 7.0, "SiO2": 3.9, "ZrO2": 25.0, "air": 1.0}
    meff_dc = {"InSb": 0.014 * m_e, "InAs": 0.026 * m_e}
    
    meff = meff_dc[NW_mat]
    eps_r = eps_r_dc[ox_mat] 
    Vwf = 0.0
    guess = [mu, Vth, Rs]
    C = _Capacitance_cyl_plane(LG, d_NW, t_ox, eps_r, hex)
    I = _IV_FET(Vgate=Vgate, Vsd=Vbias, L=L, LG=LG, C=C, mu=mu, Vth=Vth, Rs=Rs)
    return I