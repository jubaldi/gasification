#!/usr/bin/env python

"""
This script defines functions to calculate energy parameters required for
non-isothermic equilibrium simulation of gasification 
processes. It uses some predefined functions from Cantera package.

@author = Rodolfo Rodrigues
@contact = rodolfo.rodrigues@ufsm.br
@data = May, 2022

"""

#==============================================================================
# import libraries/files
#==============================================================================
import pp as ppold
import feedstock as fsold
import pp2 as pp
import feedstock2 as fs
import cantera as ct
import numpy as np
import scipy.optimize as opt

#==============================================================================
# predefine parameters
#==============================================================================
R = ct.gas_constant  # 8314.4621 Pa*m^3/K/kmol
Tn = 273.15  # K
Pn = ct.one_atm  # 101315 Pa
zero = np.zeros(1)
one = np.ones(1)

#==============================================================================
# special functions
#==============================================================================

def get_h_cp(mix, value='h,cp', duty=0.0):
    '''
    Given a mixture, return either the enthalpy (h) or the specific heat capacity (cp).
    TODO: Add duty term to enthalpy calculation

    Parameters
    ----------
    mix : Cantera 'Mixture' object
        Object containing the mole amount of each species in the dry fuel.
    value : string
        'h' for enthalpy only, 'cp' for cp only, other for both.
    duty : float
        TODO

    Returns
    -------
    h : float
        Enthalpy [J] per 1 kg of fuel
    cp : float
        Specific heat capacity [J/(kmol.K)]
    '''
    si = mix.phase_index('solid')
    gi = mix.phase_index('gas')

    sMoles = mix.phase_moles(si)
    gMoles = mix.phase_moles(gi)

    h = (sMoles * mix.phase(si).enthalpy_mole
    + gMoles * mix.phase(gi).enthalpy_mole) / sum(mix.species_moles)
    
    cp = (sMoles * mix.phase(si).cp_mole
    + gMoles * mix.phase(gi).cp_mole) / sum(mix.species_moles)

    if value == 'h':
        return h
    elif value == 'cp':
        return cp
    else:
        return h, cp

# def enthalpy_of_formation(self, hhv):
#     '''
#     Estimate the standard enthalpy of formation of fuel [J/kg] from higher 
#     heating value and species composition.
    
#     Parameters
#     ----------
#     self : ndarray

#     Returns
#     -------
#     hfo : ndarray
#         standard enthalpy of formation of fuel [J/kg]
#     '''
#     f, stoic = get_fuel_db(self)
#     mol = f.species_moles # kmol
#     Mw = sum(mol*pp.Mw)
#     # standard enthalpy of formation [J/kg]
#     return (mol[pp.i_C]*pp.Hfo_CO2 + mol[pp.i_H]/2*pp.Hfo_H2Ol \
#             + mol[pp.i_N]*pp.Hfo_N2 + mol[pp.i_S]*pp.Hfo_SO2 \
#             + mol[pp.i_Cl]*pp.Hfo_ClO + mol[pp.i_SiO2]*pp.Hfo_SiO2 \
#             + mol[pp.i_CaO]*pp.Hfo_CaO + mol[pp.i_Al2O3]*pp.Hfo_Al2O3 \
#             + mol[pp.i_Fe2O3]*pp.Hfo_Fe2O3 + mol[pp.i_Na2O]*pp.Hfo_Na2O \
#             + mol[pp.i_K2O]*pp.Hfo_K2O + mol[pp.i_MgO]*pp.Hfo_MgO \
#             + mol[pp.i_P2O5]*pp.Hfo_P2O5 + mol[pp.i_TiO2]*pp.Hfo_TiO2 \
#             + mol[pp.i_SO3]*pp.Hfo_SO3 + mol[pp.i_Cr2O3]*pp.Hfo_Cr2O3 \
#             - stoic*pp.Hfo_O2 + hhv*1e6*Mw)/mol[pp.i_C]


# def simple_equilibrate_hp(self, moisture, fuel, air=zero, steam=zero, 
#                           P=ct.one_atm, duty=0):
#     """
#     Adiabatic multi-phase equilibrium calculation holding enthalpy and 
#     pressure fixed.
    
#     Use `equilibrate_hp' function for nonconventional fuels.

#     Parameters
#     ----------
#     self : ndarray
#         Mass fraction of fuel compounds in d.b. [kg/kg]
#     moisture : float
#         Mass fraction of moisture fuel [kg/kg]
#     fuel : float
#         Mass amount of fuel in d.b. [kg]
#     air : float
#         Mass amount of air [kg]
#     steam : float
#         Mass amount of steam [kg]
#     P : float
#         Pressure [Pa] (default = 1 atm)
#     duty : float
#         Duty fraction of outlet energy (default = 0)
#         Positive value means lost heat.

#     Returns
#     -------
#     content : object
#         Reactor state
#     inlet : float
#         Mole amount of inlet species [kmol]
#     outlet : float
#         Mole amount of outlet species [kmol]
#     T : float
#         Equilibrium temperature [K]
#     """
#     f = get_feed(self, moisture, fuel, air, steam)
#     # save initial composition
#     inlet = f.species_moles
#     # get enthalpy
#     H = f.H
#     # set desired condition
#     f.P = P
#     if duty != 0: f.H = (1-duty)*H
#     # calculate equilibrium
#     f.equilibrate('HP') #, solver='vcs', max_iter=200, estimate_equil=-1)
#     # save final composition
#     outlet = f.species_moles
#     T = f.T
#     return {'content':f, 'outlet':outlet, 'T':T, 'inlet':inlet}

# def equilibrate_hp(self, hfo, fuel, mw, moisture=zero, air=zero, steam=zero, 
#                    P=ct.one_atm, duty=0, guess=None, solver=0, disp=0):
#     '''
#     Non-isothermic multi-phase equilibrium calculation holding enthalpy and 
#     pressure fixed.
    
#     Use `simple_equilibrate_hp' function for conventional fuels.

#     Parameters
#     ----------
#     self : ndarray
#         Mass fraction of fuel compounds in d.b. [kg/kg]
#     moisture : float
#         Mass fraction of moisture fuel [kg/kg]
#     fuel : float
#         Mass amount of fuel in d.b. [kg]
#     mw : float
#         Molecular weight of fuel in d.b. [kg/kmol]
#     air : float
#         Mass amount of air [kg]
#     steam : float
#         Mass amount of steam [kg]
#     P : float
#         Pressure [Pa] (default = 1 atm)
#     duty : float
#         Duty fraction of outlet energy (default = 0)
#         Positive value means lost heat.
#     guess : float
#         Guess value of temperature for equilibrium calculations [K]
#     solver : integer
#         solver = 0, default calculation
#         solver = 1, scipy calculation
#     disp : integer
#         Display status notification of calculation.
#         Default = 0, no notification.

#     Returns
#     -------
#     content : objet
#         Reactor state    
#     inlet : float
#         Mole amount of inlet species [kmol]
#     outlet : float
#         Mole amount of outlet species [kmol]
#     T : float
#         Equilibrium temperature [K]
#     '''
#     f = get_feed(self, moisture, fuel, air, steam)
#     mole_moisture, mole_steam = get_water(self, moisture, fuel, steam)
#     # save initial composition
#     inlet = f.species_moles
#     # get moles of fuel
#     mole_fuel = fuel/mw
#     # get moles of air species
#     mole_O2 = inlet[pp.i_O2]
#     mole_N2 = inlet[pp.i_N2]
#     mole_Ar = inlet[pp.i_Ar]    
#     # inlet enthalpy [J/kmol]
#     inlet_h = (mole_fuel*hfo + mole_moisture*(pp.Hfo_H2Ol + pp.H_vap) \
#                 + mole_O2*pp.Hfo_O2 + mole_N2*pp.Hfo_N2 + mole_Ar*pp.Hfo_Ar \
#                 + mole_steam*pp.H_vap)/(mole_fuel + mole_moisture + mole_O2 \
#                 + mole_N2 + mole_Ar + mole_steam)
#     # use default guess value
#     if guess == None: guess = pp.To
#     # equilibrium calculation at T and P constant
#     def equilibrate_tp(self, T, P):
#         self.T = T
#         self.P = P
#         self.equilibrate('TP')
#         return self
#     # set phases
#     f = equilibrate_tp(f, guess, P)
#     # choose solver
#     # 0: own solver (default) (adapted from CATON et al., 2009)
#     # 1: scipy solver (scipy.optimize.minimize_scalar)
#     if solver == 0:
#         # default solver (adapted from CATON et al., 2009)
#         # set parameters to iterative calculation
#         dT = 50 # temperature increment
#         tol = 0.01 # tolerance
#         iters = 0 # initial iteration
#         # first state
#         # enthalpy and specific heat of outlet species
#         outlet_h, outlet_cp  = get_enthalpy(f,'h,cp')
#         # duty
#         outlet_h = (1-duty)*outlet_h
#         outlet_cp = (1-duty)*outlet_cp
#         # define the error
#         T_err0 = (outlet_h - inlet_h)/outlet_cp
#         # iterative calculation
#         # estimate equilibrium temperature and product composition
#         while (abs(T_err0) > tol):
#             guess += dT
#             f = equilibrate_tp(f, guess, P)
#             outlet_h, outlet_cp  = get_enthalpy(f,'h,cp')
#             # duty
#             outlet_h = (1-duty)*outlet_h
#             outlet_cp = (1-duty)*outlet_cp
#             T_err = (outlet_h - inlet_h)/outlet_cp
#             if (cmp(T_err, 0) != cmp(T_err0, 0)): # verify change of sign
#                 guess -= dT # go back to previous temperature
#                 dT *= 0.5 # decrease increment
#             else:
#                 # verify change of curve inclination after highest temperature
#                 if (abs(T_err) > abs(T_err0)):
#                     dT *= -1 # change of increment sign
#                 T_err0 = T_err # update value!
#             iters += 1 # counter
#             if iters == 200:
#                 print('maximum number of iterations reached')
#                 break
#             if disp == 2: 
#                 print('T = %4.2f, T_err = %0.4g, iters = %2.0f' %(guess,
#                                                                   T_err,iters))
#         if disp == 1:
#             print('T = %4.2f, T_err = %0.4g, iters = %2.0f' %(guess,
#                                                               T_err,iters))
#         T = f.T
#         outlet = f.species_moles
#     else:
#         # alternative solver (it uses minimize_scalar method)
#         def residual(x):
#             # set phases
#             f.T = x
#             f.P = P
#             f.equilibrate('TP')
#             # outlet enthalpy [J/kmol] with duty source
#             outlet_h  = (1-duty)*get_enthalpy(f,'h')
#             return (outlet_h - inlet_h)**2
#         # estimate equilibrium temperature
#         res = opt.minimize_scalar(residual,method='bounded',bounds=(200,6000),
#                                   bracket=(residual(1200),residual(3000)))
#         # estimate equilibrium product composition
#         T = res.x[0]
#         f = equilibrate_tp(f, T, P)
#         outlet = f.species_moles
#     return {'content':f, 'outlet':outlet, 'T':T, 'inlet':inlet}
