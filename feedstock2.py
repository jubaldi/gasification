#!/usr/bin/env python

"""This script defines functions to get values from a fuel database.

The fuel database is a predefined CSV file ('fuels.csv'). The file can be 
edited using any spreadsheet editor to add new fuel compounds.

@author = Rodolfo Rodrigues
@contact = rodolfo.rodrigues@ufsm.br
@data = April, 2012, rev.: June, 2013 (adapted to use cython Cantera)
"""

# FIXME: 2013/05/25 Not possible to use a list of single value as parameter to 
# functions.
# FIXME: 2022/04/15 Ask professor if no data should return 0 or np.nan.
# TODO: 2022/04/27 Create functions: from fuelID and mass get Mix, from Mix get mass.
# TODO: 2022/04/18 Ash composition and biochemical composition should be in csv.
# TODO: 2022/05/01 Find ash composition for petroleum coke.
# TODO: 2022/05/01 Add a cantera mixture for steam + water.

#==============================================================================
# import libraries
#==============================================================================
import numpy as np
import pandas as pd
import csv
import pp

# Create a dataframe from the csv
with open('fuels.csv','r') as f1:
    fuels = pd.read_csv(f1, sep=',', header=0, index_col=0)
    f1.close()

# create a list of the fuel names from the first column of the dataframe
myID = fuels.index.values

# Description, Type, Category, etc: read directly using fuels.loc
# Moisture, FC, VM: read directly then divide by 100
# HHV, LHV: read directly
# C, H, O, N, S, Cl: read directly then divide by 100

with open('ashcomp.csv', 'r') as f2:
    ashcomp = pd.read_csv(f2, sep=',', header=0, index_col=0)
    f2.close()

def ash(fuelID):
    '''
    Get the ash fraction value for a given fuel. 
    The fuel must be available in the database (file: 'fuels.csv').
    
    Default values are used if they are not available. Those values are
    chosen from VASSILEV et al. (2013) according to type and category of fuel.

    Parameters
    ----------
    fuelID : str
        The fuel ID as given by the CSV.

    Returns
    -------
    ash : float
        Ash fraction [kg/kg]
    composition : ndarray
        Composition of ash fraction [kg/kg]
    '''
    ash = fuels.loc[fuelID]['Ash']/100
    comp = fuels.loc[fuelID]['SiO2':'Cr2O3']/100
    if pd.isnull(comp.loc['SiO2']) and pd.isnull(comp.loc['CaO']):
        if ash == 0:
            fuelType = 'Other'
            fuelCategory = 'Other'
        else:
            fuelType = fuels.loc[fuelID]['Type']
            if pd.isnull(fuels.loc[fuelID]['Category']):
                fuelCategory = 'Other'
            else:
                fuelCategory = fuels.loc[fuelID]['Category']
        a1 = ashcomp.loc[ashcomp['Type']==fuelType]
        a2 = a1.loc[a1['Category']==fuelCategory]
        composition = a2.loc['SiO2':'Cr2O3']/100
    else:
        composition = comp
    return ash, composition

print(ash('Cedar'))