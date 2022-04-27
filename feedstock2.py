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

#==============================================================================
# import libraries
#==============================================================================
import numpy as np
import pandas as pd
import csv
import pp

# Create a dataframe from the csv
with open('fuels.csv','r') as f:
    fuels = pd.read_csv(f, sep=',', header=0, index_col=0)
    f.close()

# create a list of the fuel names from the first column of the dataframe
myID = fuels.index.values

def description(fuelID):
    """
    Get the description for a single fuel or list of fuels.
    These fuels must be available in the database (file: 'fuels.csv'). 
        
    Parameters
    ----------
    fuelID : str | list
        The fuel ID or list of fuel IDs.
        
    Returns
    -------
    description : list
        A list containing the description for each fuel as a string.
    """
    if isinstance(fuelID, str):
        return [fuels.loc[fuelID]['Description']]
    else:
        return [fuels.loc[i]['Description'] for i in fuelID]

def fuelType(fuelID):
    '''
    Get the fuel type for a list of fuels. 
    These fuels must be available in the database (file: 'fuels.csv').
    
    Parameters
    ----------
    fuelID : str | list
        The fuel ID or list of fuel IDs.

    Returns
    -------
    fuelType : list
        Fuel type
    '''
    if isinstance(fuelID, str):
        return [fuels.loc[fuelID]['Type']]
    else:
        return [fuels.loc[i]['Type'] for i in fuelID]

def category(fuelID):
    """
    Get the category for a list of fuels. 
    These fuels must be available in the database (file: 'fuels.csv').

    Parameters
    ----------
    fuelID : str | list
        The fuel ID or list of fuel IDs.

    Returns
    -------
    category : list
        Category
    """
    if isinstance(fuelID, str):
        return [fuels.loc[fuelID]['Category']]
    else:
        return [fuels.loc[i]['Category'] for i in fuelID]

def moisture(fuelID):
    """
    Get the moisture fraction value for a list of fuels. 
    These fuels must be available in the database (file: 'fuels.csv').

    Parameters
    ----------
    fuelID : str | list
        The fuel ID or list of fuel IDs.

    Returns
    -------
    moist : list of float
        Moisture fraction [kg/kg]
    """
    if isinstance(fuelID, str):
        return [fuels.loc[fuelID]['Lower moisture']/100]
    else:
        return [fuels.loc[i]['Lower moisture']/100 for i in fuelID]

def ash(fuelID):
    '''
    Get the ash fraction value for a list of fuels. 
    These fuels must be available in the database (file: 'fuels.csv').
    
    Default values are used if there are not available. Those values are
    chosen from VASSILEV et al. (2013) according to type and category of fuel.

    Parameters
    ----------
    fuelID : str | list
        The fuel ID or list of fuel IDs.

    Returns
    -------
    ash : list of float
        Ash fraction [kg/kg]
    composition : ndarray
        Composition of ash fraction [kg/kg]
    '''
    if isinstance(fuelID, list):
        ash = [fuels.loc[i]['Ash']/100 for i in fuelID]
        composition = 0
    return {'fraction':ash, 'composition': composition}

# TODO: 2022/04/18 Ash composition and biochemical composition should be in csv.