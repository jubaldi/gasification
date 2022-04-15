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
        The fuel ID.
        
    Returns
    -------
    description : str | list
        The description of the fuel.
    """
    if isinstance(fuelID, str):
        return fuels.loc[fuelID]['Description']
    else:
        return [fuels.loc[i]['Description'] for i in fuelID]

def fuelType(fuelID):
    """
    Get the fuel type for a single fuel or list of fuels.
    These fuels must be available in the database (file: 'fuels.csv').
        
    Parameters
    ----------
    fuelID : str | list
        The fuel ID.
        
    Returns
    -------
    fuelType : str | list
        The fuel type.
    """
    if isinstance(fuelID, str):
        return fuels.loc[fuelID]['Type']
    else:
        return [fuels.loc[i]['Type'] for i in fuelID]

def category(fuelID):
    """
    Get the category for a single fuel or list of fuels.
    These fuels must be available in the database (file: 'fuels.csv').
            
    Parameters
    ----------
    fuelID : str | list
        The fuel ID.
            
    Returns
    -------
    category : str | list
        The category.
    """
    if isinstance(fuelID, str):
        return fuels.loc[fuelID]['Category']
    else:
        return [fuels.loc[i]['Category'] for i in fuelID]

def reference(fuelID):
    """
    Get the reference for a single fuel or list of fuels.
    These fuels must be available in the database (file: 'fuels.csv').
            
    Parameters
    ----------
    fuelID : str | list
        The fuel ID.
            
    Returns
    -------
    reference : str | list
        The reference.
    """
    if isinstance(fuelID, str):
        return fuels.loc[fuelID]['Reference']
    else:
        return [fuels.loc[i]['Reference'] for i in fuelID]

def year(fuelID):
    """
    Get the year for a single fuel or list of fuels.
    These fuels must be available in the database (file: 'fuels.csv').
            
    Parameters
    ----------
    fuelID : str | list
        The fuel ID.
            
    Returns
    -------
    year : int | list
        The year.
    """
    if isinstance(fuelID, str):
        return int(fuels.loc[fuelID]['Year'])
    else:
        return [int(fuels.loc[i]['Year']) for i in fuelID]

def DOI(fuelID):
    """
    Get the DOI for a single fuel or list of fuels.
    These fuels must be available in the database (file: 'fuels.csv').
            
    Parameters
    ----------
    fuelID : str | list
        The fuel ID.
            
    Returns
    -------
    DOI : str | list
        The DOI.
    """
    if isinstance(fuelID, str):
        return fuels.loc[fuelID]['DOI']
    else:
        return [fuels.loc[i]['DOI'] for i in fuelID]

def lowerMoisture(fuelID):
    """
    Get the lower moisture for a single fuel or list of fuels.
    These fuels must be available in the database (file: 'fuels.csv').
            
    Parameters
    ----------
    fuelID : str | list
        The fuel ID.
            
    Returns
    -------
    lowerMoisture : float | array
        The lower moisture fraction [kg/kg] for each fuel
    """
    lowerMoisture_df = fuels.loc[fuelID]['Lower moisture']
    if isinstance(fuelID, str):
        if lowerMoisture_df == '':
            return np.nan
        else:
            return float(lowerMoisture_df)/100
    else:
        lowerMoisture = []
        for value in lowerMoisture_df:
            if value == '':
                lowerMoisture.append(np.nan)
            else:
                lowerMoisture.append(float(value)/100)
        return lowerMoisture
    
def ash(fuelID):
    """
    Get the ash fraction value for a list of fuels. 
    These fuels must be available in the database (file: 'fuels.csv').
    
    Default values are used if they are not available. Those values are
    chosen from VASSILEV et al. (2013) according to type and category of fuel.

    Parameters
    ----------
    fuelID : str | list
        The fuel ID.

    Returns
    -------
    ash : float | list
        Ash fraction [kg/kg]
    composition : ndarray | list
        Composition of ash fraction [kg/kg], of SiO2, CaO, Al2O3, Fe2O3, Na2O, 
        K2O, MgO, P2O5, TiO2, SO3, Cr2O3 respectively.
    
    References
    ----------
    VASSILEV, S. V.; BAXTER, D.; ANDERSEN, L. K.; VASSILEVA, C. G. An overview 
    of the composition and application of biomass ash. Part 1. Phase-mineral 
    and chemical composition and classification. Fuel. v. 105, p. 40-76, 2013.
    """

    myFuels = fuels.loc[fuelID]
    ash_df = myFuels['Ash']/100.0
    return ash
