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
# QSTN: 2022/04/15 Ask professor if no data should return 0 or np.nan.
# TODO: 2022/04/27 Create functions: from fuelID and mass get Mix, from Mix get mass.
# TODO: 2022/04/18 Biochemical composition should be in csv.
# TODO: 2022/05/01 Find ash composition for petroleum coke.
# TODO: 2022/05/01 Add a cantera mixture for steam + water.
# TODO: 2022/05/01 Find heating values HHV and LHV from correlations.
# QSTN: 2022/05/01 Should we assume that remaining mass is fixed carbon?

# TODO LATER: 2022/05/01 Get heat of formation.

# MISTAKE: (1 - ash) ???

#==============================================================================
# import libraries
#==============================================================================
import numpy as np
import pandas as pd
import json
import pp2 as pp

#==============================================================================
# special functions
#==============================================================================

# Create a dataframe from the csv
with open('fuels-updated.csv','r') as f1:
    fuels = pd.read_csv(f1, sep=',', header=0, index_col=0)
    f1.close()

with open('fuels.json','r') as fuelFile:
    fuelData = json.load(fuelFile)
    fuelFile.close()

# create a list of the fuel names from the first column of the dataframe
myID = fuels.index.values

# Description, Type, Category, etc: read directly using fuels.loc
# Moisture, FC, VM: read directly then divide by 100
# HHV, LHV: read directly
# C, H, O, N, S, Cl: read directly then divide by 100

with open('ashes.json', 'r') as ashFile:
    ashData = json.load(ashFile)
    ashFile.close()

def moisture(fuelID, basis="db"):
    '''
    Gets the moisture content for a given fuel. 
    The fuel must be available in the database (file: 'fuels-updated.csv').

    Parameters
    ----------
    fuelID : str
        The fuel ID as given by the CSV.
    basis : str
        "db" = dry basis (default), "wb" = wet basis.
    
    Returns
    -------
    moist : float
        Moisture content [kg/kg]
    '''
    if fuelID not in myID:
        raise ValueError('Fuel ID not found in database.')
    
    # Read values from CSV
    rMoist = fuels.loc[fuelID]['Moisture']
    rMoistBasis = fuels.loc[fuelID]['Mbasis']

    # Check if empty
    if pd.isnull(rMoist):
        moist = 0
    else:
        moist = rMoist/100 # Convert % to fraction
    
    # Convert data to dry basis if not already
    if rMoistBasis == "wb":
        moist = moist / (1 - moist) # moist is now in db
    
    # Convert dry basis to desired basis for returning
    if basis == "wb":
        moist = moist / (1 + moist) # moist is now in wb
    
    return moist

def proximate(fuelID, basis="db"):
    '''
    Gets the proximate analysis for a given fuel. 
    The fuel must be available in the database (file: 'fuels-updated.csv').

    Parameters
    ----------
    fuelID : str
        The fuel ID as given by the CSV.
    basis : str
        "db" = dry basis (default), "wb" = wet basis.

    Returns
    -------
    prox : dict
        Proximate analysis of fuel [kg/kg]:
        {"Ash": Ash, "FC": Fixed Carbon, "VM": Volatile Matter}
    '''
    if fuelID not in myID:
        raise ValueError('Fuel ID not found in database.')

    # Read values from CSV
    rProxBasis = fuels.loc[fuelID]['Pbasis']

    prox = {}

    # Read values and check if empty
    for p in ["FC", "Ash", "VM"]:
        rProx = fuels.loc[fuelID][p]
        if pd.isnull(rProx):
            if p == "VM":
                prox[p] = 1 - prox["FC"] - prox["Ash"] # Assume remaining mass is VM
            else:
                prox[p] = 0
        else:
            prox[p] = rProx/100
    
    # Grab moisture in wb
    moistWB = moisture(fuelID, basis="wb")

    # Convert data to dry basis if not already
    if rProxBasis == "wb":
        for p in prox:
            prox[p] = prox[p] / (1 - moistWB)
    
    # Convert dry basis to desired basis for returning
    if basis == "wb":
        for p in prox:
            prox[p] = prox[p] * (1 - moistWB)

    return prox

def ashComp(fuelID):
    '''
    Gets the ash composition for a given fuel. 
    The fuel must be available in the database (file: 'fuels.csv').
    
    Default values are used if they are not available. Those values are
    chosen from VASSILEV et al. (2013) according to type and category of fuel.

    Parameters
    ----------
    fuelID : str
        The fuel ID as given by the CSV.

    Returns
    -------
    comp : ndarray
        Composition of ash [kg/kg]:
        [SiO2, CaO, Al2O3, Fe2O3, Na2O, K2O, MgO, P2O5, TiO2, SO3, Cr2O3]
    '''
    if fuelID not in myID:
        raise ValueError('Fuel ID not found in database.')
    
    # Read values from CSV
    ash = proximate(fuelID)["Ash"]
    rComp = fuels.loc[fuelID]['SiO2':'Cr2O3']/100

    # If ash fraction is 0, ash composition doesn't matter
    if ash == 0:
        comp = np.ones(len(rComp))/len(rComp)
    else:
        comp = np.zeros(len(rComp))

    # If ash composition is not given by csv,
    # use values from Vassilev2013 given Type and Category
        if pd.isnull(rComp['SiO2']) and pd.isnull(rComp['CaO']):
            rType = fuels.loc[fuelID]['Type']
            rCat = fuels.loc[fuelID]['Category']
            if rType in ashData: 
                fType = rType
            else:
                fType = 'Other'
            if rCat in ashData[fType]:
                fCat = rCat
            else:
                fCat = 'Other'
            compDict = ashData[fType][fCat]
            comp = np.array(list(compDict.values()))/100

    # If ash composition is given by csv, use it
        else:
            for index, species in enumerate(rComp):
                if pd.isnull(species):
                    comp[index] = 0
                else:
                    comp[index] = species

    # Normalize composition
    comp = comp/np.sum(comp)

    return comp

def ultimate(fuelID, basis="db"):
    '''
    Gets the ultimate analysis for a given fuel. 
    The fuel must be available in the database (file: 'fuels-updated.csv').

    Parameters
    ----------
    fuelID : str
        The fuel ID as given by the CSV.
    basis : str
        "db" = dry basis (default), "wb" = wet basis, "daf" = dry ash-free basis.
    
    Returns
    -------
    ult : ndarray
        Ultimate analysis of fuel [kg/kg]: 
        [C, H, O, N, S, Cl]
    '''
    if fuelID not in myID:
        raise ValueError('Fuel ID not found in database.')
    
    # Read values from CSV
    rUltBasis = fuels.loc[fuelID]['Ubasis']
    rC = fuels.loc[fuelID]['C']
    rH = fuels.loc[fuelID]['H']
    rO = fuels.loc[fuelID]['O']
    rN = fuels.loc[fuelID]['N']
    rS = fuels.loc[fuelID]['S']
    rCl = fuels.loc[fuelID]['Cl']

    ult = np.zeros(6)

    # Check if empty
    for i, elem in enumerate([rC, rH, rO, rN, rS, rCl]):
        if pd.isnull(elem):
            ult[i] = 0
        else:
            ult[i] = elem/100

    moistWB = moisture(fuelID, basis="wb")
    ashDB = proximate(fuelID, basis="db")["Ash"]
    
    # Convert data to dry basis if not already
    if rUltBasis == "wb":
        for i, elem in enumerate([rC, rH, rO, rN, rS, rCl]):
            ult[i] = ult[i] / (1 - moistWB)
    elif rUltBasis == "daf":
        for i, elem in enumerate([rC, rH, rO, rN, rS, rCl]):
            ult[i] = ult[i] * (1 - ashDB)
        
    # Convert dry basis to desired basis for returning
    if basis == "wb":
        for i, elem in enumerate([rC, rH, rO, rN, rS, rCl]):
            ult[i] = ult[i] * (1 - moistWB)
    elif basis == "daf":
        for i, elem in enumerate([rC, rH, rO, rN, rS, rCl]):
            ult[i] = ult[i] / (1 - ashDB)
    
    return ult

def fuelComp(fuelID):
    '''
    Gets the full mass composition for the dry fuel.
    The fuel must be available in the database (file: 'fuels.csv').

    Parameters
    ----------
    fuelID : str
        The fuel ID as given by the CSV.

    Returns
    -------
    fuelComp : dict
        A dictionary representing the mass fraction of each species.
    '''

    ash = proximate(fuelID)["Ash"]
    ult = ultimate(fuelID)
    ashCompos = ashComp(fuelID)

    fuelComp = {}

    for index, species in enumerate(["C(gr)", "H", "O", "N", "S", "CL"]):
        fuelComp[species] = ult[index]

    for index, species in enumerate(['SiO2(hqz)', 'CaO(s)', 'AL2O3(a)', 'Fe2O3(s)', 
    'Na2O(c)', 'K2O(s)', 'MgO(s)', 'P2O5', 'TiO2(ru)', 'SO3', 'Cr2O3(s)']):
        fuelComp[species] = ashCompos[index]*ash

    # Normalize each fuel composition
    s = sum(fuelComp.values())
    for (key, value) in fuelComp.items():
        fuelComp[key] = value/s

    return fuelComp

def fuelComp2(ash, ult, ashCompos):
    '''
    Get the full mass composition for the dry fuel.
    The fuel does not need to be in the database.

    Parameters
    ----------
    ash : float
        Ash fraction of fuel, in dry basis [kg/kg]
    ult : ndarray
        Ultimate analysis of fuel, in dry basis [kg/kg]
    ashCompos : ndarray
        Composition of ash [kg/kg]

    Returns
    -------
    fuelComp : dict
        A dictionary representing the mass fraction of each species.
    '''

    fuelComp = {}

    for index, species in enumerate(["C(gr)", "H", "O", "N", "S", "CL"]):
        fuelComp[species] = ult[index]

    for index, species in enumerate(['SiO2(hqz)', 'CaO(s)', 'AL2O3(a)', 'Fe2O3(s)', 
    'Na2O(c)', 'K2O(s)', 'MgO(s)', 'P2O5', 'TiO2(ru)', 'SO3', 'Cr2O3(s)']):
        fuelComp[species] = ashCompos[index]*ash

    # Normalize each fuel composition
    s = sum(fuelComp.values())
    for (key, value) in fuelComp.items():
        fuelComp[key] = value/s

    return fuelComp

def addToDatabase(fuelID, Info, moist, prox, HV, biochem, ult, ashC):
    '''
    Adds a new fuel to the database.
    The fuel must not be available in the database (file: 'fuels-updated.csv').

    Parameters
    ----------
    fuelID : str
        The fuel ID to be registered in the CSV.
    Info : dict
        A dictionary containing identifying information about the fuel.
        Can contain: "Description", "Type", "Category", "Reference", "Year", "DOI"
    moist : list
        Moisture content information (%m/m).
        [Mbasis, Moisture]
    prox : list
        Proximate analysis information (%m/m).
        [Pbasis, FC, VM, Ash]
    HV : list
        Heating values information (MJ/kg).
        [HHV, LHV]
    biochem : list
        Biochemical composition information (%m/m).
        [Cellulose, Hemicellulose, Lignin]
    ult : list
        Ultimate analysis information (%m/m).
        [Ubasis, C, H, O, N, S, Cl]
    ashC : list
        Ash composition information (%m/m of ash).
        [SiO2, CaO, Al2O3, Fe2O3, Na2O, K2O, MgO, P2O5, TiO2, SO3, Cr2O3]
    
    Returns
    -------
    None
    '''
    # Check if fuelID already exists
    if fuelID in fuels.index:
        raise(ValueError("Fuel ID already exists in database."))
    
    # Create a new empty row in the CSV
    fuels.loc[fuelID] = pd.Series(dtype=object)
    for key, value in Info.items():
        # Check if key is already a column
        if key in fuels.columns:
            fuels.loc[fuelID, key] = value
    fuels.loc[fuelID,'Mbasis'] = moist[0]
    fuels.loc[fuelID,'Moisture'] = moist[1]
    fuels.loc[fuelID,'Pbasis'] = prox[0]
    fuels.loc[fuelID,'FC'] = prox[1]
    fuels.loc[fuelID,'VM'] = prox[2]
    fuels.loc[fuelID,'Ash'] = prox[3]
    fuels.loc[fuelID,'HHV'] = HV[0]
    fuels.loc[fuelID,'LHV'] = HV[1]
    fuels.loc[fuelID,'Cellulose'] = biochem[0]
    fuels.loc[fuelID,'Hemicellulose'] = biochem[1]
    fuels.loc[fuelID,'Lignin'] = biochem[2]
    fuels.loc[fuelID,'Ubasis'] = ult[0]
    fuels.loc[fuelID,'C'] = ult[1]
    fuels.loc[fuelID,'H'] = ult[2]
    fuels.loc[fuelID,'O'] = ult[3]
    fuels.loc[fuelID,'N'] = ult[4]
    fuels.loc[fuelID,'S'] = ult[5]
    fuels.loc[fuelID,'Cl'] = ult[6]
    fuels.loc[fuelID,'SiO2'] = ashC[0]
    fuels.loc[fuelID,'CaO'] = ashC[1]
    fuels.loc[fuelID,'Al2O3'] = ashC[2]
    fuels.loc[fuelID,'Fe2O3'] = ashC[3]
    fuels.loc[fuelID,'Na2O'] = ashC[4]
    fuels.loc[fuelID,'K2O'] = ashC[5]
    fuels.loc[fuelID,'MgO'] = ashC[6]
    fuels.loc[fuelID,'P2O5'] = ashC[7]
    fuels.loc[fuelID,'TiO2'] = ashC[8]
    fuels.loc[fuelID,'SO3'] = ashC[9]
    fuels.loc[fuelID,'Cr2O3'] = ashC[10]
    fuels.to_csv('fuels-updated.csv')

fuelID = "Test"
Info = {"Description": "Test", "Type": "Test", "Category": "Test", "Reference": "Test", "Year": "Test", "DOI": "Test"}
moist = ["wb", 1]
prox = ["wb", 1, 1, 1]
HV = [1, 1]
biochem = [1, 1, 1]
ult = ["wb", 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
ashC = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

addToDatabase(fuelID, Info, moist, prox, HV, biochem, ult, ashC)

def HV(fuelID, type='both', moist=0.0):
    '''
    Gets the heating values (HHV or LHV) for a given fuel.
    The fuel must be available in the database (file: 'fuels.csv').

    Parameters
    ----------
    fuelID : str
        The fuel ID as given by the CSV.
    type : str
        Either 'HHV', 'LHV' or other
    moist : float
        The moisture content of the fuel, dry basis [kg/kg]

    Returns
    -------
    HV : float | dict
        The heating value [kJ/kg]
    '''
    if fuelID not in myID:
        raise ValueError('Fuel ID not found in database.')

    # Check if there is available HHV or LHV data
    if "HHV" not in fuelData[fuelID]:
        isHHV = False
    else:
        if not fuelData[fuelID]["HHV"]:
            isHHV = False
        else:
            isHHV = True
    
    if "LHV" not in fuelData[fuelID]:
        isLHV = False
    else:
        if not fuelData[fuelID]["LHV"]:
            isLHV = False
        else:
            isLHV = True
    
    # Read hydrogen content for use in correlations
    if not fuelData[fuelID]["H"]: 
        H = 0
    else: 
        H = float(fuelData[fuelID]["H"])/100 

    if isLHV and (not isHHV):
        LHV = float(fuelData[fuelID]["LHV"])
        HHV = (LHV + 2.258*moist)/(1 - moist) + 20.1790296248*H
    elif isHHV and (not isLHV):
        HHV = float(fuelData[fuelID]["HHV"])
        LHV = (HHV - 20.1790296248*H)*(1 - moist) - 2.258*moist
    elif (not isLHV) and (not isHHV):
        HHV = 20 # TODO: Use correlations for HV calculations
        LHV = 20
    else:
        HHV = float(fuelData[fuelID]["HHV"])
        LHV = float(fuelData[fuelID]["LHV"])

    LHV *= (1 - moist)
    HHV *= (1 - moist)

    # Get HHV if type is 'HHV'
    if type == 'HHV':
        HV = HHV
    # Get LHV if type is 'LHV'
    elif type == 'LHV':
        HV = LHV
    # Get both if type is 'both'
    else:
        HV = {'HHV': HHV, 'LHV': LHV}

    return HV
