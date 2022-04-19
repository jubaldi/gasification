import pandas as pd
import numpy as np

# create a new empty dataframe
ashcomp = pd.DataFrame(columns=['Type', 'Category',
 'SiO2', 'CaO', 'Al2O3', 'Fe2O3', 'Na2O', 'K2O', 
 'MgO', 'P2O5', 'TiO2', 'SO3', 'Cr2O3'])

# add a line to dataframe:
ashcomp.loc[len(ashcomp)] = ['Coal', 'Subbituminous',
                54.74, 7.05, 22.86, 5.30,
                1.09, 1.67,  2.14, 0.08,
                1.00, 4.07,     0]

ashcomp.loc[len(ashcomp)] = ['Coal', 'Other',
                54.06, 6.57, 23.18, 6.85,
                0.82, 1.60,  1.83, 0.50,
                1.05, 3.54,     0]

print(ashcomp.head(3))
