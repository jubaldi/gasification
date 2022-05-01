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

ashcomp.loc[len(ashcomp)] = ['Biomass', 'Wood',
                22.22, 43.03, 5.09, 3.44,
                2.85, 10.75, 6.07, 3.48,
                0.29,  2.78, 0]

ashcomp.loc[len(ashcomp)] = ['Biomass', 'Charcoal',
                22.22, 43.03, 5.09, 3.44,
                2.85, 10.75, 6.07, 3.48,
                0.29,  2.78, 0]

ashcomp.loc[len(ashcomp)] = ['Biomass', 'Straw',
                43.94, 14.13, 2.71, 1.42,
                1.35, 24.49, 4.66, 4.13,
                0.16,  3.01, 0]

ashcomp.loc[len(ashcomp)] = ['Biomass', 'Grass',
                46.18, 11.23, 1.39, 0.98,
                1.25, 24.59, 4.02, 6.62,
                0.08,  3.66, 0]

ashcomp.loc[len(ashcomp)] = ['Biomass', 'Agricultural',
                33.39, 14.86, 3.66, 3.26,
                2.29, 26.65, 5.62, 6.48,
                0.18,  3.61, 0]

ashcomp.loc[len(ashcomp)] = ['Biomass', 'Other',
                29.76, 25.27, 5.51, 4.00,
                2.48, 17.91, 5.42, 5.71,
                0.66,  3.28, 0]

ashcomp.loc[len(ashcomp)] = ['Other', 'Other',
                50, 20, 10, 10,
                10, 0, 0, 0,
                0,  0, 0]

ashcomp.to_csv('ashcomp.csv')

#print(ashcomp.head(3))
