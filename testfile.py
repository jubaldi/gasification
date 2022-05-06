import numpy as np
import pandas as pd
import feedstock as fs1
import feedstock2 as fs
import pp as pp1
import pp2 as pp
import gasifier2 as g
import fuel as fu
import energy as en
import outputs as op

fuel1 = 'SMC'
fuel2 = 'RS'
T = 940 + 273.15
P = 1
ER = 0.32
SR = 0.49/(pp.Mw['H2O']/pp.Mw['C'])
m1 = 2.4
m2 = 0.6
fm1 = fs.getFuelMix(fuel1, m1)
fm2 = fs.getFuelMix(fuel2, m2)
b = fs.blend(fm1, fm2)
air = fs.ERtoair(b, ER)
steam = fs.SRtosteam(b, SR)
f = g.getFeed(b, 0, air, steam)
#print(fs.getFuelMass(f)/(air+steam+m1+m2))

a = g.isotCogasification(fuel1, fuel2, m1, m2/m1, 0, T, P, air, steam, 'air', 'steam')
print(a.species_moles[pp.i['CO2']]/sum(a.species_moles)*100)
print(a.species_moles[pp.i['H2']]/sum(a.species_moles)*100)
print(a.species_moles[pp.i['N2']]/sum(a.species_moles)*100)
print(a.species_moles[pp.i['CO']]/sum(a.species_moles)*100)
print(a.species_moles[pp.i['CH4']]/sum(a.species_moles)*100)