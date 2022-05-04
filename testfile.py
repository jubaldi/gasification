import numpy as np
import pandas as pd
import feedstock as fs1
import feedstock2 as fs2
import pp
import pp2
import gasifier as g

#print(0, pp.Hfo_air == pp2.Hfo_air)
#print(1, pp.H_vap_mass == pp2.H_vap_mass)
#print(2, pp.H_vap == pp2.H_vap)
#print(3, pp.Mw_f == pp2.Mw_f)
#print(4, pp.Mw_air == pp2.Mw_air)
#print(5, pp.i_SiO2 == pp2.i['SiO2(hqz)'])
#print(6, pp.Hfo_Cr2O3 == pp2.Hfo['Cr2O3(s)'])
#print(7, pp.Hfo_H2Ol == pp2.Hfo['H2O(l)'])

fuel = ['Bagasse1','Bagasse0']
frac = fs1.fraction(fuel)
a = g.get_fuel_db(frac['mass'][0])
mix = fs2.getFuelMix(fuel[0], 1)
b = fs2.stoichO2(mix)

print(a[1])
print(b)