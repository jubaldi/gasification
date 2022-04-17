import numpy as np
import pandas as pd
import feedstock as fs1
import feedstock2 as fs2
import pp
import pp2

print(0, pp.Hfo_air == pp2.Hfo_air)
print(1, pp.H_vap_mass == pp2.H_vap_mass)
print(2, pp.H_vap == pp2.H_vap)
print(3, pp.Mw_f == pp2.Mw_f)
print(4, pp.Mw_air == pp2.Mw_air)
print(5, pp.i_SiO2 == pp2.i['SiO2(hqz)'])
print(6, pp.Hfo_Cr2O3 == pp2.Hfo['Cr2O3(s)'])
print(7, pp.Hfo_H2Ol == pp2.Hfo['H2O(l)'])