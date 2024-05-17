import numpy       as np
import matplotlib as mpl
from HBVMod import HBVMod
import pymc3 as pm

forcing=pd.read_excel('E:/Learning Materials/TUD/hydrologic modeling/final assignment/data_hesperange.xls')
forcing = forcing.iloc[:,1:4]
forcing = forcing.as_matrix()
          #      Imax Ce Sumax beta Pmax   Tlag   Kf  Ks
Par = np.array([2,   .5,  100,   2,  .01,    5,   .05,  .001])
              #Si, Su,   Sf, Ss
Sin= np.array([0,  100,  0,  5  ])

Qm = HBVMod(Par,forcing,Sin, hydrograph='TRUE')


# MCMC method in HBV model
with pm.Model() as model:
	Imax = pm.Uniform('Imax', lower=0,upper=1)
	Ce = pm.Uniform('Ce', lower=0,upper=1)
	Sumax = pm.Uniform('Sumax', lower=0,upper=200)
	beta = pm.Uniform('beta', lower=0,upper=1)
	Pmax = pm.Uniform('Pmax', lower=0,upper=1)
	Tlag = pm.Uniform('Tlag', lower=0,upper=5)
	Kf = pm.Uniform('Kf', lower=0,upper=1)
	Ks = pm.Uniform('Ks', lower=0,upper=1)
	
