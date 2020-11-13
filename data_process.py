import happi
import matplotlib.pyplot as plt
from numpy import array, pi

path         = "/home/men/files/test/"
S            = happi.Open(path,verbose=False)
dt           = S.namelist.Main.timestep
#print(S.namelist.Main.timestep)




Nit          = array(S.Scalar("Dens_ion2").getData())
time         = array(S.Scalar("Dens_ion2").getTimes())

Ex           = array(S.Field.Field0.Ex(timesteps=2000).getData()[0])
Nix          = array(S.Field.Field0.Rho_ion(timesteps=20000).getData()[0])
Nix0         = array(S.Field.Field0.Rho_ion(timesteps=1).getData()[0])

Ne           = array(S.Field.Field0.Rho_eon(timesteps=2000).getData()[0])
ion2         = S.ParticleBinning(0)
ele          = S.ParticleBinning(1)

##ion2.plot(figure=2)
#ele.plot(figure=1)
fig = plt.figure(figsize=(4,3),dpi=300);
plt.plot(Nix,"r")
plt.plot(Nix0,"k")
plt.xlabel('x')
plt.ylabel('rho_i')

fig = plt.figure(figsize=(4,3),dpi=300);
plt.plot(time,Nit,"r")

