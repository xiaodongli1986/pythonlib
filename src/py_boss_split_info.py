
import stdA_py3

SixBinRedshift_NS = [0.2154098242, 0.3156545036, 0.3869159269, 0.4785988874, 0.5408209467, 0.6186561000]


SixBinRedshift_N = [0.21600642244568982, 0.31646577619767519, 0.38721946672537055, 0.47897499723114539, 0.54059313570594303, 0.61850917062883259]

SixBinRedshift_S = [0.21410525733734107, 0.31388039435706311, 0.38625224442260986, 0.47757238122842527, 0.54144257236794946, 0.61905709858835356]

Seven_Edges = [436.07,     775.25,     975.38,    1172.74,    1367.47,    1507.22,    1772.34 ]
Seven_Edges_z = [stdA_py3.get_z(0.26, -1, 0.7, xx) for xx in Seven_Edges]
Seven_Edges_Om0p3071 = [stdA_py3.comov_r(0.3071, -1, 0.7, xx) for xx in Seven_Edges_z]

Lpicola_npar_nhalo = ['npar','nhalo(1e4)',256,104,512,904,768,3438]
Cola_lc_npar_nhalo = ['npar','nhalo(1e4)',256,78,512,778,768,1954]

print("In Redshifts, the BOSS simulation is split according to: \n\t", [round(xx,5) for xx in Seven_Edges_z])
print("\nIn Omega=0.26 LCDM, the BOSS simulation is split according to: \n\t", [round(xx,3) for xx in Seven_Edges])
print("\nIn Omega=0.3071 LCDM, the BOSS simulation is split according to: \n\t", [round(xx,3) for xx in Seven_Edges_Om0p3071])

print('\nLpicola nbar to nhalo -- nhalo ~ nbar**(3.2) / box: 512 Mpc/h:\n\t')
for i in range(0,8,2):
    print(str(Lpicola_npar_nhalo[i])+'\t'+str(Lpicola_npar_nhalo[i+1])+'\n')

print('\nCola_lc nbar to nhalo -- nhalo ~ nbar**(3) / box: 512 Mpc/h:\n\t')
for i in range(0,8,2):
    print(str(Cola_lc_npar_nhalo[i])+'\t'+str(Cola_lc_npar_nhalo[i+1])+'\n')

print('''\n\n# Python script for split in any wCDM cosmology :
# run exportana3 in terminal
import stdA_py3
om = ...
w = ...

zs = [0.1500042520023839, 0.27410190459552325, 0.3510294396279122, 0.4300002403586951, 0.51130764221859, 0.5719528753056964, 0.6929005379550364]
print([round(stdA_py3.comov_r(om, w, 0.7, xx),3) for xx in zs])

''')
