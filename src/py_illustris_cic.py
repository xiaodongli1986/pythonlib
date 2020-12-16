
print('''
from illustris_subboxes import *

nsplit=10
overlaprat=0.0
#runbashfile = sim+'_nc'+str(nc)+'_nsplit'+str(nsplit)+'_%.3f'%overlaprat+'_parttype'+str(parttype)+'.sh'
#print(' * create runbashfile: ', runbashfile)
#runbashf = open(runbashfile, 'w')

for sim_nc in [
       # ['Illustris-1', 50],
       # ['Illustris-1', 100],
       # ['Illustris-1', 200],
        ['Illustris-3', 50],
        ['Illustris-3', 100],
        ]:
    sim, nc = sim_nc
    for parttype in [0,1,4,5]:
        print('=================================================')
        print(' * sim, nc, parttype = ', sim, nc, parttype)

        runbashfile = sim+'_nc'+str(nc)+'_nsplit'+str(nsplit)+'_%.3f'%overlaprat+'_parttype'+str(parttype)+'.sh'
        print(' * create runbashfile: ', runbashfile)
        runbashf = open(runbashfile, 'w')

        illu = illustris_subbox(sim=sim, nsplit=nsplit, overlaprat=overlaprat, parttype=parttype)
        #illu.check_files(0,0,0)
        bashfile = illu.do_cic(nc); runbashf.write('sh '+bashfile+'\\n')
        runbashf.close()
        ''')
