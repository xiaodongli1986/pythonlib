#### import numpy as np
printstr = "# Example: py_illustris_hdf5_to_binary   > py_illustris_hdf5_to_binary.py"
print(printstr)

source_str = '''
import h5py
import os
import numpy as np
import struct
import sys
make_plot = True
if make_plot:
    import matplotlib.pyplot as ply

#----------------------------------
# 1. Basic names
#----------------------------------
sim, snap, snapdir = "Illustris-3", "135", "Snapshot"
#simname = sim+"_snap"+snap
#cic_nmesh = 1200

#----------------------------------
# 2. Output options
#----------------------------------
# 0: gas; 1: dm; 4: star; 5: BH
output_columns = {"PartType0": ["Coordinates", "Velocities", "Masses", "StarFormationRate", 
               "GFM_Metallicity", "NeutralHydrogenAbundance",
               "StarFormationRate", "InternalEnergy"],
 "PartType1": ["Coordinates", "Velocities"],
 "PartType4": ["Coordinates", "Velocities", "Masses"],
 "PartType5": ["Coordinates", "Velocities", "BH_Mass", "BH_Progs", 
               "BH_Mdot", "Masses", "HostHaloMass" ]
}

#----------------------------------
# 3. Files
#----------------------------------
files = os.popen("ls ../../data/Illustris/"+sim+"/"+snapdir+"/snap_"+snap+".*.hdf5").read().split()

#----------------------------------
# get some important numbers
#----------------------------------
##

snapf = h5py.File(files[0], "r")
header = dict(snapf["Header"].attrs.items()); 

boxsize = header["BoxSize"]

#----------------------------------
# Convert to binary
#----------------------------------
if True:
    for nowfile in files:
        snapf = h5py.File(nowfile, "r")
        for part_type in output_columns.keys():
            data = snapf[part_type]
            output_data = np.column_stack([np.array(data[col]).reshape(len(data[col]), -1) for col in output_columns[part_type]])
            output_file = nowfile.replace("data/Illustris", "data/Illustris_cic").replace(".hdf5", "").replace("snap_"+str(snap), "snap_"+str(snap)+"."+part_type)
            print(" -- writing ",output_data.shape," array to: \\n\\t", output_file)
            nowf = open(output_file, "wb")
            nowf.write(struct.pack("3d", len(output_data)+0., len(output_data[0])+0., boxsize))
            nowf.write(struct.pack(str(output_data.size)+"f", *(output_data.reshape(-1))))
            nowf.close()
'''
print(source_str)
nuisance_codes = '''

#----------------------------------
# 2. Split options
#----------------------------------
nsplit = 10

overlap = 0.1


#----------------------------------
# nuisance functions
#----------------------------------
def xyzminmax(xyz):
    return  xyz[:,0].min(), xyz[:,0].max(), xyz[:,1].min(), xyz[:,1].max(), xyz[:,2].min(), xyz[:,2].max()
def get_split_edges(xyz, split_edges):
    xmin, xmax, ymin, ymax, zmin, zmax = \
        xyz[:,0].min(), xyz[:,0].max(), xyz[:,1].min(), xyz[:,1].max(), xyz[:,2].min(), xyz[:,2].max()
    xsplit_edges, ysplit_edges, zsplit_edges = [], [], []
    for isplit, split_edge in enumerate(split_edges):
        for edge in split_edge:
            xyz1, xyz2 = edge
            if xyz1 < xmin < xyz2 or xyz1 < xmax < xyz2 or (xmin < xyz1 and xyz2 < xmax):
                xsplit_edges.append(isplit)
            if xyz1 < ymin < xyz2 or xyz1 < ymax < xyz2 or (ymin < xyz1 and xyz2 < ymax):
                ysplit_edges.append(isplit)
            if xyz1 < zmin < xyz2 or xyz1 < zmax < xyz2 or (zmin < xyz1 and xyz2 < zmax):
                zsplit_edges.append(isplit)
    xsplit_edges = list(set(xsplit_edges))
    ysplit_edges = list(set(ysplit_edges))
    zsplit_edges = list(set(zsplit_edges))
    return xmin, xmax, ymin, ymax, zmin, zmax, xsplit_edges, ysplit_edges, zsplit_edges
if False:
    bounds = np.linspace(0, boxsize, nsplit+1); delta = bounds[1]-bounds[0]

    edges = np.column_stack([bounds[:-1],bounds[1:]])

    split_edges = []
    for edge in edges:
        split_edges.append([ [edge[0]-delta * overlap, edge[1]+delta * overlap] ])
    #print(split_edges)
    def conv_split_edge(edge):
        print(edge); i1, i2 = edge[0]
        rlt = [ [max(i1,0), min(i2, boxsize)] ]
        if i1 <0: rlt.append([boxsize+i1, boxsize])
        if i2 > boxsize: rlt.append([0, i2-boxsize])
        return rlt
    split_edges = [conv_split_edge(X) for X in split_edges]
    for iX, X in enumerate(split_edges):
        print(iX, "-th split: \\n\\tsplit-ranges = ", X)          '''
