

print("We will create a fake headinfo file:\n\tfake_headinfo.txt")

nowf = open("fake_headinfo.txt", 'w')

nowf.write("# ntotal+0., headinfo.boxsize, nonzeromass, headinfo.redshift, headinfo.Omega0, headinfo.HubbleParam\n")
nowf.write("0 256  1.0E+8  -1  0.3071 0.6787 0  ")

nowf.close()


