
print("Creating a bash file to conduct rockstar halo searching for multiple snapshots..\n\t 1.py_rockstar_multisnap.sh")

nowf = open("1.py_rockstar_multisnap.sh", 'w')
nowf.write('''py_rockstar \*snap\*a.\*
py_rockstar \*snap\*b.\*
py_rockstar \*snap\*c.\*
py_rockstar \*snap\*d.\*
py_rockstar \*snap\*e.\*
py_rockstar \*snap\*f.\*
py_rockstar \*snap\*g.\*
py_rockstar \*snap\*h.\*
py_rockstar \*snap\*i.\*
py_rockstar \*snap\*j.\*
py_rockstar \*snap\*k.\*
py_rockstar \*snap\*l.\*
py_rockstar \*snap\*m.\*
py_rockstar \*snap\*n.\*
py_rockstar \*snap\*o.\*
py_rockstar \*snap\*p.\*
py_rockstar \*snap\*q.\*
py_rockstar \*snap\*r.\*
py_rockstar \*snap\*s.\*
py_rockstar \*snap\*t.\*
py_rockstar \*snap\*u.\*
py_rockstar \*snap\*v.\*
py_rockstar \*snap\*w.\*
py_rockstar \*snap\*x.\*
''')
nowf.close()
