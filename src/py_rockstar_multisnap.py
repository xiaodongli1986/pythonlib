
print("Creating a bash file to conduct rockstar halo searching for multiple snapshots..\n\t 1.py_rockstar_multisnap.sh")

nowf = open("1.py_rockstar_multisnap.sh", 'w')
nowf.write('''
ranseed=01810
basename=sim
py_rockstar ${basename}\*snap${ranseed}a.\*
py_rockstar ${basename}\*snap${ranseed}b.\*
py_rockstar ${basename}\*snap${ranseed}c.\*
py_rockstar ${basename}\*snap${ranseed}d.\*
py_rockstar ${basename}\*snap${ranseed}e.\*
py_rockstar ${basename}\*snap${ranseed}f.\*
py_rockstar ${basename}\*snap${ranseed}g.\*
py_rockstar ${basename}\*snap${ranseed}h.\*
py_rockstar ${basename}\*snap${ranseed}i.\*
py_rockstar ${basename}\*snap${ranseed}j.\*
py_rockstar ${basename}\*snap${ranseed}k.\*
py_rockstar ${basename}\*snap${ranseed}l.\*
py_rockstar ${basename}\*snap${ranseed}m.\*
py_rockstar ${basename}\*snap${ranseed}n.\*
py_rockstar ${basename}\*snap${ranseed}o.\*
py_rockstar ${basename}\*snap${ranseed}p.\*
py_rockstar ${basename}\*snap${ranseed}q.\*
py_rockstar ${basename}\*snap${ranseed}r.\*
py_rockstar ${basename}\*snap${ranseed}s.\*
py_rockstar ${basename}\*snap${ranseed}t.\*
py_rockstar ${basename}\*snap${ranseed}u.\*
py_rockstar ${basename}\*snap${ranseed}v.\*
py_rockstar ${basename}\*snap${ranseed}w.\*
py_rockstar ${basename}\*snap${ranseed}x.\*
''')
nowf.close()
