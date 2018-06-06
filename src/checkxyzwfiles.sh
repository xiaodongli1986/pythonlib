#echo 'Please type the files you want to check & plot (e.g., \*.xyzw to check all xyzw files) ...'
#read files
#echo \*.xyzw
echo 'Please type your randrat...'
read randrat
echo $randrat
py_Plot scatter3d \*.xyzw  -randrat $randrat
py_Plot scatter \*.xyzw -xcol 1 -ycol 2 -randrat $randrat
py_Plot scatter \*.xyzw -xcol 1 -ycol 3 -randrat $randrat
py_Plot scatter \*.xyzw -xcol 2 -ycol 3 -randrat $randrat
py_Plot hist \*.xyzw -xcol 1
py_Plot hist \*.xyzw -xcol 2
py_Plot hist \*.xyzw -xcol 3
py_Plot hist \*.xyzw -xcol 4
