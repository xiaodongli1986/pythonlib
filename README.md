
1. What are things here

	a. A collection of python routines 
	b. A powerful, convenient programme py_Plot 

2. How to install


Install anaconda to have all libraries needed

	http://continuum.io/downloads

Then add something like this to your bashrc:

	export pythonlibPATH=/home/xiaodongli/software/pythonlib
	export PYTHONPATH=${pythonlibPATH}:${PYTHONPATH}
	export PATH=${pythonlibPATH}/bin:${PATH}

###############################
 python library
###############################

The most useful one is stdA.
I also wrote something for 2pcf and boss data.

To use them

	execfile('/home/xiaodongli/software/pythonlib/stdA.py')

Or 
	import stdA as stdA


###############################
 python bin
###############################

Python EXEs which can be used in terminal.

To compile 
	
	cd src
	Modify the first line of Makefile "pythonbin = \#!/home/xiaodongli/software/anaconda/bin/python"
	make

Add your own EXE:

	cd src
	./Add_EXE YOURBINNAME

A file named 'src/py_YOURBINNAME.py' will be created and its information will be added to Makefile. Edit that .py file, then you can call it by shell command 'py_YOURBINNAME'
