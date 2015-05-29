
Here I have a collection of python routines, and a powerful, convenient programme py_Plot 

To install, first install anaconda

	http://continuum.io/downloads

Then add something like this to your bashrc:

	export pythonlibPATH=/home/xiaodongli/software/pythonlib
	export PYTHONPATH=${pythonlibPATH}:${PYTHONPATH}
	export PATH=${pythonlibPATH}/bin:${PATH}

To use the libraries, type this in your source

	import stdA as stdA
	execfile('/home/xiaodongli/software/pythonlib/stdA.py')

To install py_Plot
	
	cd src
	Modify the first line of Makefile "pythonbin = \#!/home/xiaodongli/software/anaconda/bin/python"
	make

Add your own EXE to src:

	cd src
	./Add_EXE YOURBINNAME

		A file named 'src/py_YOURBINNAME.py' will be created and its information will be added to Makefile. Edit that .py file, then you can call it by shell command 'py_YOURBINNAME'
