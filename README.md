
A collection of python routines, and a powerful programme py_Plot.

To use it, firstly download and install anaconda

	http://continuum.io/downloads

Then add something like this to ~/.bashrc:

	export pythonlibPATH=/home/xiaodongli/software/pythonlib
	export PYTHONPATH=${pythonlibPATH}:${PYTHONPATH}
	export PATH=${pythonlibPATH}/bin:${PATH}

To use the libraries, type this in your source

	import stdA as stdA
		or you can type "execfile('/home/xiaodongli/software/pythonlib/stdA.py')"

To install py_Plot
	
	cd src
	Modify the first line of Makefile "pythonbin = \#!/home/xiaodongli/software/anaconda/bin/python"
	make

To see how to use py_Plot

	py_Plot

Here is an example
	
	py_Plot scatter3d \*.txt -xcol 1 -ycol 2 -zcol 3 -randrat 0.1 -savefig T -figfmt png -showfig F 
		3d scatter plot 1th, 2th, 3th columns of all files end with .txt
		randomly select 10% of the file and plot
		all plottings saved as png files, no display on the screen


Write your own EXE

	cd src
	./Add_EXE YOUREXENAME

		A file named 'src/py_YOUREXENAME.py' will be created and its information will be added to Makefile.
		Edit that .py file, then make
