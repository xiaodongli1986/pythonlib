pythonbin = \#!/home/xiaodongli/software/anaconda/bin/python
EXEs = ../bin/py_helloworld  ../bin/py_Plot  ../bin/uws_download ../bin/uws_id  ../bin/uws_show  ../bin/uws_submit

all: $(EXEs)
default: $(EXEs)

../bin/%: %.py
	rm -rf ../bin/$*
	echo '${pythonbin}' >> ../bin/$*
	cat  $*.py >> ../bin/$*
	chmod a+x ../bin/$*
		
clean :
	rm -rf ${EXEs}
	
