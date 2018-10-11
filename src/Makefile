pythonbin = \#!/home/miao/anaconda2/bin/python
EXEs = ../bin/py_helloworld  ../bin/py_Plot  ../bin/uws_download ../bin/uws_id  ../bin/uws_show  ../bin/uws_submit ../bin/uws_mergefiles ../bin/checkxyzwfiles ../bin/py_ran ../bin/uws_checkrun ../bin/py_stat ../bin/py_MIAOEXE

all: $(EXEs)
default: $(EXEs)

../bin/%: %.py
	rm -rf ../bin/$*
	echo '${pythonbin}' >> ../bin/$*
	cat  $*.py >> ../bin/$*
	chmod a+x ../bin/$*

../bin/%: %.sh
	cp $*.sh ../bin/$*.sh
	chmod a+x ../bin/$*.sh
		
clean :
	rm -rf ${EXEs}
	
