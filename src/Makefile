pythonbin = \#!/home/xiaodongli/software/anaconda2/bin/python
#EXEs = ../bin/py_helloworld  ../bin/py_Plot  ../bin/uws_download ../bin/uws_id  ../bin/uws_show  ../bin/uws_submit ../bin/uws_mergefiles ../bin/checkxyzwfiles ../bin/py_ran ../bin/uws_checkrun ../bin/py_stat ../bin/py_MIAOEXE ../bin/py_CUTE
EXEs =   ../bin/py_Plot  ../bin/uws_download ../bin/uws_id  ../bin/uws_show  ../bin/uws_submit ../bin/uws_mergefiles ../bin/checkxyzwfiles ../bin/py_ran ../bin/uws_checkrun ../bin/py_stat ../bin/py_MIAOEXE ../bin/py_CUTE ../bin/py_ran_shell ../bin/py_sort_file_acto_size ../bin/py_merge_files ../bin/py_LightconeHalo_halo_for_iboxes ../bin/py_LightconeHalo_halo_for_iboxes_merge

all: $(EXEs)
default: $(EXEs)

../bin/%: %.py
	/usr/bin/rm -rf ../bin/$*
	echo '${pythonbin}' >> ../bin/$*
	cat  $*.py >> ../bin/$*
	chmod a+x ../bin/$*

../bin/%: %.sh
	cp $*.sh ../bin/$*.sh
	chmod a+x ../bin/$*.sh
		
clean :
	/usr/bin/rm -rf ${EXEs}
	
