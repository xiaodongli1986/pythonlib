pythonbin = \#!/home/ubuntu/archiconda3/bin/python
python3bin = \#!/home/ubuntu/archiconda3/bin/python
#EXEs = ../bin/py_helloworld  ../bin/py_Plot  ../bin/uws_download ../bin/uws_id  ../bin/uws_show  ../bin/uws_submit ../bin/uws_mergefiles ../bin/checkxyzwfiles ../bin/py_ran ../bin/uws_checkrun ../bin/py_stat ../bin/py_MIAOEXE ../bin/py_CUTE
EXEs =   ../bin/py_Plot  ../bin/uws_download ../bin/uws_id  ../bin/uws_show  ../bin/uws_submit ../bin/uws_mergefiles ../bin/checkxyzwfiles ../bin/py_ran ../bin/uws_checkrun ../bin/py_stat ../bin/py_MIAOEXE ../bin/py_CUTE ../bin/py_ran_shell ../bin/py_sort_file_acto_size ../bin/py_merge_files ../bin/py_LightconeHalo_halo_for_iboxes ../bin/py_LightconeHalo_halo_for_iboxes_merge ../bin/py_xi ../bin/py_binary_ran ../bin/py_binary_ascii_conv  ../bin/exec_all_sh

EXEpy3s = ../bin_py3/py_cic ../bin_py3/py_rockstar_to_fmt3 ../bin_py3/py_lpicola ../bin_py3/py_rockstar ../bin_py3/py_lpicola_postprocess ../bin_py3/py_add1 ../bin_py3/py_cola_lc  ../bin_py3/py_boss_split_info  ../bin_py3/py_data_mock  ../bin_py3/py_rockstar_multisnap   ../bin_py3/py_selenium  ../bin_py3/py_fake_headinfo
#../bin_py3/py_cola

all: $(EXEs) $(EXEpy3s)
default: $(EXEs) $(EXEpy3s)

../bin/%: %.py
	#rm  ../bin/$*
	echo '${pythonbin}' >> ../bin/$*
	cat  $*.py >> ../bin/$*
	chmod a+x ../bin/$*
../bin_py3/%: %.py
	#rm  ../bin_py3/$*
	echo '${python3bin}' >> ../bin_py3/$*
	cat  $*.py >> ../bin_py3/$*
	chmod a+x ../bin_py3/$*
	cp ../bin_py3/$* ../bin/

../bin/%: %.sh
	cp $*.sh ../bin/$*.sh
	chmod a+x ../bin/$*.sh
		
clean :
	rm -rf ${EXEs} ${EXEpy3s}
	
