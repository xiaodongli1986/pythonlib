#!/bin/bash
for(( i=0; i<1; ))
do
	echo 'ls and found all .sh files:'
	nowfiles=$(ls *.sh)
	echo $nowfiles
	for nowfile in $nowfiles
	do
		sleep 3
		echo 'now run '$nowfile
		sleep 3
		sh $nowfile
	done
	sleep 3
done
