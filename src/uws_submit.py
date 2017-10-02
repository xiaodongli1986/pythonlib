
#for ((snapnum=80;snapnum>=80;snapnum--))
#do
#for ((iz=1;iz<=10;iz++))
for iz in 7 8 9
#for ((iz=11;iz<=250;iz++))
do
#continue
z2=$(($iz*10))
z1=$(($z2-10))
echo $iz $z1 $z2
uws $cstr job new queue="long" query="SELECT * from BigMDPL.Rockstar_z0 where z<${z2} and z>${z1}" table="BigMDPL_z${z1}to${z2}" --run
done
#done

#uws $cstr job new queue="long" query="SELECT DISTINCT * FROM BigMDPL.Redshifts ORDER BY snapnum DESC" table="BigMDPL_Redshifts" --run
#uws $cstr job new queue="long" query="SELECT * from BigMDPL.Rockstar where snapnum=79  and z<1 and z>0" table="BigMDPL_test3" --run
#uws $cstr job new queue="long" query="SELECT * from BigMDPL.Rockstar_z0 where z<1 and z>0" table="BigMDPL_test4" --run

#http --auth xiaodongli:lx821118 --form --follow POST https://www.cosmosim.org/uws/query query="SELECT DISTINCT * FROM BigMDPL.Redshifts ORDER BY snapnum DESC" table="BigMDPL_reds"
