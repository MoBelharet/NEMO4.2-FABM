
configuration=3D
compile_with_fabm=false
include_top=true
offline_mode=false

config_name=nemo_${configuration}
oce_lab=OCE
#-------------- INCLUDE TOP ------------------------
if [ $include_top = true ]; then
	top_lab1=TOP
	top_lab2=key_top
	config_name="${config_name}_pisces"

	if [ $compile_with_fabm = true ] ; then
		config_name="${config_name}_fabm"
	fi
else
	top_lab1=""
	top_lab2=""
fi


#-------------- OFFLINE MODE -----------------------
if [ $offline_mode = true ]; then
	#oce_lab=""
	off_lab1="OFF"
	off_lab2=key_offline
	config_name="${config_name}_offline"
else
	off_lab1=""
	off_lab2=""
fi
#---------------------------------------------------

#PRJD=/home/ext/mr/smer/belharetm/PROJECT_NEMO_FABM_PISCES
#FABM_INC=$PRJD/local/fabm/nemo/include ; export FABM_INC
#FABM_LIB=$PRJD/local/fabm/nemo/lib64 ; export FABM_LIB

rm -f cfg/$config_name/cpp_$config_name.fcm

if [ $compile_with_fabm = true ] ; then
     if [ $configuration = 1D ]; then
    	./makenemo -d "$oce_lab $top_lab1 $off_lab1"  -m belenos_MERCATOR_fabm -r C1D_PAPA -n $config_name -j 12 -k 0 add_key "$top_lab2 key_fabm" $1
     elif [ $configuration = 3D ]; then
        ./makenemo -m belenos_MERCATOR_fabm -r ORCA2_ICE_PISCES -n $config_name -j 12 -k 0 add_key "key_fabm" $1
     fi
else
  #  ./makenemo -d "$oce_lab $top_lab1 $off_lab1"  -m belenos_MERCATOR -r C1D_PAPA -n $config_name  add_key "$top_lab2" $1
    if [ $configuration = 1D ]; then
    	./makenemo -d "OCE TOP"  -m belenos_MERCATOR -r C1D_PAPA -n $config_name -j 12 -k 0 add_key "key_nosignedzero key_top key_xios" $1
    elif [ $configuration = 3D ]; then
        ./makenemo -m belenos_MERCATOR -r ORCA2_ICE_PISCES -n $config_name -j 12 -k 0  $1
    fi

    
fi

#key_nosignedzero

