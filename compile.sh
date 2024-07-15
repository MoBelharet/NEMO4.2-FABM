
compile_with_fabm=true
include_top=true
offline_mode=false

config_name=test_nemo
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

FABM_INC=/home/mbelharet/FABM_PROJECT/PROJECT_FABM_NEMO4.2_PISCES/local/fabm/nemo/include ; export FABM_INC
FABM_LIB=/home/mbelharet/FABM_PROJECT/PROJECT_FABM_NEMO4.2_PISCES/local/fabm/nemo/lib ; export FABM_LIB

rm -f cfg/$config_name/cpp_$config_name.fcm

if [ $compile_with_fabm = true ] ; then
    #./makenemo -m X64_PX_GCC_OMPI_fabm -r C1D_PAPA_FABM -n test_nemo_fabm $1

    ./makenemo -d "$oce_lab $top_lab1 $off_lab1"  -m X64_PX_GCC_OMPI_fabm -r C1D_PAPA -n $config_name -j 12 -k 0 add_key "$top_lab2 key_fabm" $1
    #rm -f /home/mbelharet/FABM_PROJECT/PROJECT_FABM_NEMO4.2_PISCES/run_nemo_fabm_pisces/nemo.exe
    #ln cfgs/test_nemo_fabm/BLD/bin/nemo.exe  /home/mbelharet/FABM_PROJECT/PROJECT_FABM_NEMO4.2_PISCES/run_nemo_fabm_pisces/
else
    ./makenemo -d "$oce_lab $top_lab1 $off_lab1"  -m X64_PX_GCC_OMPI -r C1D_PAPA -n $config_name -j 12 -k 0 add_key "$top_lab2" $1
    #rm -f /home/mbelharet/FABM_PROJECT/PROJECT_FABM_NEMO4.2_PISCES/run_nemo_pisces/nemo.exe
    #ln cfgs/$config_name/BLD/bin/nemo.exe /home/mbelharet/FABM_PROJECT/PROJECT_FABM_NEMO4.2_PISCES/run_nemo_pisces/
fi

#key_nosignedzero

