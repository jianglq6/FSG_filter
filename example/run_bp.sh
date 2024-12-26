#!/bin/bash

#set -x
set -e

date

#-- source intel lib
#source /apps/intel/bin/compilervars.sh intel64
#export FI_PROVIDER=tcp
#MPIDIR=/apps/mpich-3.3

#-- program related dir
CUR_DIR=`pwd`
EXEC_DIR=$CUR_DIR/../
EXEC_WAVE=$EXEC_DIR/main_curv_fsg_2d

#-- project dir for output
PROJDIR=$CUR_DIR/bp_model/fsg_nofilter
PAR_FILE=${PROJDIR}/test.json
GRID_DIR=${PROJDIR}/output
MEDIA_DIR=${PROJDIR}/output
SOURCE_DIR=${PROJDIR}/output
OUTPUT_DIR=${PROJDIR}/output

#-- input dir for input files
IN_MEDIA_LAY_FILE=$CUR_DIR/data/BPTTI_Headpartlg_nx*nz_1016*916.txt

#-- create dir
mkdir -p $PROJDIR
mkdir -p $OUTPUT_DIR
mkdir -p $GRID_DIR
mkdir -p $MEDIA_DIR

#----------------------------------------------------------------------
#-- create main conf
#----------------------------------------------------------------------
cat << ieof > $PAR_FILE
{
  "number_of_total_grid_points_x" : 500,
  "number_of_total_grid_points_z" : 450,

  "size_of_time_step" : 0.006,
  "number_of_time_steps" : 1000,
  "#time_window_length" : 15,
  "#check_stability" : 1,

  "grid_generation_method" : {
      "#import" : "$GRID_DIR",
      "cartesian" : {
        "origin"  : [ 0.0  , 0.0 ],
        "inteval" : [ 60.0,  60.0 ]
      }
  },
  "is_export_grid" : 0,
  "grid_export_dir"   : "$GRID_DIR",

  "medium" : {
      "input_way" : "infile_ascii",
      "input_file" : "${IN_MEDIA_LAY_FILE}"
  },
  "is_export_media" : 1,
  "media_export_dir"  : "$MEDIA_DIR",

  "===is_filter": "0or1",
  "is_filter" : 0,

  "source" : {
     "#abs_4_slides":"",
     "index" : [ 300, 85],
     "ricker_center_frequency" : 2.0,
     "ricker_peak_time" : 0.5,
     "force_vector" : [ 1e16, 1e16]
  },
  "is_export_source" : 1,
  "source_export_dir"  : "$SOURCE_DIR",

  "output_dir" : "$OUTPUT_DIR",

  "check_nan_every_nummber_of_steps" : 0,
  "output_all" : 0
}
ieof

echo "+ created $PAR_FILE"

#-------------------------------------------------------------------------------
#-- Performce simulation
#-------------------------------------------------------------------------------
#

#-- gen run script
cat << ieof > ${PROJDIR}/cgfd2d_wave_sim.sh
#!/bin/bash

#-- simulation
printf "\nStart simualtion ...\n";

time $EXEC_WAVE $PAR_FILE 100
if [ $? -ne 0 ]; then
    printf "\nSimulation fail! stop!\n"
    exit 1
fi

ieof

#-------------------------------------------------------------------------------
#-- start run
#-------------------------------------------------------------------------------

chmod 755 ${PROJDIR}/cgfd2d_wave_sim.sh
${PROJDIR}/cgfd2d_wave_sim.sh
if [ $? -ne 0 ]; then
    printf "\nSimulation fail! stop!\n"
    exit 1
fi

date

# vim:ft=conf:ts=4:sw=4:nu:et:ai:
