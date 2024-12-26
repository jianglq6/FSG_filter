#!/bin/bash

#set -x
set -e

date

#-- source intel lib
#source /apps/intel-2019.3/intel/bin/compilervars.sh intel64
#export FI_PROVIDER=tcp
#MPIDIR= /apps/intel-2019.3/mpich-3.3

#-- program related dir
CUR_DIR=`pwd`
EXEC_DIR=$CUR_DIR/../
EXEC_WAVE=$EXEC_DIR/main_curv_fsg_2d

#-- project dir for output
PROJDIR=$CUR_DIR/layer_model/fsg_filter_tti
PAR_FILE=${PROJDIR}/test.json
GRID_DIR=${PROJDIR}/output
MEDIA_DIR=${PROJDIR}/output
SOURCE_DIR=${PROJDIR}/output
OUTPUT_DIR=${PROJDIR}/output

#-- input dir for input files
IN_MEDIA_LAY_FILE=$CUR_DIR/data/test_2aniso_layer.md2lay

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
  "number_of_total_grid_points_z" : 500,

  "size_of_time_step" : 0.01,
  "number_of_time_steps" : 700,
  "#time_window_length" : 15,
  "check_stability" : 1,

  "grid_generation_method" : {
      "#import" : "$GRID_DIR",
      "cartesian" : {
        "origin"  : [ 0.0  ,  0.0 ],
        "inteval" : [ 100.0,  100.0 ]
      }
  },
  "is_export_grid" : 0,
  "grid_export_dir"   : "$GRID_DIR",

  "#===if_filtering===#": "0 or 1",
  "is_filter" : 1,

  "medium" : {
      "input_way" : "infile_layer",
      "input_file" : "${IN_MEDIA_LAY_FILE}",
      "#==if_using_local_points==" : "loc",
      "#==if_sm_calculus==" : "tti",
      "equivalent_medium_method" : "tti"
  },
  "is_export_media" : 1,
  "media_export_dir"  : "$MEDIA_DIR",

  "source" : {
         "index" : [ 250, 275 ],
         "ricker_center_frequency" : 2.0,
         "ricker_peak_time" : 0.5,
         "force_vector" : [ 1e16, 1e16]
   },

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
