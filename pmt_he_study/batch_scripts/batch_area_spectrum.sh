#!/usr/bin/env bash

# ------------------------------------------------------------------------------------ #
# Dave Waters, 24.01.20
#
# This is a script for launching batch jobs for the analysis of PMT data.
# You should first test the macro ("root_macro" defined below) interactively and then
# set the parameters here to run in batch mode.
#
# It is assumed that the root macro itself takes care of storing the output in a
# suitable location and format.
#
# History
# =======
# [1] 24.01.20 : initial version tested and submitted to GitHub. Not yet tested on a
#     large scale, so problems cannot be excluded. For example what happens if jobs
#     running concurrently clash when appending results to a common output file ?
#
# USER INPUTS:
# ============
# Descriptive tag (will be used for naming of directory, scripts & log files):
# tag="AS_26.03.20"
while getopts "t:V:" option;
do
  case "${option}" in
      t )
        tag=${OPTARG}
        ;;
      V )
        voltage=${OPTARG}
        ;;
  esac
done
score="_"
volt="V"

file_filenames="filenames_$voltage.ascii"

root_dir="/unix/nemo4/PMT_He_Study_nemo4"
script_log_dir=batch_log_$tag
mkdir $root_dir/PMT-ShapeAnalysis/pmt_he_study/batch_files/$script_log_dir
mkdir $root_dir/PMT-ShapeAnalysis/pmt_he_study/batch_files/$script_log_dir/scripts
mkdir $root_dir/PMT-ShapeAnalysis/pmt_he_study/batch_files/$script_log_dir/logs

ls $root_dir/data/raw_xml_files/*/*$voltage* /unix/nemo3/PMT_He_Study_nemo3/data/raw_xml_files/*/*$voltage* > $root_dir/PMT-ShapeAnalysis/pmt_he_study/batch_files/$script_log_dir/$file_filenames
data_file_list=`cat $root_dir/PMT-ShapeAnalysis/pmt_he_study/batch_files/$script_log_dir/$file_filenames`
#
# ROOT macro to use:
root_macro="$root_dir/PMT-ShapeAnalysis/pmt_he_study/pmt_charge_analysis.cpp"

echo $voltage
echo ""
echo $tag
echo ""
#
# Batch queue to use (almost certainly the short queue):
queue="short"
# ------------------------------------------------------------------------------------ #

ifile=0
for file in $data_file_list;
do
  echo "Processing file: $file"
  # Extract the date from the file/directory:
  date=`echo $file | awk -F '/' '{print $7}'`
  xml_file=${file##*/}
  base=${xml_file%.*}
  temp=`echo $base | awk -F '_t' '{print $2}'`
  time=`echo $temp | awk -F '_' '{print $1}'`
  root_file="$root_dir/ROOT_files/Area_Spectrum_Files/$voltage$volt/$date$score$time$score$base.root"

  if test -f "$root_file"; then
    echo "$root_file already exists"
  else
    # Prepare command in stages:
    macro_plus_arguments="\"$root_macro(\\\"$file\\\",\\\"$date\\\",$voltage)\""
    echo $macro_plus_arguments
    #
    interactive_command="root.exe -q -b $macro_plus_arguments >& $root_dir/PMT-ShapeAnalysis/pmt_he_study/batch_files/$script_log_dir/logs/batch_file_$ifile.log"
    echo $interactive_command
    # $interactive_command
    #
    # Prepare batch script in a temporary directory (not sure if this is robust)
    touch $root_dir/PMT-ShapeAnalysis/pmt_he_study/batch_files/$script_log_dir/scripts/batch_file_$ifile.sh
    chmod a+rwx $root_dir/PMT-ShapeAnalysis/pmt_he_study/batch_files/$script_log_dir/scripts/batch_file_$ifile.sh
    echo "echo 'Running batch_file_$ifile.sh ...' >> $root_dir/PMT-ShapeAnalysis/pmt_he_study/batch_files/$script_log_dir/logs/batch_file_$ifile.log" > $root_dir/PMT-ShapeAnalysis/pmt_he_study/batch_files/$script_log_dir/scripts/batch_file_$ifile.sh
    echo $interactive_command >> $root_dir/PMT-ShapeAnalysis/pmt_he_study/batch_files/$script_log_dir/scripts/batch_file_$ifile.sh
    batch_command="qsub -q $queue $root_dir/PMT-ShapeAnalysis/pmt_he_study/batch_files/$script_log_dir/scripts/batch_file_$ifile.sh"
    echo $batch_command
    # For testing (without submitting to the batch farm), comment out the next line:
    $batch_command
  fi

  echo ""
  let ifile=$ifile+1
done
