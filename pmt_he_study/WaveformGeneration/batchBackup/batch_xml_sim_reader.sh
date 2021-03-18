#!/usr/bin/env bash

#Launches batch job for analysing simulated data
#Essentially want to repeat the command './xml_reader --i input_file --c config_file --o output_file


#tag = RandomApPlacementsSim
while getopts :t: option;
do
  case "${option}" in
      t )
        tag=${OPTARG}
        ;;
  esac
done

root_dir="/home/tomopavicic/PMT_Helium_WD/PMT-ShapeAnalysis/pmt_he_study/WaveformGeneration/$tag"
#script_log_dir=batch_log_$tag
mkdir $root_dir/batch_files
mkdir $root_dir/batch_files/scripts
mkdir $root_dir/batch_files/logs

ls $root_dir/raw_xml/*.xml > $root_dir/batch_files/filenames.ascii

data_file_list=`cat $root_dir/batch_files/filenames.ascii`

config_file="$root_dir/ATC40Copy.conf"

# ROOT macro to use:
source_executable="/home/tomopavicic/PMT_Helium_WD/PMT-ShapeAnalysis/pmt_he_study/build/xml_reader"
#
underscore="_"
dates="110011"
# Batch queue to use:
queue="short"

#
ifile=0
for file in $data_file_list;
do
  #date=`echo $file | awk -F '/' '{print $7}'`
  xml_file=${file##*/}
  base=${xml_file%.*}
  root_file="$root_dir/RootFiles/$dates$underscore$base.root"
  if test -f "$root_file"; then
    echo "$root_file already exists"
  else
    echo "$root_file does not exists"
    echo "Processing file: $file"
    # Prepare command in stages:
    executable_plus_arguments="$source_executable --i $file --c $config_file --o $root_file"
    #
    interactive_command="$executable_plus_arguments >& $root_dir/batch_files/logs/batch_file_$ifile.log"
    echo $interactive_command
    #
    # Prepare batch script in a temporary directory (not sure if this is robust)
    touch $root_dir/batch_files/scripts/batch_file_$ifile.sh
    chmod a+rwx $root_dir/batch_files/scripts/batch_file_$ifile.sh
    echo "source /home/tomopavicic/PMT_Helium_WD/PMT-ShapeAnalysis/pmt_he_study/setup.sh" >> $root_dir/batch_files/scripts/batch_file_$ifile.sh
    echo "echo 'Running batch_file_$ifile.sh ...' >> $root_dir/batch_files/logs/batch_file_$ifile.log" >> $root_dir/batch_files/scripts/batch_file_$ifile.sh
    echo $interactive_command >> $root_dir/batch_files/scripts/batch_file_$ifile.sh
    batch_command="qsub -q $queue $root_dir/batch_files/scripts/batch_file_$ifile.sh"
    echo $batch_command
    # For testing (without submitting to the batch farm), comment out the next line:
    $batch_command
  fi

  echo ""
  let ifile=$ifile+1
done