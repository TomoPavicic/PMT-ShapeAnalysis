"""
This is the python file that runs over the He permeation .xml files
By default it does not do any matched filtering analysis unless you tell it to do so

Run this file at the terminal by:
    - python3 analysis.py -i $INPUT_DATA_FILE_PATH$ [-c $INPUT_CONFIG_FILE$] [-sweep True/False]
                                                [-f True/False] [-r True/False]
where [] arguments are optional, the defaults will be used.
The default settings can be changed in PMT_Object.py in the _inti_() function
"""
import sys
sys.path.insert(1, '..')

from functions.other_functions import parse_arguments, get_date_time, fit_bismuth_function_from_file
from functions.data_reader_functions import process_xml_file
from scr.PMT_Array import PMT_Array


def main():

    # Handle the input arguments:
    ##############################
    args = parse_arguments()
    input_data_file_name = args.i
    config_file_name = args.c
    sweep_bool = args.sweep
    recreate_bool = args.r
    bismuth_bool = args.f
    ##############################

    # Do some string manipulation to get the date and time from the file name
    #################################################
    date, time = get_date_time(input_data_file_name)
    #################################################

    # Check to see if the data file exists
    try:
        print(">>> Reading data from file: {}".format(input_data_file_name))
        date_file = open(input_data_file_name, 'r')
    except FileNotFoundError as fnf_error:
        print(fnf_error)
        raise Exception("Error opening data file {}".format(input_data_file_name))

    # Define the name of the output root file you want to create
    ##########################################################
    temp = input_data_file_name.split("/")
    temp1 = temp[-1].split(".")
    output_file_name = date + "_" + temp1[0] + "_output.root"
    ##########################################################

    if recreate_bool == "False" or recreate_bool is None:
        fit_bismuth_function_from_file(output_file_name)
    elif recreate_bool == "True":

        # Create the object to store all the pmt information
        topology = [2, 1]
        pmt_array = PMT_Array(topology, date + "_" + time)
        pmt_array.set_pmt_id("GAO607_" + date + "_" + time, 0)
        pmt_array.set_pmt_id("GAO612_" + date + "_" + time, 1)

        # Set the cuts you wish to apply
        # If you don't do this the defaults are used
        if config_file_name is not None:
            pmt_array.apply_setting(config_file_name)

        if sweep_bool == "True":
            pmt_array.set_sweep_bool(True)
            pmt_array.set_pmt_templates("/Users/willquinn/Documents/PhD/PMT_Permeation_Project/data/21.06.19_A1400_B1400_t1130_templates.root", "Template_Waveform_Channel0_A1400_B1400_t1130")

        try:
            root_file = open(output_file_name, 'r')
            print(">>> Output file already exists: {}".format(output_file_name))
            print(">>> Recreating...")
            root_file.close()
            del root_file
        except:
            print(">>> Output file doesn't exist, creating {}".format(output_file_name))

        # To clean up the code I put the xml parsing code into a separate function
        # Pass this function your pmt_array and it will fill the data directly into it
        process_xml_file(input_data_file_name, pmt_array)

        # This is another way of processing the file but is slower
        '''waveform_list = process_xml_file_new(input_data_file_name)
        # Then loop over this list
        event_counter = 0
        percentage_counter = 1
        print(">>> Looping over events")
        temp_start = TIME.time()
        for i_waveform in range(len(waveform_list)):

            if event_counter == int(len(waveform_list)/10):
                event_counter = 0
                intermediate = TIME.time()
                time_length = intermediate - temp_start
                print(">>>\n>>>  %.3f s.\n" % (intermediate - temp_start))
                temp_start = intermediate
                print("Processed {}% of data...".format(percentage_counter * 10))
                estimate = time_length * (10 - percentage_counter)
                print(">>> Estimated time till termination %.3f seconds\n\n" % estimate)
                percentage_counter += 1

            pmt_waveform = PMT_Waveform(waveform_list[i_waveform][1], pmt_array.get_pmt_object_number(waveform_list[i_waveform][0]))

            if pmt_waveform.get_pulse_trigger():
                pmt_waveform.fill_pmt_hists()

            del pmt_waveform
            event_counter += 1'''

        # This is not OOP, it is very specific to the PME Permeation project
        if bismuth_bool == "True":
            pmt_array.fit_bismuth_function()

        # Save the results to a file
        pmt_array.save_to_file(output_file_name)


if __name__ == '__main__':
    main()
