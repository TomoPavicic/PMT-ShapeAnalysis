"""
This file is an example file to show the user how to use the PMT class

Without uncommenting this file reads sample data and produces two files:
    - test_output_root_file.root
    - test_output_root_file_array_version.root

William Quinn 2020
"""

import sys
sys.path.insert(1, '..')

from scr.PMT_Array import PMT_Array
from scr.PMT_Object import PMT_Object
from scr.PMT_Waveform import PMT_Waveform


def read_test_file():
    try:
        pmt_data_file = open("./test_data.dat", 'r')
    except FileNotFoundError as fnf_error:
        print(fnf_error)
        raise Exception("Error opening data file.")
    data_list = []
    for line_index, line in enumerate(pmt_data_file.readlines()):
        waveform = []
        for token_index, token in enumerate(line.split()):
            waveform.append(token.strip())
        data_list.append(waveform)
    return data_list


def main():
    # Create an arbitrary experiment with a single PMT
    # The topology can then just be [1,1]
    topology = [1, 1]

    # Create a PMT_Array which automatically creates a series of PMT_Objects
    # equal to the pmts in your grid
    pmt_array = PMT_Array(topology, "test")

    # Set the cuts from a config file. This also sets up the histograms in the array format
    pmt_array.apply_setting("example_config_file.txt")

    # Or you can create a PMT_Oject directly. The Array class is a simple object
    # to contain all the pmts and has a few functions to make things easier
    pmt_object = PMT_Object("PMT_1", "test")

    # Set up the histograms using the default settings
    pmt_object.set_up_histograms()

    # Now you need data. I have supplied a sample of waveforms in ./Data/run_214
    # Real data will be in various different formats
    # For now I will use a function that I have created for this to output the data
    # into the format we want.... a list of pmt data points

    # Ignore this function as it is a simple read ascii function.
    # The output is a 2D list with each [i][0:1024] being a different waveform
    data_list = read_test_file()

    # Loop over this data
    for i_waveform in range(len(data_list)):

        # Put each waveform into its own dedicated container
        pmt_waveform = PMT_Waveform(data_list[i_waveform], pmt_object)

        # or using the array, note the pmt_number convention - 1 PMT -> PMT_number = 0
        pmt_waveform_array_version = PMT_Waveform(data_list[i_waveform], pmt_array.get_pmt_object_number(0))

        # From this point you can do whatever you like. All the info is stored in pmt_waveform
        # Uncomment and try the following

        # Save the waveform to a file
        '''test_root_file_output_1 = ROOT.TFile(pmt_object.get_pmt_id() + "_" + pmt_waveform.get_pmt_trace_id() + ".root", "RECREATE")
        pmt_waveform.save_pmt_waveform_histogram(test_root_file_output_1)
        test_root_file_output_1.Close()
        del test_root_file_output_1'''

        # Store the waveform results into the PMT object
        # Although there are cuts in the PMT classes I'd recommend at this point applying additional cuts
        # to filter your data
        ######################
        #      APPLY CUTS    #
        ######################
        # Fill results. You don't need to input any info to tell what to fill where as all this info is
        # handled in the classes
        # Only default hists have bee created. Thoroughly read the class to see how to change from defaults
        pmt_waveform.fill_pmt_hists()
        pmt_waveform_array_version.fill_pmt_hists()

        # Note here that if you store the waveform into memory delete it
        del pmt_waveform
        del pmt_waveform_array_version

    # Once you are done with looping over the data you can look at the results
    pmt_object.save_to_file("test_output_root_file.root")
    pmt_array.save_to_file("test_output_root_file_array_version.root")
    print(">>> Success")


if __name__ == '__main__':
    main()
