from scr.PMT_Array import PMT_Array
from scr.PMT_Waveform import PMT_Waveform
from xml.dom import minidom
import time as TIME
import ROOT


def process_xml_file(input_data_file_name: str, pmt_array: PMT_Array):
    print(">>> Parsing the data file...")
    processing_start = TIME.time()

    # parse an xml file by name
    xml_file = minidom.parse(input_data_file_name)
    events = xml_file.getElementsByTagName('event')
    parse_time = TIME.time() - processing_start
    print(">>> File is good. Parse time: %.3f s" % parse_time)
    print(">>> Number of Events: {}".format(len(events)))

    event_counter = 0
    percentage_counter = 1
    print(">>> Looping over events")
    temp_start = TIME.time()
    for event_index, event in enumerate(events):

        if event_counter == int(len(events)/20):
            event_counter = 0
            intermediate = TIME.time()
            time_length = intermediate - temp_start
            print(">>>\n>>>  %.3f s.\n" % (intermediate - temp_start))
            temp_start = intermediate
            print("Processed {}% of data...".format(percentage_counter * 5))
            estimate = time_length * (20 - percentage_counter)
            print(">>> Estimated time till termination %.3f seconds\n\n" % estimate)
            percentage_counter += 1

        traces = event.getElementsByTagName('trace')
        for trace_index, trace in enumerate(traces):
            # Channel refers to the pmt number
            channel_id = int(trace.attributes['channel'].value)

            # Important Code:
            # This is where you pass the data to the OOP code which does all the analysis for you
            # Ideally this would be in the main analysis.py file
            ##########################################################################################################
            pmt_waveform = PMT_Waveform(trace.firstChild.data.split(" "), pmt_array.get_pmt_object_number(channel_id))
            # check waveform to see if you want to fill histograms
            # The pulse trigger logic is in PMT_Waveform.py
            if pmt_waveform.get_pulse_trigger():
                pmt_waveform.fill_pmt_hists()
            else:
                pass
            del pmt_waveform
            ##########################################################################################################

        event_counter += 1


def process_xml_file_new(input_data_file_name: str):
    print(">>> Parsing the data file...")
    processing_start = TIME.time()

    # parse an xml file by name
    xml_file = minidom.parse(input_data_file_name)
    events = xml_file.getElementsByTagName('event')
    parse_time = TIME.time() - processing_start
    intermediate_time = TIME
    print(">>> File is good. Parse time: %.3f s" % parse_time)
    print(">>> Number of Events: {}".format(len(events)))

    waveform_data_list = []
    for event_index, event in enumerate(events):

        traces = event.getElementsByTagName('trace')
        for trace_index, trace in enumerate(traces):
            waveform_data_list.append([int(trace.attributes['channel'].value), trace.firstChild.data.split(" ")])
    return waveform_data_list


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


def process_crd_file(input_data_file_name: str, pmt_array: PMT_Array, waveform_output_file: ROOT.TFile):

    try:
        pmt_data_file = open(input_data_file_name, 'r')
    except FileNotFoundError as fnf_error:
        print(fnf_error)
        raise Exception("Error opening data file. Skip to the next file...")

    new_waveform_bool = False
    line_number_int = 0

    # print(pmt_array.get_pmt_object_number(0).get_histogram_dict().keys())

    for pmt_data_index, pmt_data_line in enumerate(pmt_data_file.readlines()[10:]):
        pmt_data_line_tokens = pmt_data_line.split(" ")

        if pmt_data_line_tokens[0] == "=" and pmt_data_line_tokens[1] == "HIT":
            new_waveform_bool = True
            line_number_int = 0
        else:
            pass

        if new_waveform_bool and line_number_int == 1:
            pmt_slot_number = int(pmt_data_line_tokens[1]) # Column
            pmt_channel_number = int(pmt_data_line_tokens[3]) # Row

            pmt_number = int(pmt_slot_number + pmt_array.get_pmt_topology()[1]*pmt_channel_number)
            '''if pmt_slot_number == 0:
                print(pmt_slot_number, pmt_channel_number)
                print(pmt_number)'''
            event_id_LTO = int(pmt_data_line_tokens[5])
            event_id_HT = int(pmt_data_line_tokens[7])
            pmt_waveform_peak_cell = int(pmt_data_line_tokens[27])
            pmt_waveform_charge = float(pmt_data_line_tokens[29])
            pmt_waveform_rise_time = float(pmt_data_line_tokens[39])
            if int(event_id_HT) != 0:
                pass
            elif int(event_id_LTO) != 0:
                pass
            else:
                pass

        elif new_waveform_bool and line_number_int == 2:
            if event_id_HT != 0:

                pmt_adc_values = []
                for i_adc in range(len(pmt_data_line_tokens)):
                    pmt_adc_values.append(pmt_data_line_tokens[i_adc])

                pmt_waveform = PMT_Waveform(pmt_adc_values, pmt_array.get_pmt_object_position([pmt_channel_number, pmt_slot_number]))
                if pmt_waveform.get_pulse_trigger():
                    '''print("results: ", pmt_waveform.get_results_dict())
                    print("")'''
                    '''print("MF Amplitude: ",pmt_waveform.get_pmt_pulse_mf_amp())
                    print("MF Shape: ", pmt_waveform.get_pmt_pulse_mf_shape())
                    print("Amplitude: ", pmt_waveform.get_pmt_pulse_peak_amplitude())
                    print("")'''
                    #print("Slot: ", pmt_channel_number, "Channel: ", pmt_channel_number)
                    pmt_waveform.fill_pmt_hists()

                if pmt_waveform.get_pmt_apulse_trigger():
                    print("pre_pulse")
                    pmt_waveform.save_pmt_waveform_histogram(waveform_output_file)
                    temp_hist = ROOT.TH1I(pmt_waveform.get_pmt_trace_id(),
                                          pmt_waveform.get_pmt_trace_id(),
                                          pmt_waveform.get_pmt_waveform_length(),
                                          0,
                                          pmt_waveform.get_pmt_waveform_length())

                    for i_value in range(pmt_waveform.get_pmt_waveform_length()):
                        temp_hist.SetBinContent(i_value + 1, pmt_waveform.get_pmt_waveform()[i_value])

                    waveform_output_file.cd()
                    temp_hist.Write()
                    del temp_hist

                del pmt_waveform

            new_waveform_bool = False

        line_number_int += 1
    pmt_data_file.close()


def read_filenames(input_file_name: str):
    filenames_list = []
    try:
        file = open(input_file_name, 'r')
    except FileNotFoundError as fnf_error:
        print(fnf_error)
        raise Exception("Error opening data file. Skip to the next file...")
    for index, line in enumerate(file.readlines()):
        if line is not "" or line is not None:
            filenames_list.append(line.split()[0])
    return filenames_list

