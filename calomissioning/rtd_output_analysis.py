import sys
sys.path.insert(1, '..')

import ROOT
import tqdm
import numpy as np
import matplotlib.pyplot as plt
from functions.other_functions import *
from src.PMT_Array import PMT_Array
from src.PMT_Waveform import PMT_Waveform


def get_rise_time(waveform, ratio):
    # waveform need baseline of zero
    peak = int(np.argmin(waveform))
    point = 0
    for i in range(peak, waveform.size):
        if waveform[i] > np.amin(waveform)*ratio:
            point = i
            break

    return point - peak


def write_average_waveforms(average_waveforms: list, write_file: str):
    print("")
    print(">>> Writing templates to file: ", write_file)
    print("")

    output_file = ROOT.TFile(write_file, "RECREATE")
    output_file.cd()

    for i in tqdm.tqdm(range(len(average_waveforms))):
        temp_hist = ROOT.TH1D("Template_"+str(i), "Template_"+str(i),
                              len(average_waveforms[i]), 0, len(average_waveforms[i]))

        norm = np.sqrt(np.dot(average_waveforms[i], average_waveforms[i]))

        for j in range(len(average_waveforms[i])):
            temp_hist.SetBinContent(j, average_waveforms[i][j]/norm)

        temp_hist.Write()
        del temp_hist

    output_file.Close()

    return


def read_average_waveforms(temp_file_name: str, num_pmts: int):
    print("")
    print(">>> Reading templates from file: ", temp_file_name)
    print("")

    temp_file = ROOT.TFile(temp_file_name, "READ")
    temp_file.cd()

    templates = [[] for k in range(num_pmts)]

    for i in range(num_pmts):
        temp_hist = temp_file.Get("Template_"+str(i))

        for j in range(1, int(temp_hist.GetNbinsX()) + 1):
            templates[i].append(temp_hist.GetBinContent(j))

        assert len(templates[i]) == 400
    return templates


def read_cable_lengths():
    cable_lengths = []
    return cable_lengths


def main():

    # Handle the file inputs
    args = sncd_parse_arguments()
    input_data_filename = args.i
    output_data_filename = args.o
    config_file_name = args.c
    main_wall = args.w

    if main_wall == 'fr':
        main_wall = 1
    elif main_wall == 'it':
        main_wall = 0

    num_events = 1000000
    num_bins = 100
    max_amp = 1000

    # isolate the run number for naming convienience
    run_number = input_data_filename.split("_")[1]

    # Check to see if file exists - standard
    try:
        data_file = open(input_data_filename, 'r')
    except FileNotFoundError as fnf_error:
        print(fnf_error)
        raise Exception("Error opening data file")

    # Set up the pmt array
    topology = [13, 20]
    num_pmts = topology[0]*topology[1]
    pmt_array = PMT_Array(topology, run_number)
    for i in range(topology[0]):
        for j in range(topology[1]):
            num = i + topology[0]*j
            pmt_array.set_pmt_id("M:{}.{}.{}".format(main_wall, j, i), num)

    # Configure the array of PMTs - not as important here just yet
    if config_file_name is not None:
        pmt_array.apply_setting(config_file_name)

    cable_lengths = read_cable_lengths()

    # Open file
    file = ROOT.TFile(input_data_filename, "READ")
    file.cd()

    # The tree inside is called "T"
    tree = file.T

    # Counter for testing
    event_counter = [0 for i in range(num_pmts)]

    # raw_amplitudes = [[] for i in range(num_pmts)]

    '''templates = read_average_waveforms("~/Desktop/test_template.root", num_pmts)
    amp_hists = ROOT.TList()

    output_file = ROOT.TFile(output_data_filename, "RECREATE")

    for i_om in tqdm.tqdm(range(num_pmts)):
        temp_hist = ROOT.TH1D(pmt_array.get_pmt_object_number(i_om).get_pmt_id() + "_amplitude_index_spectrum",
                              pmt_array.get_pmt_object_number(i_om).get_pmt_id() + "_amplitude_index_spectrum",
                              num_bins, 0, max_amp)
        amp_hists.Add(temp_hist)'''

    h_map = ROOT.TH2I("event_map", "event_map", topology[1], 0, topology[1], topology[0], 0, topology[0])
    
    # create_log("run_167_output.txt")

    x = []

    # Run over file to create containers
    for event in tqdm.tqdm(tree):
        
        event_num = event.event_num
        row = event.row
        col = event.column
        OM_ID = event.OM_ID
        charge = -1*event.charge
        baseline = event.baseline
        amplitude = -1*event.amplitude
        fall_time = event.fall_time
        rise_time = event.rise_time
        peak_time = event.peak_time
        calo_hit_num = event.calo_hit_num
        calo_tdc = event.calo_tdc
        run_num = event.run_num
        wall_num = event.wall_num
        trig_num = event.trig_id

        if col == 9 and row == 7:
            x.append(trig_num)
        
        '''output_log("run_167_output.txt", "{},{},{},{},{},"
                                         "{},{},{},{},{},"
                                         "{},{},{},{},{}".format(event_num, row, col, OM_ID, charge,
                                                                 baseline, amplitude, fall_time, rise_time,
                                                                 peak_time, calo_hit_num, calo_tdc, run_num,
                                                                 wall_num, trig_num))'''

        #print("calo_hit_num:", calo_hit_num, "event_num:", event_num, "trig_num:", trig_num)
        # raw_amplitudes[OM_ID].append(amplitude)
        
        h_map.Fill(col, row, 1)

        '''pmt_waveform = PMT_Waveform(event.waveform, pmt_array.get_pmt_object_number(OM_ID))
        peak = pmt_waveform.get_pmt_pulse_peak_position()

        if peak > 51:
            try:
                amplitude_index = np.dot(pmt_waveform.get_pmt_waveform_reduced()[peak - 50:peak + 350], templates[OM_ID])
                # raw_amplitudes[OM_ID].append(amplitude_index)
                amp_hists[OM_ID].Fill(amplitude_index)

            except:
                print("Waveform", len(pmt_waveform.get_pmt_waveform_reduced()[peak - 50:peak + 350]))
                print("Template", len(templates[OM_ID]))'''

        '''if event_counter[0] == num_events:
            break
        event_counter[OM_ID] += 1'''

        # del pmt_waveform

    sel_evnts = []
    for i, val in enumerate(x):
        if val in sel_evnts:
            pass
        else:
            sel_evnts.append(val)

    y = [[] for i in range(len(sel_evnts))]
    sel_evnts = np.array(sel_evnts)

    for event in tqdm.tqdm(tree):
        event_num = event.event_num
        row = event.row
        col = event.column
        OM_ID = event.OM_ID
        charge = -1 * event.charge
        baseline = event.baseline
        amplitude = -1 * event.amplitude
        fall_time = event.fall_time
        rise_time = event.rise_time
        peak_time = event.peak_time
        calo_hit_num = event.calo_hit_num
        calo_tdc = event.calo_tdc
        run_num = event.run_num
        wall_num = event.wall_num
        trig_num = event.trig_id
        pulse_time = event.pulse_time

        val = calo_tdc*6.25 - 400 + pulse_time - cable_lengths[OM_ID]

        if len(np.where(sel_evnts == trig_num)[0]) > 0:
            t = np.where(sel_evnts == trig_num)[0][0]
            y[t].append([trig_num, val])

    '''amp_bins = [i*max_amp/num_bins for i in range(num_bins)]
    amp_cuts = []

    # output_file.cd()
    for i_om in tqdm.tqdm(range(num_pmts)):
        #amp_hists[i_om].SaveAs("/Users/willquinn/Desktop/PDFs/amplitude_index_"+str(i_om)+".pdf")
        #amp_hists[i_om].Write()
        plt.hist(raw_amplitudes[i_om], bins=amp_bins, color='g')
        plt.xlabel("Amplitude /mV")
        plt.ylabel("Counts")
        plt.title(pmt_array.get_pmt_object_number(i_om).get_pmt_id() + "_amplitude_spectrum Events: " + str(event_counter[i_om]))
        plt.savefig("/Users/willquinn/Desktop/PDFs/om_" + str(i_om) + "_amplitude_spectrum_HT")
        plt.grid()
        plt.yscale('log')
        plt.close()
    # output_file.Close()'''

    ROOT.gStyle.SetOptStat(0)
    canvas = ROOT.TCanvas()
    canvas.cd()
    h_map.Draw("colztext")
    canvas.Update()
    canvas.SaveAs("~/Desktop/map.png")
    
    return

    '''raw_amplitudes_array = []
    for i_om in tqdm.tqdm(range(num_pmts)):

        plt.hist(raw_amplitudes[i_om], bins=amp_bins, color='g', log=True)
        plt.xlabel("amplitude /mV")
        plt.ylabel("counts")
        plt.grid()
        plt.xlim(1500, max_amp)
        plt.title("OM " + pmt_array.get_pmt_object_number(i_om).get_pmt_id() + " Events: " + str(event_counter[i_om]))
        plt.savefig("/Users/willquinn/Desktop/PDFs/amplitude_index_spectrum_"+str(i_om))
        plt.close()

        pass'''

    '''temp_array, bin_edges = np.histogram(raw_amplitudes[i_om], bins=amp_bins)
        peak_amp_pos = np.argmax(temp_array)
        # peak_amp = np.amax(temp_array)
        amp_cuts.append(peak_amp_pos * max_amp/num_bins)
        # raw_amplitudes_array.append(temp_array)'''

    '''template_waveforms = [np.zeros(400) for j in range(num_pmts)]
    average_counter = [0 for j in range(num_pmts)]

    for event in tqdm.tqdm(tree):
        OM_ID = event.OM_ID
        amplitude = -1 * event.amplitude / 8

        lower_cut = amp_cuts[OM_ID] - 10*max_amp/num_bins
        higher_cut = amp_cuts[OM_ID] + 10 * max_amp / num_bins

        if average_counter[OM_ID] < 100:
            if lower_cut < amplitude < higher_cut:
                pmt_waveform = PMT_Waveform(event.waveform, pmt_array.get_pmt_object_number(OM_ID))
                peak = pmt_waveform.get_pmt_pulse_peak_position()
                template_waveforms[OM_ID] += np.array(pmt_waveform.get_pmt_waveform_reduced()[peak-50:peak+350])
                average_counter[OM_ID] += 1
        else:
            pass

    x = np.array([i/2.56 for i in range(400)])

    for i_om in tqdm.tqdm(range(num_pmts)):
        if average_counter[i_om] > 0:
            template_waveforms[i_om] = template_waveforms[i_om]/average_counter[i_om]
        else:
            template_waveforms[i_om] = np.zeros(400)

        plt.plot(x, template_waveforms[i_om], 'r-')
        plt.grid()
        plt.xlabel("Timestamp /ns")
        plt.ylabel("Averaged ADC count /mV")
        plt.title("OM " + pmt_array.get_pmt_object_number(i_om).get_pmt_id() + " Template")
        plt.savefig("/Users/willquinn/Desktop/PDFs/template_"+str(i_om))
        plt.close()

    write_average_waveforms(template_waveforms, "~/Desktop/test_template.root")'''


if __name__ == '__main__':
    main()
