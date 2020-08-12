import sys

sys.path.insert(1, '..')

# import ROOT and bash commands
import ROOT
import tqdm

import numpy as np

# import custom made classes
from functions.other_functions import pmt_parse_arguments
from src.PMT_Array import PMT_Array


def print_settings(pmt_array: PMT_Array):
    pass


def read_tree(root_file_name: str, pmt_array: PMT_Array, output_file_location: str, output_file_name: str):

    file = ROOT.TFile(root_file_name, "READ")
    file.cd()

    date = root_file_name.split("/")[-1].split("_")[0]
    voltage = int(root_file_name.split("/")[-1].split("_")[1].split("A")[1])

    if voltage == 1000:
        max_amp = 400
    else:
        max_amp = 1000
    max_charge = 60

    tree = file.T

    # Create the histograms that we will want to store
    # these will be used to extract the resolution
    charge_hists = []
    amp_hists = []
    baselines = []
    apulse_nums_hists = []
    apulse_times_hists = []
    apulse_amplitudes_hists = []
    for i_om in range(pmt_array.get_pmt_total_number()):
        nbins = pmt_array.get_pmt_object_number(i_om).get_setting("nbins")

        charge_hist = ROOT.TH1D(date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() + "_charge_spectrum_" + str(voltage) + "V",
                                date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() + "_charge_spectrum_" + str(voltage) + "V",
                                nbins, 0, max_charge)
        amp_hist = ROOT.TH1D(date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() + "_amplitude_spectrum_" + str(voltage) + "V",
                             date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() + "_amplitude_spectrum_" + str(voltage) + "V",
                             nbins, 0, max_amp)
        baseline = ROOT.TH1D(date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() + "_baseline_distribution_" + str(voltage) + "V",
                             date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() + "_baseline_distribution_" + str(voltage) + "V",
                             nbins, 978, 981)
        apulse_num = ROOT.TH1I(date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() + "_apulse_num_" + str(voltage) + "V",
                               date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() + "_apulse_num_" + str(voltage) + "V",
                               20, 0, 20)
        apulse_time = ROOT.TH1I(date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() + "_apulse_times_" + str(voltage) + "V",
                                date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() + "_apulse_times_" + str(voltage) + "V",
                                7000, 0, 7000)
        apulse_amplitude = ROOT.TH1D(date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() + "_apulse_amplitudes_" + str(voltage) + "V",
                                     date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() + "_apulse_amplitudes_" + str(voltage) + "V",
                                     nbins, 0, 500)
        charge_hists.append(charge_hist)
        amp_hists.append(amp_hist)
        baselines.append(baseline)
        apulse_nums_hists.append(apulse_num)
        apulse_times_hists.append(apulse_time)
        apulse_amplitudes_hists.append(apulse_amplitude)

    for event in tree:

        # Access the information inside the NTuple
        OM_ID = event.OM_ID
        if OM_ID == 0 or OM_ID == 1:
            pass
        else:
            continue

        pulse_amplitude = int(event.pulse_amplitude)
        pulse_charge = event.pulse_charge
        pulse_baseline = event.pulse_baseline
        apulse_num = event.apulse_num
        apulse_times = np.array(event.apulse_times)
        apulse_amplitudes = np.array(event.apulse_amplitudes)
        apulse_shapes = np.array(event.apulse_shapes)

        charge_hists[OM_ID].Fill(pulse_charge)
        amp_hists[OM_ID].Fill(pulse_amplitude)
        baselines[OM_ID].Fill(pulse_baseline)

        # Now apply new amplitude and shape cuts
        new_apulse_num = 0
        filer_list = []
        try:
            for i_apulse in range(apulse_num):
                if apulse_shapes[i_apulse] > pmt_array.get_pmt_object_number(OM_ID).get_setting("mf_shape_threshold")\
                        and apulse_amplitudes[i_apulse] > pmt_array.get_pmt_object_number(OM_ID).get_setting("mf_amp_threshold")\
                        and apulse_times[i_apulse] > pmt_array.get_pmt_object_number(OM_ID).get_setting("sweep_range")[0]:
                    filer_list.append(True)
                    new_apulse_num += 1
                    apulse_amplitudes_hists[OM_ID].Fill(apulse_amplitudes[i_apulse])
                    apulse_times_hists[OM_ID].Fill(apulse_times[i_apulse])
                else:
                    filer_list.append(False)
            apulse_nums_hists[OM_ID].Fill(new_apulse_num)
        except:
            pass

    output_file = ROOT.TFile(output_file_location + "/ROOT_files/" + str(voltage) + "V/" + output_file_name, "RECREATE")
    output_file.cd()

    for i_om in range(pmt_array.get_pmt_total_number()):
        if charge_hists[i_om].GetEntries() > 0:
            charge_hists[i_om].Write()
            amp_hists[i_om].Write()
            baselines[i_om].Write()
            apulse_nums_hists[i_om].Write()
            apulse_times_hists[i_om].Write()
            apulse_amplitudes_hists[i_om].Write()
        else:
            pass

    file.Close()
    output_file.Close()


def main():
    # Handle the input arguments:
    ##############################
    args = pmt_parse_arguments()
    input_directory = args.i
    config_file_name = args.c
    output_directory = args.o
    ##############################

    filenames_txt = input_directory + "/filenames.txt"

    try:
        print(">>> Reading data from file: {}".format(filenames_txt))
        date_file = open(filenames_txt, 'r')
    except FileNotFoundError as fnf_error:
        print(fnf_error)
        raise Exception("Error opening data file {}".format(filenames_txt))

    filenames = np.loadtxt(filenames_txt, delimiter=',', dtype={
        'names': ['filename'],
        'formats': ['S100']}, unpack=True)

    topology = [2, 1]
    pmt_array = PMT_Array(topology, "summary")
    pmt_array.set_pmt_id("GAO607", 0)
    pmt_array.set_pmt_id("GAO612", 1)

    # Set the cuts you wish to apply
    # If you don't do this the defaults are used
    if config_file_name is not None:
        pmt_array.apply_setting(config_file_name)
        print_settings(pmt_array)

    for i_file in tqdm.tqdm(range(filenames.size)):
        file = filenames[i_file][0].decode("utf-8")

        read_tree(input_directory + "/" + file, pmt_array, output_directory, file.split(".root")[0] + "_output.root")


if __name__ == '__main__':
    main()
