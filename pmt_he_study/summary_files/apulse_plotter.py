import sys

sys.path.insert(1, '..')

# import ROOT and bash commands
import ROOT
import tqdm

# import python plotting and numpy modules
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks

# import stats module
from scipy.optimize import curve_fit

# import custom made classes
from functions.other_functions import pmt_parse_arguments, fit, chi2, process_date, linear, gaus
from scr.PMT_Array import PMT_Array


def write_to_file(fit_parameters):
    pass


def get_error_a_divide_b(da, a, db, b):
    c = a/b
    dc = c * np.sqrt((da/a)**2 + (db/b)**2)
    return dc


def get_resolution(mu: float, mu_err: float, sig: float, sig_err: float):
    res = 0
    res_err = 0
    if mu == 0:
        pass
    else:
        res = sig/mu
        res_err = get_error_a_divide_b(sig_err, sig, mu_err, mu)
    return res*100, res_err*100


def read_file(date: str, voltage: int, root_file_name: str, pmt_array: PMT_Array, output_file_location: str):

    file = ROOT.TFile(root_file_name, "READ")
    file.cd()

    apulse_info = [[] for i in range(pmt_array.get_pmt_total_number())]

    for i_om in range(pmt_array.get_pmt_total_number()):
        apulse_num_hist = file.Get(date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() +
                                   "_apulse_num_" + str(voltage) + "V")
        apulse_time_hist = file.Get(date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() +
                                    "_apulse_times_" + str(voltage) + "V")
        apulse_amplitude_hist = file.Get(date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() +
                                          "_apulse_amplitudes_" + str(voltage) + "V")

        try:
            apulse_num_hist.GetEntries()
            apulse_time_hist.GetEntries()
            apulse_amplitude_hist.GetEntries()
        except:
            continue

        apulse_rate = 0
        for i_bin in range(2, apulse_num_hist.GetNbinsX()):
            apulse_rate += apulse_num_hist.GetBinContent(i_bin)

        apulse_rate = apulse_rate/apulse_num_hist.GetEntries() * 100
        '''c1 = ROOT.TCanvas()
        charge_hist.Draw()
        bi_fit.Draw("same")
        c1.SetGrid()
        c1.Update()
        ROOT.gStyle.SetOptFit(1)
        c1.SaveAs(name)'''

        pars = {
            "apulse_rate": apulse_rate
        }
        apulse_info[i_om].append(pars)

    file.Close()
    return apulse_info


def main():
    # Handle the input arguments:
    ##############################
    args = pmt_parse_arguments()
    input_directory = args.i
    # config_file_name = args.c
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

    '''# Set the cuts you wish to apply
    # If you don't do this the defaults are used
    if config_file_name is not None:
        pmt_array.apply_setting(config_file_name)
        # print_settings(pmt_array)'''

    # Set up the containers for the summary
    apulse_rates = [[] for i in range(pmt_array.get_pmt_total_number())]
    apulse_rates_err = [[] for i in range(pmt_array.get_pmt_total_number())]
    dates = [[] for i in range(pmt_array.get_pmt_total_number())]

    for i_file in tqdm.tqdm(range(filenames.size)):
        file = filenames[i_file][0].decode("utf-8")

        date = file.split("_")[0]
        voltage = int(file.split("_")[1].split("A")[1])

        if voltage == 1400:
            pass
        else:
            continue

        apulse_info = read_file(date, voltage, input_directory + "/" + file, pmt_array, output_directory)
        write_to_file(apulse_info)

        for i_om in range(pmt_array.get_pmt_total_number()):

            if len(apulse_info[i_om]) > 0:
                pass
            else:
                continue

            apulse_rate = apulse_info[i_om][0]["apulse_rate"]

            apulse_rates[i_om].append(apulse_rate)
            apulse_rates_err[i_om].append(0)
            dates[i_om].append(int(date))

    # Plot individual summaries
    for i_om in range(pmt_array.get_pmt_total_number()):

        if len(apulse_rates[i_om]) > 0:
            pass
        else:
            continue
        date = process_date(dates[i_om])

        try:
            start = np.where(date == 0)[0][0]
        except:
            start = np.where(date == 1)[0][0]
        mid = np.where(date == 98)[0][0]

       # print("start:",start)

        plt.plot(date[:start + 1], np.array(apulse_rates[i_om][:start + 1]),
                 "g.", label="Atmospheric He")
        plt.plot(date[start+1:mid + 1], np.array(apulse_rates[i_om][start+1:mid + 1]),
                 "b.", label="1% He")
        plt.plot(date[mid+1:], np.array(apulse_rates[i_om][mid+1:]),
                 "r.", label="10% He")
        plt.axvline(date[start], 0, 100, ls='--', color='k')
        plt.axvline(date[mid], 0, 100, ls='--', color='k')
        plt.xlabel("exposure days relative to 190611")
        plt.ylabel("Afterpulse rate /%")
        plt.title(pmt_array.get_pmt_object_number(i_om).get_pmt_id() + " afterpulse rate vs exposure time")
        plt.grid()
        #plt.ylim(150,300)
        plt.legend(loc='lower right')
        plt.savefig(output_directory + "/summary_plots/" +
                    pmt_array.get_pmt_object_number(i_om).get_pmt_id() + "_apulse_rate_vs_time")
        plt.close()

    # Plot ratio
    x_date = []
    ratio = []
    ratio_err = []
    gain_ratio = []
    for i in range(len(dates[0])):
        for j in range(len(dates[1])):
            if dates[0][i] == dates[1][j]:
                x_date.append(dates[0][i])
                ratio.append(apulse_rates[0][i] / apulse_rates[1][j])
                break

    x_date = process_date(x_date)

    plt.plot(x_date, ratio, "k.")
    plt.axvline(98, color="r", ls="--")
    plt.axvline(0, color="b", ls="--")
    plt.xlabel("exposure days relative to 190611")
    plt.ylabel("Ratio apulse rate Ch0/Ch1")
    plt.title("Ratio of after pulse rates of CH 0 & 1 vs time")
    plt.grid()
    plt.xlim(np.amin(np.array(x_date)), np.amax(np.array(x_date)))
    #plt.ylim(0, 2)
    plt.savefig(output_directory + "/summary_plots/apulse_rate_ratio_vs_time")
    plt.close()


if __name__ == '__main__':
    main()
