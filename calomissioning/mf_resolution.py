import sys
sys.path.insert(1, '..')

import ROOT
import tqdm
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from functions.other_functions import sncd_parse_arguments, gaussian, chi2, inner_product
from scr.PMT_Array import PMT_Array
from scr.PMT_Waveform import PMT_Waveform


def tl_func(x, par):
    f = par[0] / 2 * (1 + ROOT.TMath.Erf((par[1] - x) / (par[2] * np.sqrt(2)) ))
    return f


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

    # isolate the run number for naming convienience
    run_number = input_data_filename.split("_")[1]

    # Check to see if file exists - standard
    try:
        data_file = open(input_data_filename, 'r')
    except FileNotFoundError as fnf_error:
        print(fnf_error)
        raise Exception("Error opening data file")

    # Set up the pmt array
    topology = [20, 13]
    num_pmts = topology[0]*topology[1]
    pmt_array = PMT_Array(topology, run_number)
    for i in range(topology[0]):
        for j in range(topology[1]):
            num = i + topology[0]*j
            pmt_array.set_pmt_id("M:{}.{}.{}".format(main_wall, j, i), num)

    # Configure the array of PMTs - not as important here just yet
    if config_file_name is not None:
        pmt_array.apply_setting(config_file_name)

    # Open file
    file = ROOT.TFile(input_data_filename, "READ")
    file.cd()

    for i_om in tqdm.tqdm(range(num_pmts)):
        temp_hist = file.Get(pmt_array.get_pmt_object_number(i_om).get_pmt_id()
                             + "_amplitude_index_spectrum")

        fit = ROOT.TF1("fit", "[0] / 2 * (1 + TMath::Erf(([1] - x) / ([2] * TMath::Sqrt(2)) ))", 4500, 8000)
        fit.SetParNames("A", "mu", "sig")

        #fit.SetParLimits(0, 0, 400)
        #fit.SetParLimits(1, bi207_1MeV_peak_position - 2, bi207_1MeV_peak_position + 1)
        #fit.SetParLimits(2, 0.8, 2)
        fit.SetParameters(1000, 5000, 10000)

        temp_hist.Fit(fit, "0Q")
        temp_hist.SetXTitle("Amplitude index")
        temp_hist.SetYTitle("Counts")

        name = "/Users/willquinn/Desktop/test_PDFs/om_" + str(i_om) + "_amplitude_spectrum.png"

        c1 = ROOT.TCanvas()
        temp_hist.Draw()
        fit.Draw("same")
        c1.SetLogy()
        c1.SetGrid()
        c1.Update()
        ROOT.gStyle.SetOptFit(1)
        c1.SaveAs(name)

        print(">>>")
        print("Fit output:")
        print("[0]:", fit.GetParameter(0), "+/-", fit.GetParError(0))
        print("[1]:", fit.GetParameter(1), "+/-", fit.GetParError(1))
        print("[2]:", fit.GetParameter(2), "+/-", fit.GetParError(2))
        print("Chi2/NDoF:", fit.GetChisquare(), "/", fit.GetNDF(), "=", fit.GetChisquare() / fit.GetNDF())
        print(">>>")
        break


if __name__ == '__main__':
    main()
