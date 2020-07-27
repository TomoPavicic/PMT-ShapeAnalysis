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
from src.PMT_Array import PMT_Array


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

    fit_parameter = [[] for i in range(pmt_array.get_pmt_total_number())]

    for i_om in range(pmt_array.get_pmt_total_number()):
        charge_hist = file.Get(date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() +
                               "_charge_spectrum_" + str(voltage) + "V")
        amp_hist = file.Get(date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() +
                            "_amplitude_spectrum_" + str(voltage) + "V")
        baseline_hist = file.Get(date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() +
                                 "_baseline_distribution_" + str(voltage) + "V")

        try:
            charge_hist.GetEntries()
            amp_hist.GetEntries()
            baseline_hist.GetEntries()
        except:
            continue

        mu_guess = charge_hist.GetMaximumBin() * charge_hist.GetBinWidth(0)
        lower_range = mu_guess - 2
        higher_range = mu_guess + 6

        bi_fit = ROOT.TF1("fit", "[0]*(7.08*TMath::Gaus(x,[1],[2])"
                                 " + 1.84*TMath::Gaus(x,[1]*(1 + 72.144/975.651),[2]*1.036)"
                                 " + 0.44*TMath::Gaus(x,[1]*(1 + 84.154/975.651),[2]*1.042))"
                                 " + [3]*(exp([4]*x)/(1 + exp((x - [5])/[6])))",
                          lower_range, higher_range)
        bi_fit.SetParNames("A", "mu", "sig", "B", "e", "c_e", "s")

        bi_fit.SetParLimits(0, 0, 500)
        bi_fit.SetParLimits(1, mu_guess - 1, mu_guess + 1)
        bi_fit.SetParLimits(2, 0.8, 1.2)
        bi_fit.SetParLimits(3, 0, 500)
        bi_fit.SetParLimits(4, 0.1, 10)
        bi_fit.SetParLimits(5, mu_guess - 6, mu_guess)
        bi_fit.SetParLimits(6, 1, 3)
        bi_fit.SetParameters(100, mu_guess, 1, 400, 0.5, mu_guess - 2, 2)

        name = output_file_location + "/plots/" + date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() + "_charge_spectrum_fit.pdf"

        charge_hist.Fit(bi_fit, "0Q","", lower_range, higher_range)
        charge_hist.SetXTitle("Charge /pC")
        charge_hist.SetYTitle("Counts")
        charge_hist.SetTitle(date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() + "_charge_spectrum_fit")

        c1 = ROOT.TCanvas()
        charge_hist.Draw()
        bi_fit.Draw("same")
        c1.SetGrid()
        c1.Update()
        ROOT.gStyle.SetOptFit(1)
        c1.SaveAs(name)

        '''print(">>>")
        print("Fit output:")
        for i in range(3):
            print(bi_fit.GetParName(i), bi_fit.GetParameter(i), "+/-", bi_fit.GetParError(i))
        print("Chi2/NDoF:", bi_fit.GetChisquare(), "/", bi_fit.GetNDF(), "=", bi_fit.GetChisquare() / bi_fit.GetNDF())
        print(">>>")'''
        if bi_fit.GetNDF() == 0:
            pass
        else:
            pars = {
                "mu": bi_fit.GetParameter(1),
                "mu_err": bi_fit.GetParError(1),
                "sig": bi_fit.GetParameter(2),
                "sig_err": bi_fit.GetParError(2),
                "base_mu": baseline_hist.GetMean(),
                "base_sig": baseline_hist.GetStdDev(),
                "gain": amp_hist.GetMaximumBin(),
                "chi2": bi_fit.GetChisquare() / bi_fit.GetNDF()
            }
            fit_parameter[i_om].append(pars)

    file.Close()
    return fit_parameter


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
    resolutions = [[] for i in range(pmt_array.get_pmt_total_number())]
    resolutions_err = [[] for i in range(pmt_array.get_pmt_total_number())]
    dates = [[] for i in range(pmt_array.get_pmt_total_number())]
    gains = [[] for i in range(pmt_array.get_pmt_total_number())]
    gains_err = [[] for i in range(pmt_array.get_pmt_total_number())]
    baseline_means = [[] for i in range(pmt_array.get_pmt_total_number())]
    baseline_sigs = [[] for i in range(pmt_array.get_pmt_total_number())]
    fit_chi2 = [[] for i in range(pmt_array.get_pmt_total_number())]

    for i_file in tqdm.tqdm(range(filenames.size)):
        file = filenames[i_file][0].decode("utf-8")

        date = file.split("_")[0]
        voltage = int(file.split("_")[1].split("A")[1])

        if voltage == 1000:
            pass
        else:
            continue

        fit_parameters = read_file(date, voltage, input_directory + "/" + file, pmt_array, output_directory)
        write_to_file(fit_parameters)

        for i_om in range(pmt_array.get_pmt_total_number()):

            if len(fit_parameters[i_om]) > 0:
                pass
            else:
                continue

            mu = fit_parameters[i_om][0]["mu"]
            mu_err = fit_parameters[i_om][0]["mu_err"]
            sig = fit_parameters[i_om][0]["sig"]
            sig_err = fit_parameters[i_om][0]["sig_err"]
            chi_2 = fit_parameters[i_om][0]["chi2"]
            gain = fit_parameters[i_om][0]["gain"]
            baseline_mean = fit_parameters[i_om][0]["base_mu"]
            baseline_sig = fit_parameters[i_om][0]["base_sig"]

            res, res_err = get_resolution(mu, mu_err, sig, sig_err)

            resolutions[i_om].append(res)
            resolutions_err[i_om].append(res_err)
            dates[i_om].append(int(date))
            gains[i_om].append(gain)
            baseline_means[i_om].append(baseline_mean)
            baseline_sigs[i_om].append(baseline_sig)
            fit_chi2[i_om].append(chi_2)

    # Plot individual summaries
    for i_om in range(pmt_array.get_pmt_total_number()):

        if len(resolutions[i_om]) > 0:
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

        plt.plot(date[:start + 1], np.array(gains[i_om][:start + 1])*2,
                 "g.", label="Atmospheric He")
        plt.plot(date[start+1:mid + 1], np.array(gains[i_om][start+1:mid + 1])*2,
                 "b.", label="1% He")
        plt.plot(date[mid+1:], np.array(gains[i_om][mid+1:])*2,
                 "r.", label="10% He")
        plt.axvline(date[start], 0, 100, ls='--', color='k')
        plt.axvline(date[mid], 0, 100, ls='--', color='k')
        plt.xlabel("exposure days relative to 190611")
        plt.ylabel("PMT gain at 1Mev /mV")
        plt.title(pmt_array.get_pmt_object_number(i_om).get_pmt_id() + " Gain at 1MeV vs exposure time")
        plt.grid()
        plt.ylim(150,300)
        plt.legend(loc='lower right')
        plt.savefig(output_directory + "/summary_plots/" +
                    pmt_array.get_pmt_object_number(i_om).get_pmt_id() + "_gains_vs_time")
        plt.close()

        plt.errorbar(date, baseline_means[i_om], yerr=baseline_sigs[i_om], fmt='k.-', ecolor='r')
        plt.grid()
        plt.xlabel("exposure days relative to 190611")
        plt.ylabel("Baseline mean /mV")
        plt.title(pmt_array.get_pmt_object_number(i_om).get_pmt_id() + " Baseline mean vs exposure time")
        plt.savefig(output_directory + "/summary_plots/" +
                    pmt_array.get_pmt_object_number(i_om).get_pmt_id() + "_baseline_mean_vs_time")

        plt.close()

        '''plt.plot(date, baseline_sigs[i_om])
        plt.xlabel("exposure days relative to 190611")
        plt.ylabel("Baseline std-dev")
        plt.title(pmt_array.get_pmt_object_number(i_om).get_pmt_id() + " Baseline std-dev vs exposure time")
        plt.savefig(output_directory + "/summary_plots/" +
                    pmt_array.get_pmt_object_number(i_om).get_pmt_id() + "_baseline_sigma_vs_time")
        plt.close()'''

        res_filter = []
        for i_date in range(len(date)):
            if 1 < resolutions[i_om][i_date] < 3.75 and fit_chi2[i_om][i_date] < 10:
                res_filter.append(True)
            else:
                res_filter.append(False)

        popt, pcov = curve_fit(linear, np.array(date[start:])[res_filter[start:]], np.array(resolutions[i_om][start:])[res_filter[start:]],
                               sigma=np.array(resolutions_err[i_om][start:])[res_filter[start:]], p0=[0.001, 2], bounds=[[0, 0], [0.002, 3.5]], maxfev=500000)
        x_array = np.linspace(date[start], np.amax(date), 2)
        chi_2 = chi2(np.array(resolutions[i_om][start:])[res_filter[start:]], np.array(resolutions_err[i_om][start:])[res_filter[start:]], linear(date[start:][res_filter[start:]], *popt), len(popt))

        plt.errorbar(date[:start + 1], resolutions[i_om][:start + 1], yerr=resolutions_err[i_om][:start + 1],
                     fmt="g.", label="Atmospheric He")
        plt.plot(np.array(date[start:])[res_filter[start:]], np.array(resolutions[i_om][start:])[res_filter[start:]], 'ko',
                 label="used values")
        plt.errorbar(date[start+1:mid + 1], resolutions[i_om][start+1:mid + 1], yerr=resolutions_err[i_om][start+1:mid + 1],
                     fmt="b.", label="1% He")
        plt.errorbar(date[mid+1:], resolutions[i_om][mid+1:], yerr=resolutions_err[i_om][mid+1:],
                     fmt="r.", label="10% He")
        plt.plot(x_array, linear(x_array, *popt), 'k-')
        plt.xlabel("exposure days relative to 190611 \n $y = (${:.1e}$ ± ${:.0e})$x + (${:.1e}$ ± ${:.0e}$)$ $\chi^2_R = {:.2}$".format(popt[0], np.sqrt(pcov[0, 0]),popt[1], np.sqrt(pcov[1, 1]), chi_2))
        plt.ylabel("Resolution at 1MeV /% $\sigma / \mu$")
        plt.title(pmt_array.get_pmt_object_number(i_om).get_pmt_id() + " Resolution vs exposure time")
        plt.axvline(date[start], 0, 100, ls='--', color='k')
        plt.axvline(date[mid], 0, 100, ls='--', color='k')
        plt.grid()
        plt.ylim(2, 4.5)
        plt.xlim(np.amin(date), np.amax(date))
        plt.legend(loc='upper right')
        plt.tight_layout()
        plt.savefig(output_directory + "/summary_plots/" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() +
                    "_resolution_vs_time")
        plt.close()

        plt.plot(np.array(date[start:])[res_filter[start:]], np.array(fit_chi2[i_om][start:])[res_filter[start:]],
                 'ko', label="used values")
        plt.plot(date[:start + 1], fit_chi2[i_om][:start + 1], "g.", label="Atmospheric He")
        plt.plot(date[start + 1:mid + 1], fit_chi2[i_om][start + 1:mid + 1], "b.", label="1% He")
        plt.plot(date[mid + 1:], fit_chi2[i_om][mid + 1:], "r.", label="10% He")
        plt.xlabel("exposure days relative to 190611")
        plt.ylabel("$\chi^2_R$")
        plt.title(pmt_array.get_pmt_object_number(i_om).get_pmt_id() + " Resolution fit $\chi^2_R$ vs exposure time")
        plt.grid()
        plt.axvline(date[start], 0, 100, ls='--', color='k')
        plt.axvline(date[mid], 0, 100, ls='--', color='k')
        plt.legend(loc='upper right')
        plt.xlim(np.amin(date), np.amax(date))
        plt.ylim(0, 10)
        plt.tight_layout()
        plt.savefig(output_directory + "/summary_plots/" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() +
                    "_resolution_vs_time_chi2")
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
                ratio.append(resolutions[0][i] / resolutions[1][j])
                ratio_err.append(resolutions[0][i] / resolutions[1][j] * np.sqrt(
                    (resolutions_err[0][i] / resolutions[0][i]) ** 2 + (resolutions_err[1][j] / resolutions[1][j]) ** 2))
                gain_ratio.append(gains[0][i]/gains[1][j])
                break

    popt, pcov = curve_fit(linear, x_date, ratio,
                           sigma=ratio_err, p0=[1, 1], bounds=[[0, 0], [10, 10]])
    x_date = process_date(x_date)
    x_array = np.linspace(np.amin(x_date), np.amax(x_date), 2)
    chi_2 = chi2(ratio, ratio_err, linear(x_date, *popt), 2)

    plt.errorbar(x_date, ratio, yerr=ratio_err, fmt="k.")
    plt.plot(x_array, linear(x_array, *popt), "g-",
             label="$y = (${:.1e}$ ± ${:.0e})$\\times x + (${:.1e}$ ± ${:.0e}$) \chi^2_R = {:.2}$".format(popt[0], np.sqrt(pcov[0, 0]),popt[1], np.sqrt(pcov[1, 1]), chi_2))
    plt.axvline(98, color="r", ls="--")
    plt.axvline(0, color="b", ls="--")
    plt.xlabel("exposure days relative to 190611")
    plt.ylabel("Ratio res_Ch0/res_Ch1")
    plt.title("Ratio of resolution of CH 0 & 1 vs time")
    plt.grid()
    plt.xlim(np.amin(np.array(x_date)), np.amax(np.array(x_date)))
    plt.ylim(0, 2)
    plt.savefig(output_directory + "/summary_plots/resolution_ratio_vs_time")
    plt.close()

    plt.plot(x_date, gain_ratio, "k.")
    plt.axvline(98, color="r", ls="--")
    plt.axvline(0, color="b", ls="--")
    plt.xlabel("exposure days relative to 190611")
    plt.ylabel("Ratio gain_Ch0/gain_Ch1")
    plt.title("Ratio of gain of CH 0 & 1 vs time")
    plt.grid()
    plt.xlim(np.amin(np.array(x_date)), np.amax(np.array(x_date)))
    plt.ylim(0.7, 1)
    plt.savefig(output_directory + "/summary_plots/gain_ratio_vs_time")
    plt.close()


if __name__ == '__main__':
    main()
