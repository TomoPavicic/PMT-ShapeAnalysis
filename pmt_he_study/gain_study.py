import sys
sys.path.insert(1, '..')

import numpy as np
from functions.data_reader_functions import read_filenames
from functions.other_functions import io_parse_arguments, process_date, process_exposure, get_date_time, fit, chi2
import matplotlib.pyplot as plt
# import stats module
from scipy.optimize import curve_fit
import ROOT


def main():
    args = io_parse_arguments()
    input_data_filename = args.i
    data_file_list = read_filenames(input_data_filename)
    PATH = input_data_filename.split("filenames")[0]

    print(input_data_filename)

    date_array = [[], []]
    mu_array = [[], []]
    err_array = [[], []]

    for data_file_index, data_file_name in enumerate(data_file_list):
        try:
            _file = open(PATH + data_file_name, 'r')
        except:
            continue

        file = ROOT.TFile(PATH + data_file_name, "READ")

        area_ch0 = file.Get("Area_Spectra_Ch0;1")
        amp_ch0 = file.Get("Amplitude_Spectra_Ch0;1")
        area_ch1 = file.Get("Area_Spectra_Ch1;1")
        amp_ch1 = file.Get("Amplitude_Spectra_Ch1;1")

        area_hists = [area_ch0, area_ch1]
        amp_hists = [amp_ch0, amp_ch1]

        date, time = get_date_time(data_file_name)

        try:
            area_ch0.GetEntries()
            area_ch1.GetEntries()
            print(data_file_index)
            print("CH0 Events: ", area_ch0.GetEntries())
            print("CH1 Events: ", area_ch1.GetEntries())

        except:
            continue

        for i_pmt in range(2):

            if amp_hists[i_pmt].GetEntries() > 0:

                date_array[i_pmt].append(date)

                x_0 = []
                y_0 = []
                y_err_0 = []
                for ibin in range(1, amp_ch0.GetEntries()):
                    x_0.append(ibin * amp_ch0.GetBinWidth(5))
                    y_0.append(amp_ch0.GetBinContent(ibin))
                    if amp_ch0.GetBinContent(ibin) > 0:
                        y_err_0.append(np.sqrt(amp_ch0.GetBinContent(ibin)))
                    else:
                        y_err_0.append(0)
                x_0 = np.array(x_0)
                y_0 = np.array(y_0)
                y_err_0 = np.array(y_err_0)

                range_lower = np.argmax(y_0 - int(5 * amp_ch0.GetBinWidth(5)))
                range_higher = np.argmax(y_0) + int(10 * amp_ch0.GetBinWidth(5))

                y_array = y_0[range_lower:range_higher]
                x_array = x_0[range_lower:range_higher]

                y_err = y_err_0[range_lower:range_higher]

                est_mu = x_array[np.argmax(y_array)]
                p_guess = [est_mu, 1, 50, 0.5, 0.2, est_mu, 2]

                print(p_guess)
                p_bounds = [[est_mu - 10, 0, 0, 0, 0, est_mu - 10, 0], [est_mu + 10, 10, 1000, 1000, 1000, est_mu + 10, 1000]]

                popt, pcov = curve_fit(fit, x_array, y_array, p0=p_guess, bounds=p_bounds, maxfev=5000)

                string = '{},'.format(i_pmt)

                print("\n>>> The optimised fitted parameters are: ", popt)
                print(">>> The covariance matrix is: ", pcov)
                for i in range(len(popt)):
                    string += '{},{},'.format(popt[i], np.sqrt(pcov[i][i]))
                    print(">>> Error on parameter {}: {} is {}".format(i, popt[i], np.sqrt(pcov[i][i])))

                chi_2 = chi2(y_array, y_err, fit(x_array, *popt), 3)
                print(">>> The reduced chi2 is: ", chi_2)
                string += '{}\n'.format(chi_2)

                x = np.linspace(np.min(x_array), np.max(x_array), 10000)

                # plt.hist(pulse_charges[i_pmt], bins=150, range=(0, 50), facecolor='b', alpha=0.25)
                plt.errorbar(x_array, y_array, yerr=y_err, fmt='.')
                plt.plot(x, fit(x, *popt))
                # plt.text(38, 1000, "$\chi^2$ = {}".format(round(chi_2, 2)))
                plt.grid(True)
                plt.title('Ch' + str(i_pmt) + '_Amplitude_Spectrum')
                plt.ylabel("Counts")
                plt.xlabel("Amplitude /mV")
                plt.show()
                # plt.xlim(est_mu - 25, est_mu + 14)
                #plt.yscale('log')
                #plt.savefig('/home/wquinn/pmt_analysis/plots/' + pmt_array.get_pmt_object_number(i_pmt).get_pmt_id() + '_Charge_Spectrum.png')

        del area_ch0
        del area_ch1
        del amp_ch0
        del amp_ch1
        del file


if __name__ == '__main__':
    main()
