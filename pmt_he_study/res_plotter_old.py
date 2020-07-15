import sys

sys.path.insert(1, '..')

import matplotlib.pyplot as plt
import numpy as np
from functions.other_functions import io_parse_arguments, process_date, process_exposure, linear, chi2
from scipy.optimize import curve_fit


def main():
    args = io_parse_arguments()
    data_file_Ch0_name = args.i + "resolution_vs_date_Ch0.txt"
    data_file_Ch1_name = args.i + "resolution_vs_date_Ch1.txt"

    try:
        data_file_Ch0 = open(data_file_Ch0_name, 'r')
        data_file_Ch1 = open(data_file_Ch1_name, 'r')
    except FileNotFoundError as fnf_error:
        print(fnf_error)
        raise Exception("Error opening data file")

    data_Ch0 = np.loadtxt(data_file_Ch0_name, delimiter=',', dtype={
        'names': ('unix', 'date', 'timestamp', 'mu', 'mu_err', 'sigma', 'sigma_err', 'res', 'chi_err'),
        'formats': ('i4', 'i4', 'i4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4')}, unpack=True)

    data_Ch1 = np.loadtxt(data_file_Ch1_name, delimiter=',',
                          dtype={'names': (
                              'unix', 'date', 'timestamp', 'mu', 'mu_err', 'sigma', 'sigma_err', 'res', 'chi_err'),
                              'formats': ('i4', 'i4', 'i4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4')}, unpack=True)

    data = [data_Ch0, data_Ch1]
    dates = [data_Ch0[1], data_Ch1[1]]
    day_nums = [process_date(data_Ch0[1]), process_date(data_Ch1[1])]
    exposures = [process_exposure(data_Ch0[1]), process_exposure(data_Ch1[1])]
    mus = [data_Ch0[3], data_Ch1[3]]
    mu_errs = [data_Ch0[4], data_Ch1[4]]
    sigs = [data_Ch0[5], data_Ch1[5]]
    sig_errs = [data_Ch0[6], data_Ch1[6]]
    res = [data_Ch0[7], data_Ch1[7]]
    res_errs = [res[0] * np.sqrt((sig_errs[0] / sigs[0]) ** 2 + (mu_errs[0] / mus[0]) ** 2),
               res[1] * np.sqrt((sig_errs[1] / sigs[1]) ** 2 + (mu_errs[1] / mus[1]) ** 2)]
    chis = [data_Ch0[8], data_Ch1[8]]

    p_guess = [0.02420914, 3]
    p_bounds = [[0, 0], [1, 4]]

    filter_lists = [[], []]

    for i in range(2):
        x = day_nums[i]
        filter_list = []
        for chi_2 in chis[i]:
            if chi_2 < 10:
                filter_list.append(True)
            else:
                filter_list.append(False)
        for index,element in enumerate(x):
            if element < 0:
                filter_list[index] = False
            else:
                pass
        filter_lists[i] = filter_list
        x_array = np.linspace(np.amin(x[filter_list]), np.amax(x[filter_list]), 10000)

        popt, pcov = curve_fit(linear, x[filter_list], res[i][filter_list],
                               sigma=res_errs[i][filter_list],
                               p0=p_guess, bounds=p_bounds, maxfev=5000)

        chi_2 = chi2(linear(x[filter_list], *popt),
                     res_errs[i][filter_list],
                     res[i][filter_list],
                     2)

        print(popt, pcov, chi_2)

        fig, ax1 = plt.subplots()
        ax1.set_title("1MeV PMT Resolution Ch{}".format(i))
        ax1.grid(True)
        color = 'tab:red'
        plt.errorbar(x[filter_list],res[i][filter_list],
                     res_errs[i][filter_list],
                     fmt='.',
                     color=color)
        plt.plot(x_array, linear(x_array, *popt), 'g')
        if i ==0:
            plt.axvline(98, color='r', ls='--')
            plt.axvline(0, color='b', ls='--')
        ax1.text(0, 3, 'linear fit: y = ({:.4f}±{:.4f})*x + ({:.3f}±{:.3f}) chi2R: {:.3f}'.format(popt[0],
                                                                                                      np.sqrt(
                                                                                                          pcov[0][0]),
                                                                                                      popt[1],
                                                                                                      np.sqrt(
                                                                                                          pcov[1][1]),
                                                                                                      chi_2))

        ax1.set_xlabel("Exposure Time /days")
        ax1.set_ylabel("1Mev Resolution /%", color=color)
        ax1.set_ylim(2, 4.25)
        ax1.tick_params(axis='y', labelcolor=color)

        color = 'tab:blue'
        ax2 = ax1.twinx()
        ax2.plot(x[filter_list], chis[i][filter_list], '.', color=color)
        ax2.set_ylabel('Chi2 of bismuth 1MeV fit', color=color)  # we already handled the x-label with ax1
        ax2.tick_params(axis='y', labelcolor=color)
        ax2.set_ylim(0, 10)
        fig.tight_layout()

        plt.savefig("/Users/willquinn/Desktop/resolution_vs_date_time_" + str(i))
        plt.close()

    for i in range(2):
        x = exposures[i]
        filter_list = []
        for chi_2 in chis[i]:
            if chi_2 < 10:
                filter_list.append(True)
            else:
                filter_list.append(False)
        for index,element in enumerate(x):
            if element < 0.01:
                filter_list[index] = False
            else:
                pass
        filter_lists[i] = filter_list
        x_array = np.linspace(np.amin(x[filter_list]), np.amax(x[filter_list]), 10000)

        popt, pcov = curve_fit(linear, x[filter_list], res[i][filter_list],
                               sigma=res_errs[i][filter_list],
                               p0=p_guess, bounds=p_bounds, maxfev=5000)

        chi_2 = chi2(linear(x[filter_list], *popt),
                     res_errs[i][filter_list],
                     res[i][filter_list],
                     2)

        print(popt, pcov, chi_2)

        fig, ax1 = plt.subplots()
        ax1.set_title("1MeV PMT Resolution Ch{}".format(i))
        ax1.grid(True)
        color = 'tab:red'
        plt.errorbar(x[filter_list],res[i][filter_list],
                     res_errs[i][filter_list],
                     fmt='.',
                     color=color)
        plt.plot(x_array, linear(x_array, *popt), 'g')
        if i == 0:
            plt.axvline(1/10 + 97/100, color='r', ls='--')
            plt.axvline(29/1000000, color='b', ls='--')
        ax1.text(0, 3, 'linear fit: y = ({:.4f}±{:.4f})*x + ({:.3f}±{:.3f}) chi2R: {:.3f}'.format(popt[0],
                                                                                                      np.sqrt(
                                                                                                          pcov[0][0]),
                                                                                                      popt[1],
                                                                                                      np.sqrt(
                                                                                                          pcov[1][1]),
                                                                                                      chi_2))

        ax1.set_xlabel("Exposure Time /atm-days")
        ax1.set_ylabel("1Mev Resolution /%", color=color)
        ax1.set_ylim(2, 4.25)
        ax1.tick_params(axis='y', labelcolor=color)

        color = 'tab:blue'
        ax2 = ax1.twinx()
        ax2.plot(x[filter_list], chis[i][filter_list], '.', color=color)
        ax2.set_ylabel('Chi2 of bismuth 1MeV fit', color=color)  # we already handled the x-label with ax1
        ax2.tick_params(axis='y', labelcolor=color)
        ax2.set_ylim(0, 10)
        fig.tight_layout()

        plt.savefig("/Users/willquinn/Desktop/resolution_vs_date_exp_" + str(i))
        plt.close()

    x_t = []
    x_e = []
    y = []
    y_err = []
    for i in range(day_nums[0].size):
        for j in range(day_nums[1].size):
            if day_nums[0][i] == day_nums[1][j]:
                x_t.append(day_nums[0][i])
                x_e.append(exposures[0][i])
                y.append(res[0][i]/res[1][j])
                y_err.append(res[0][i]/res[1][j] * np.sqrt( (res_errs[0][i]/res[0][i])**2 + (res_errs[1][j]/res[1][j])**2 ))
                break

    plt.errorbar(x_t, y, yerr=y_err, fmt="k.")
    plt.axvline(98, color="r", ls="--")
    plt.axvline(0, color="b", ls="--")
    plt.xlabel("Exposure /days")
    plt.ylabel("Ratio res_Ch0/res_Ch1")
    plt.title("Ratio of resolution of CH 0 & 1 vs time")
    plt.grid()
    plt.xlim(np.amin(np.array(x_t)), np.amax(np.array(x_t)))
    plt.ylim(0,2)
    plt.savefig("/Users/willquinn/Desktop/resolution_ratio_time")


    '''x_ch0 = process_date(data_Ch0[1])
    x_ch1 = process_date(data_Ch1[1])
    x_array = np.linspace(-30, 145, 10000)'''

    '''plt.errorbar(x_ch0, data_Ch0[3], data_Ch0[4], fmt='.', color='b')
    plt.title("1MeV PMT Charge Ch0")
    plt.xlabel("Exposure Time atm-days")
    plt.ylabel("1Mev Charge /pC")
    plt.ylim(32, 42)
    #plt.xscale('log')
    plt.grid(True)
    plt.show()

    plt.errorbar(x_ch1, data_Ch1[3], data_Ch1[4], fmt='.', color='b')
    plt.title("1MeV PMT Charge Ch1")
    plt.xlabel("Exposure Time /days")
    plt.ylabel("1Mev Charge /pC")
    plt.ylim(32, 42)
    #plt.xscale('log')
    plt.grid(True)
    plt.show()'''

    '''popt_0, pcov_0 = curve_fit(linear, x_ch0, data_Ch0[7],
                               sigma=(data_Ch0[7] * np.sqrt(
                                   (data_Ch0[6] / data_Ch0[5]) ** 2 + (data_Ch0[4] / data_Ch0[3]) ** 2)),
                               p0=p_guess, maxfev=5000)
    popt_1, pcov_1 = curve_fit(linear, x_ch1, data_Ch1[7],
                               sigma=(data_Ch1[7] * np.sqrt(
                                   (data_Ch1[6] / data_Ch1[5]) ** 2 + (data_Ch1[4] / data_Ch1[3]) ** 2)),
                               p0=p_guess, maxfev=5000)

    chi2_1 = chi2(linear(x_ch1, *popt_1),
                  (data_Ch1[7] * np.sqrt((data_Ch1[6] / data_Ch1[5]) ** 2 + (data_Ch1[4] / data_Ch1[3]) ** 2)),
                  data_Ch1[7], 2)

    print(popt_0)
    print(popt_1)

    fig, ax1 = plt.subplots()
    ax1.set_title("1MeV PMT Resolution Ch0")
    ax1.grid(True)
    color = 'tab:red'
    plt.errorbar(x_ch0, data_Ch0[7],
                 data_Ch0[7] * np.sqrt((data_Ch0[6] / data_Ch0[5]) ** 2 + (data_Ch0[4] / data_Ch0[3]) ** 2), fmt='.',
                 color=color)
    plt.plot(x_array, linear(x_array, *popt_0), 'g')
    ax1.text(-20, 2.1, 'linear fit: y = ({:.4f}±{:.4f})*x + ({:.3f}±{:.3f}) chi2R: {:.3f}'.format(popt_0[0],
                                                                                                  np.sqrt(pcov_0[0][0]),
                                                                                                  popt_0[1],
                                                                                                  np.sqrt(pcov_0[1][1]),
                                                                                                  chi2_0))

    ax1.set_xlabel("Exposure Time atm-days")
    ax1.set_ylabel("1Mev Resolution /%", color=color)
    # ax1.set_ylim(2, 4.25)
    ax1.tick_params(axis='y', labelcolor=color)

    color = 'tab:blue'
    ax2 = ax1.twinx()
    ax2.plot(x_ch0, data_Ch0[7], '.', color=color)
    ax2.set_ylabel('Chi2', color=color)  # we already handled the x-label with ax1
    ax2.tick_params(axis='y', labelcolor=color)
    # ax2.set_ylim(0, 20)
    fig.tight_layout()

    plt.show()

    fig, ax1 = plt.subplots()
    ax1.set_title("1MeV PMT Resolution Ch1")
    ax1.grid(True)
    color = 'tab:red'
    plt.errorbar(x_ch1, data_Ch1[7],
                 data_Ch1[7] * np.sqrt((data_Ch1[6] / data_Ch1[5]) ** 2 + (data_Ch1[4] / data_Ch1[3]) ** 2), fmt='.',
                 color=color)
    plt.plot(x_array, linear(x_array, *popt_1), 'g')
    ax1.text(-20, 2.25,
             'linear fit: y = ({:.4f}±{:.4f})*x + ({:.3f}±{:.3f}) chi2R: {:.3f}'.format(popt_1[0],
                                                                                        np.sqrt(pcov_1[0][0]),
                                                                                        popt_1[1],
                                                                                        np.sqrt(pcov_1[1][1]),
                                                                                        chi2_1))

    ax1.set_xlabel("Exposure Time /days")
    ax1.set_ylabel("1Mev Resolution /%", color=color)
    # ax1.set_ylim(2, 4.25)
    ax1.tick_params(axis='y', labelcolor=color)

    color = 'tab:blue'
    ax2 = ax1.twinx()
    ax2.plot(x_ch1, data_Ch1[8], '.', color=color)
    ax2.set_ylabel('Chi2', color=color)  # we already handled the x-label with ax1
    ax2.tick_params(axis='y', labelcolor=color)
    # ax2.set_ylim(0, 20)
    fig.tight_layout()

    plt.show()'''


if __name__ == '__main__':
    main()
