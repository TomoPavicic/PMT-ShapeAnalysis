import sys
sys.path.insert(1, '..')

import matplotlib.pyplot as plt
import numpy as np
from functions.other_functions import io_parse_arguments, process_date, process_exposure, linear, chi2
from scipy.optimize import curve_fit


def main():
    args = io_parse_arguments()
    data_file_Ch0_name = args.i + "apulseNUM_vs_date_Ch0.txt"
    data_file_Ch1_name = args.i + "apulseNUM_vs_date_Ch1.txt"

    try:
        data_file_Ch0 = open(data_file_Ch0_name, 'r')
        data_file_Ch1 = open(data_file_Ch1_name, 'r')
    except FileNotFoundError as fnf_error:
        print(fnf_error)
        raise Exception("Error opening data file")

    apulse_data_Ch0 = np.loadtxt(data_file_Ch0_name, delimiter=',', dtype={'names': ('voltage', 'unix', 'date', 'timestamp', 'perc_l', 'std_l', 'perc_m', 'std_m', 'perc_h', 'std_h'),
                                 'formats': ('i4', 'i4', 'i4', 'i4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4')}, unpack=True)

    apulse_data_Ch1 = np.loadtxt(data_file_Ch1_name, delimiter=',', dtype={'names': ('voltage', 'unix', 'date', 'timestamp', 'perc_l', 'std_l', 'perc_m', 'std_m', 'perc_h', 'std_h'),
                                 'formats': ('i4', 'i4', 'i4', 'i4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4')}, unpack=True)

    print(np.sort(apulse_data_Ch0[2]))
    x_ch0 = process_date(apulse_data_Ch0[2])
    x_ch0_1 = process_exposure(apulse_data_Ch0[2])
    x_ch1 = process_date(apulse_data_Ch1[2])
    x_array = np.linspace(-30, 145, 10000)

    #plt.errorbar(x_ch0, apulse_data_Ch0[3], apulse_data_Ch0[4], fmt='.', color='b')
    plt.plot(x_ch0, apulse_data_Ch0[4], '.', color='r', label='low cut 0.85 5')
    plt.plot(x_ch0, apulse_data_Ch0[6], '.', color='g', label='med cut 0.9 10')
    plt.plot(x_ch0, apulse_data_Ch0[8], '.', color='b', label='hgh cut 0.95 15')
    plt.title("Average number of afterpulses per waveform ch0")
    plt.xlabel("Exposure Time/days")
    plt.ylabel("Percentage number of afterpulses /%")
    plt.axvline(98, color='r', ls='--')
    plt.axvline(0, color='b', ls='--')
    plt.legend()
    plt.ylim(0, 100)
    #plt.xscale('log')
    plt.grid(True)
    plt.savefig("/Users/willquinn/Desktop/apulse_vs_date_time_ch0")
    plt.close()

    plt.plot(x_ch0_1, apulse_data_Ch0[4], '.', color='r', label='low cut 0.85 5')
    plt.plot(x_ch0_1, apulse_data_Ch0[6], '.', color='g', label='med cut 0.9 10')
    plt.plot(x_ch0_1, apulse_data_Ch0[8], '.', color='b', label='hgh cut 0.95 15')
    plt.axvline(1 / 10 + 97 / 100, color='r', ls='--')
    plt.axvline(29 / 1000000, color='b', ls='--')
    plt.title("Average number of afterpulses per waveform ch0")
    plt.xlabel("Exposure Time/atm-days")
    plt.ylabel("Percentage number of afterpulses /%")
    plt.legend()
    plt.ylim(0, 100)
    # plt.xscale('log')
    plt.grid(True)
    plt.savefig("/Users/willquinn/Desktop/apulse_vs_date_exp_ch0")
    plt.close()

    # plt.errorbar(x_ch0, apulse_data_Ch0[3], apulse_data_Ch0[4], fmt='.', color='b')
    plt.plot(x_ch1, apulse_data_Ch1[4], '.', color='r', label='low cut 0.85 5')
    plt.plot(x_ch1, apulse_data_Ch1[6], '.', color='g', label='med cut 0.9 10')
    plt.plot(x_ch1, apulse_data_Ch1[8], '.', color='b', label='hgh cut 0.95 15')
    plt.title("Average number of afterpulses per waveform ch1")
    plt.xlabel("Exposure Time /days")
    plt.ylabel("Percentage number of afterpulses /%")
    plt.legend()
    plt.ylim(0, 100)
    # plt.xscale('log')
    plt.grid(True)
    plt.savefig("/Users/willquinn/Desktop/apulse_vs_date_time_ch1")
    plt.close()

    plt.plot(x_ch1, apulse_data_Ch0[8]/apulse_data_Ch1[8], '.', color='b', label='high cut')
    plt.title("Ratio of average afterpulse number CH 0 & 1")
    plt.xlabel("Exposure Time /days")
    plt.ylabel("ratio of avg apule rate ch0/ch1 ")
    plt.axvline(0, color='k', ls='--')
    plt.axvline(98, color='k', ls='--')
    # plt.legend()
    # plt.ylim(20, 60)
    # plt.xscale('log')
    plt.grid(True)
    plt.savefig("/Users/willquinn/Desktop/apulse_vs_date_ratio")
    plt.close()


if __name__ == '__main__':
    main()
