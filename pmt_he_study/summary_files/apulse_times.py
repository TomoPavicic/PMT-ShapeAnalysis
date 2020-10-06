import sys

sys.path.insert(1, '../..')

import matplotlib.pyplot as plt
import numpy as np
import ROOT
import tqdm

from functions.other_functions import io_parse_arguments, process_date, process_exposure


names = ['GAO607', 'GAO612']


def main():
    plt.rcParams["figure.figsize"] = (20, 5)

    args = io_parse_arguments()
    directory = args.i
    output_filename = args.o

    try:
        filenames_file = open(directory + "filenames.txt", 'r')
    except FileNotFoundError as fnf_error:
        print(fnf_error)
        raise Exception("Error opening data file")

    file_names = np.loadtxt(directory + "filenames.txt", delimiter=',',
                            dtype={'names': ['filename'], 'formats': ['S100']}, unpack=True)

    success_rate = [0, 0]
    num_average = 7
    re_bin = 10
    number_total_weeks = 52
    last_temp = [[],[]]
    x_og = []
    num_in_week = [[0 for i in range(number_total_weeks+1)],[0 for i in range(number_total_weeks+1)]]
    apulse_num_y = [[[] for i in range(10)],[[] for i in range(10)]]
    day_nums = [[],[]]

    color = 'g-'
    #for i_file in tqdm.tqdm(range(file_names.size)):
    for i_file in range(file_names.size):
        #print(directory + "/" + file_names[i_file][0].decode("utf-8"))
        file = ROOT.TFile(directory + file_names[i_file][0].decode("utf-8"), "READ")
        date = file_names[i_file][0].decode("utf-8").split("_")[0]
        exposure = process_exposure(np.array([int(date)]))[0]
        day_num = process_date(np.array([int(date)]))[0]
        week_num = int(day_num/num_average)
        time = file_names[i_file][0].decode("utf-8").split("_")[1]

        if i_file == 0:
            name = date + "_" + 'GAO607_apulse_times_1400V'
            hist = file.Get(name)
            temp = []
            the_bin = 0

            #print(hist.GetNbinsX())

            for i_bin in range(hist.GetNbinsX()):
                if the_bin == 0:
                    temp0 = 0

                temp0 += hist.GetBinContent(i_bin)

                the_bin += 1

                if the_bin == re_bin:
                    temp.append(temp0)
                    x_og.append(i_bin * hist.GetBinWidth(i_bin))
                    the_bin = 0

            for i in range(number_total_weeks+1):
                last_temp[0].append(np.zeros_like(np.array(temp)))
                last_temp[1].append(np.zeros_like(np.array(temp)))


        if day_num < 0:
            continue

        for i_channel in range(2):
            #print(date + "_" + names[i_channel] + "_apulse_times_1400")
            hist = file.Get(date + "_" + names[i_channel] + "_apulse_times_1400V")
            hist_0 = file.Get(date + "_" + names[i_channel] + "_apulse_num_1400V")

            try:
                hist.GetNbinsX()
            except:
                del hist
                continue

            temp = []
            the_bin = 0
            x = []

            for i_bin in range(hist.GetNbinsX()):
                if the_bin == 0:
                    temp0 = 0

                temp0 += hist.GetBinContent(i_bin)

                the_bin += 1

                if the_bin == re_bin:
                    temp.append(temp0)
                    x.append(i_bin*hist.GetBinWidth(i_bin))
                    the_bin = 0
            if i_channel == 0:
                if 0 < day_num < 98:
                    color = 'b-'
                if day_num >= 98:
                    color = 'r-'
                plt.plot(x, 100 * np.array(temp)/hist_0.GetEntries(), color)
                plt.plot(0, 0, 'g-', label="Atmosphere He")
                plt.plot(0, 0, 'b-', label="1% He")
                plt.plot(0, 0, 'r-', label="10% He")
                plt.axvline(1580, ls='--', color='r', label=r"$He^+$")
                plt.axvline(2680, ls='--', color='g', label=r"$H_2O^+$")
                plt.axvline(1290, ls='--', color='k', label=r"$H_2^+ He^{2+}$")
                plt.axvline(2560, ls='--', color='b', label=r"$CH_4^+$")
                plt.axvline(3850, ls='--', color='c', label=r"$CO_2^+$")
                plt.legend(loc="upper right")
                plt.ylim(0, 6)
                plt.xlim(800, 7000)
                plt.xlabel("afterpulse time in waveform /ns")
                plt.ylabel("normalised counts /%")
                plt.title("Ch{} Date:{} Exposure:{} atm-day DayNUM:{}".format(i_channel, date, round(exposure,2), day_num))
                plt.grid()
                #plt.yscale('log')
                plt.show(block=False)
                plt.pause(0.01)
                plt.close()

            last_temp[i_channel][week_num] += 100 * np.array(temp)/hist.GetEntries()
            num_in_week[i_channel][week_num] += 1
            success_rate[i_channel] += 1
            del hist

        for i_channel in range(2):
            hist = file.Get(date + "_" + names[i_channel] + "_apulse_num_1400V")

            try:
                hist.GetNbinsX()
            except:
                del hist
                continue

            for i_bin in range(len(apulse_num_y[i_channel])):
                if hist.GetEntries() > 0:
                    apulse_num_y[i_channel][i_bin].append(hist.GetBinContent(i_bin)/hist.GetEntries())
                else:
                    pass

            day_nums[i_channel].append(day_num)
            del hist

        file.Close()
        del file

    for i_week in range(len(last_temp[0])):
        if i_week < 14:
            color = 'b-'
        elif i_week >= 14:
            color = 'r-'
        else:
            color = 'g-'
        plt.plot(x_og, last_temp[0][i_week]/num_in_week[0][i_week], color)
        plt.title("Ch:{} Week num:{}".format(0, i_week))
        plt.ylim(0,6)
        plt.plot(0, 0, 'g-', label="Atmosphere He")
        plt.plot(0, 0, 'b-', label="1% He")
        plt.plot(0, 0, 'r-', label="10% He")
        plt.axvline(1580, ls='--', color='r', label=r"$He^+$")
        plt.axvline(2680, ls='--', color='g', label=r"$H_2O^+$")
        plt.axvline(1290, ls='--', color='k', label=r"$H_2^+ He^{2+}$")
        plt.axvline(2560, ls='--', color='b', label=r"$CH_4^+$")
        plt.axvline(3850, ls='--', color='c', label=r"$CO_2^+$")
        plt.xlim(800, 7000)
        plt.xlabel("afterpulse time in waveform /ns")
        plt.ylabel("normalised counts /%")
        plt.legend(loc="upper right")
        plt.grid()
        plt.show(block=False)
        plt.pause(0.1)
        plt.close()
    #plt.rcParams["figure.figsize"] = (20, 20)'''
    '''for i_bin in range(1,len(apulse_num_y[0])):
        ylim = 0
        if i_bin == 1:
            ylim = 100
        else:
            ylim = 100/ (i_bin*2)

        plt.subplot(3,1,1)
        plt.plot(day_nums[0],100*np.array(apulse_num_y[0][i_bin]), 'k.')
        plt.axvline(98, color='r', ls='--')
        plt.axvline(0, color='b', ls='--')
        plt.title("Apulse num:{}".format(i_bin-1))
        plt.grid()
        plt.ylim(0,ylim)

        plt.subplot(3,1,2)
        plt.plot(day_nums[1],100*np.array(apulse_num_y[1][i_bin]), 'k.')
        plt.axvline(98, color='r', ls='--')
        plt.axvline(0, color='b', ls='--')
        plt.ylabel("Percentage apulses /%")
        plt.grid()
        plt.ylim(0, ylim)

        plt.subplot(3,1,3)
        plt.plot(day_nums[0], np.array(apulse_num_y[0][i_bin])/np.array(apulse_num_y[1][i_bin]), 'k.')
        plt.axvline(98, color='r', ls='--')
        plt.axvline(0, color='b', ls='--')
        plt.xlabel("Exposure time /days")
        plt.grid()
        plt.ylabel("Ratio")
        plt.savefig('/Users/willquinn/Desktop/apulse_num_bin_'+str(i_bin))
        plt.close()'''

    print("Success rate: {} Ch0 {} Ch1".format(success_rate[0] / file_names.size, success_rate[1] / file_names.size))

    return


if __name__ == '__main__':
    main()
