import sys

sys.path.insert(1, '..')

# import ROOT and bash commands
import ROOT
import tqdm

# import python plotting and numpy modules
import matplotlib.pyplot as plt
import numpy as np
import random as rand

# import stats module
from scipy.optimize import curve_fit

# import custom made classes
from functions.other_functions import sncd_parse_arguments, get_date_time, get_voltage, fit, chi2
from src.PMT_Array import PMT_Array
from src.PMT_Waveform import PMT_Waveform


def main():

    # Handle the input arguments:
    ##############################
    args = sncd_parse_arguments()
    input_data_file_name = args.i
    config_file_name = args.c
    output_file_name = args.o
    sweep_bool = "True"
    ##############################

    # Do some string manipulation to get the date and time from the file name
    #################################################
    date, time = "0000", "0000"
    voltage = 1400
    #################################################

    # Check to see if the data file exists
    try:
        print(">>> Reading data from file: {}".format(input_data_file_name))
        date_file = open(input_data_file_name, 'r')
    except FileNotFoundError as fnf_error:
        print(fnf_error)
        raise Exception("Error opening data file {}".format(input_data_file_name))

    file = ROOT.TFile(input_data_file_name, "READ")
    output_file = ROOT.TFile(output_file_name, "RECREATE")
    file.cd()

    tree = file.T
    tree.Print()

    num_events = 10000

    topology = [4, 1]
    pmt_array = PMT_Array(topology, date + "_" + time)
    pmt_array.set_pmt_id("GAO607_" + date + "_" + time, 0)
    pmt_array.set_pmt_id("GAO612_" + date + "_" + time, 1)
    pmt_array.set_pmt_id("Injected_GAO607_" + date + "_" + time, 2)
    pmt_array.set_pmt_id("Injected_GAO612_" + date + "_" + time, 3)

    # Set the cuts you wish to apply
    # If you don't do this the defaults are used
    if config_file_name is not None:
        pmt_array.apply_setting(config_file_name)

    if sweep_bool == "True":
        pmt_array.set_sweep_bool(True)
        if voltage == 1000:
            pmt_array.set_pmt_templates(
                "~/Desktop/191008_A1000_B1000_templates.root",
                ["A1000_B1000_Ch0_Template", "A1000_B1000_Ch1_Template", "A1000_B1000_Ch0_Template", "A1400_B1400_Ch1_Template"])
        elif voltage == 1400:
            pmt_array.set_pmt_templates(
                "~/Desktop/190621_A1400_B1400_templates.root",
                ["A1400_B1400_Ch0_Template", "A1400_B1400_Ch1_Template", "A1400_B1400_Ch0_Template", "A1400_B1400_Ch1_Template"])

    amp_cut = pmt_array.get_pmt_object_number(0).get_setting("mf_amp_threshold")
    shape_cut = pmt_array.get_pmt_object_number(0).get_setting("mf_shape_threshold")

    t_injected = [[], []]
    a_injected = [[], []]
    a_injected_failure = [[], []]
    t_injected_success = [[], []]
    t_injected_failure = [[], []]
    mf_a_injected = [[], []]
    mf_s_injected = [[], []]

    num = [0, 0]
    enum = [0, 0]
    iterator = 0

    ran_range = 52

    injected_array = [[[], []], [[], []]]
    best_array = [[[], []], [[], []]]
    random_array = [[[], []], [[], []]]
    x_best_array = [[], []]

    for event in tqdm.tqdm(tree):
        OM_ID = event.OM_ID
        waveform = np.array(event.waveform)

        # Get template for injecting into data
        template_pulse = pmt_array.get_pmt_object_number(OM_ID).get_template_pmt_pulse()

        # Get the amplitude of the template so we can normalise to 1
        temp_amplitude = np.amin(template_pulse)

        # Create random numbers
        random_amp = rand.randrange(0, ran_range - 1, 2)
        random_time = rand.randrange(800, pmt_array.get_pmt_object_number(OM_ID).get_waveform_length() - template_pulse.size)

        '''print("A: ",random_amp)
        print("t: ",random_time)'''

        factor = random_amp / temp_amplitude

        injected_data = np.zeros_like(waveform)

        j = 0
        for i in range(injected_data.size):
            if i >= random_time:
                injected_data[i] = int(template_pulse[j] * factor)
                j += 1
                if j == template_pulse.size:
                    break

        pmt_waveform = PMT_Waveform(waveform, pmt_array.get_pmt_object_number(OM_ID))
        pmt_waveform_injected = PMT_Waveform(waveform - injected_data, pmt_array.get_pmt_object_number(OM_ID + 2))
        plt.plot(waveform - injected_data, 'r-', label="injected")
        plt.plot(waveform, 'b-', label="raw data")
        plt.xlabel("time /ns")
        plt.ylabel("ADC unit /mV")
        plt.legend(loc="lower right")
        plt.grid(True)
        plt.show()
        check = False

        if pmt_waveform.get_pmt_pulse_trigger():
            i_shape = 0
            i_amp = 0
            r_shape = 0
            r_amp = 0

            for index, value in enumerate(pmt_waveform_injected.get_pmt_pulse_times()):
                if random_time == value + 800:
                    num[OM_ID] += 1
                    i_shape = pmt_waveform_injected.pmt_waveform_sweep_shape[value]
                    i_amp = pmt_waveform_injected.pmt_waveform_sweep_amp[value]
                    r_shape = pmt_waveform.pmt_waveform_sweep_shape[value]
                    r_amp = pmt_waveform.pmt_waveform_sweep_amp[value]
                    check = True

            # Store random numbers
            t_injected[OM_ID].append(random_time)

            if check:
                a_injected[OM_ID].append(random_amp)
                t_injected_success[OM_ID].append(random_time)
                mf_a_injected[OM_ID].append(i_amp)
                mf_s_injected[OM_ID].append(i_shape)
            else:
                a_injected_failure[OM_ID].append(random_amp)
                t_injected_failure[OM_ID].append(random_time)

            injected_array[OM_ID][0].append(i_shape)
            injected_array[OM_ID][1].append(i_amp)
            random_array[OM_ID][0].append(r_shape)
            random_array[OM_ID][1].append(r_amp)

            best_shape = np.amax(pmt_waveform.pmt_waveform_sweep_shape)
            i_best = np.argmax(pmt_waveform.pmt_waveform_sweep_shape)
            best_amplitude = pmt_waveform.pmt_waveform_sweep_amp[i_best]

            x_best_array[OM_ID].append(i_best)
            best_array[OM_ID][0].append(best_shape)
            best_array[OM_ID][1].append(best_amplitude)

            '''fig, ax1 = plt.subplots()

            color = 'tab:red'
            ax1.set_xlabel('timestamp (ns)')
            ax1.set_ylabel('shape', color=color)
            ax1.plot(pmt_waveform.pmt_waveform_sweep_amp, color=color)
            ax1.plot(pmt_waveform.get_pmt_pulse_times(), pmt_waveform.pmt_waveform_sweep_amp[pmt_waveform.get_pmt_pulse_times()], "x", color='tab:green')
            ax1.tick_params(axis='y', labelcolor=color)
            ax1.plot(np.zeros_like(pmt_waveform.pmt_waveform_sweep_shape), "--", color="gray")
            plt.show(block=False)
            plt.pause(0.5)
            plt.close()'''

            enum[OM_ID] += 1
        iterator += 1
        if iterator == num_events:
            break
        #bar.update(iterator)

    #bar.finish()

    print("Percentage injected found: ",num[0]/enum[0], "% ch0", num[1]/enum[1], "% ch1")

    for i in range(2):
        plt.plot(t_injected[i], injected_array[i][0], 'r.', label="injected")
        plt.plot(t_injected[i], random_array[i][0], 'g.', label="random")
        plt.plot(t_injected[i], best_array[i][0], 'b.', label="best")
        plt.xlabel("time /ns")
        plt.ylabel("shape index")
        plt.title("Channel: {} events: {}".format(i, num_events))
        plt.grid(True)
        plt.legend(loc="lower right")
        plt.axhline(y=shape_cut, color='k', linestyle='-')
        plt.savefig("/Users/willquinn/Desktop/pmt_sim_results/shape_index_vs_time_ch{}".format(i))
        plt.show()

        plt.hist(injected_array[i][0], bins=100, range=(0,1), alpha=0.5, label="injected")
        plt.hist(random_array[i][0], bins=100, range=(0,1), alpha=0.5, label="random")
        plt.hist(best_array[i][0], bins=100, range=(0,1), alpha=0.5, label="best")
        plt.xlabel("shape index")
        plt.ylabel("counts")
        plt.title("Channel: {} events: {}".format(i, num_events))
        plt.grid(True)
        plt.yscale('log')
        plt.legend(loc="upper left")
        plt.axvline(x=shape_cut, color='r', linestyle='-')
        plt.savefig("/Users/willquinn/Desktop/pmt_sim_results/shape_index_vs_time_hist_ch{}".format(i))
        plt.show()

        plt.plot(t_injected_failure[i], a_injected_failure[i], 'r.', label="failures")
        plt.plot(t_injected_success[i], a_injected[i], 'g.', label="success")
        plt.xlabel("time /ns")
        plt.title("Channel: {} events: {}".format(i, num_events))
        plt.ylabel("amplitude injected")
        plt.grid(True)
        plt.legend(loc="lower right")
        plt.savefig("/Users/willquinn/Desktop/pmt_sim_results/amplitude_success_vs_failures_ch{}".format(i))
        plt.show()

        plt.hist(a_injected_failure[i], bins=int(ran_range/2), range=(0, ran_range), alpha=0.5, label="failures")
        plt.hist(a_injected[i], bins=int(ran_range/2), range=(0, ran_range), alpha=0.5, label="success")
        plt.xlabel("amplitude injected")
        plt.title("Channel: {} events: {}".format(i, num_events))
        plt.ylabel("counts")
        plt.grid(True)
        plt.legend(loc="upper left")
        plt.text(0, 0, 'Success percentage: {}%'.format(round(num[0]/enum[0] * 100), 4))
        plt.savefig("/Users/willquinn/Desktop/pmt_sim_results/amplitude_success_vs_failures_hist_ch{}".format(i))
        plt.show()

        new_mf_a_y = []
        new_mf_a_y_err = []
        new_mf_a_x = np.histogram(a_injected[i], bins=int(ran_range/2), range=(0,ran_range))
        x = []
        for j in range(new_mf_a_x[0].size):

            temp = []
            for k in range(len(a_injected[i])):
                if a_injected[i][k] == j * (ran_range / new_mf_a_x[0].size):
                    temp.append(mf_a_injected[i][k])
            x.append(j * (ran_range / new_mf_a_x[0].size))
            new_mf_a_y.append(np.average(temp))
            new_mf_a_y_err.append(np.std(temp))

        plt.plot(a_injected[i], mf_a_injected[i], 'b.')
        plt.errorbar(x, new_mf_a_y, yerr=new_mf_a_y_err, fmt='ro')
        plt.xlabel("random amplitude injected /mV")
        plt.title("Channel: {} events: {}".format(i, num_events))
        plt.ylabel("amplitude index")
        plt.axhline(y=amp_cut, color='g', linestyle='-')
        plt.grid(True)
        plt.savefig("/Users/willquinn/Desktop/pmt_sim_results/inj_vs_mf_amplitudes_ch{}".format(i))
        plt.show()

        new_mf_s_y = []
        new_mf_s_y_err = []
        new_mf_s_x = np.histogram(a_injected[i], bins=int(ran_range / 2), range=(0, ran_range))
        x = []
        for j in range(new_mf_s_x[0].size):
            temp = []
            for k in range(len(a_injected[i])):
                if a_injected[i][k] == j * (ran_range / new_mf_a_x[0].size):
                    temp.append(mf_s_injected[i][k])
            x.append(j * (ran_range / new_mf_a_x[0].size))
            new_mf_s_y.append(np.average(temp))
            new_mf_s_y_err.append(np.std(temp))

        plt.plot(a_injected[i], mf_s_injected[i], 'b.')
        plt.errorbar(x, new_mf_s_y, yerr=new_mf_s_y_err, fmt='ro')
        plt.xlabel("random amplitude injected /mV")
        plt.title("Channel: {} events: {}".format(i, num_events))
        plt.ylabel("shape index")
        plt.grid(True)
        plt.savefig("/Users/willquinn/Desktop/pmt_sim_results/inj_amp_vs_shape_ch{}".format(i))
        plt.show()

    output_file.Close()


if __name__ == '__main__':
    main()
