"""
    This class performs all the analysis on a PMT waveform
    Input the relevent PMT and waveform
        - PMT Object ; PMT_Object
        - Waveform ; waveform_array
"""

import numpy as np
from scipy.signal import find_peaks
from src.PMT_Object import PMT_Object
import matplotlib.pyplot as plt
import ROOT


class PMT_Waveform:

    def __init__(self, pmt_waveform_list: list, pmt_object: PMT_Object):

        # Store the PMT_Object into memory
        self.pmt_object = pmt_object

        self.pmt_trace_id = self.pmt_object.get_data_id() + "_ch" + self.pmt_object.get_pmt_id() + "_tr" + str(self.pmt_object.get_event_number())

        # If the input waveform list is of strings we need to check each element
        if type(pmt_waveform_list[0]) == str:
            temp_list_1 = []
            for i in range(len(pmt_waveform_list)):
                if pmt_waveform_list[i] == '' or pmt_waveform_list[i] == '\n':
                    pass
                else:
                    temp_list_1.append(pmt_waveform_list[i].strip())
        # If not a string we don't need to worry
        else:
            temp_list_1 = pmt_waveform_list

        # Store waveform
        self.pmt_waveform = np.array(temp_list_1, dtype='float')

        # Check whether the pulse lies where it is expected i.e above the 100 ns timestamp
        self.pmt_pulse_peak_amplitude = 0.0
        self.pmt_pulse_peak_position = np.argmin(self.pmt_waveform)

        self.pmt_pulse_start = 0
        self.pmt_pulse_end = 0
        self.pmt_baseline = 0.0
        self.pmt_pulse_charge = 0.0
        self.pmt_apulse_charge = 0.0
        self.pmt_pulse_mf_shape = 0.0
        self.pmt_pulse_mf_amp = 0.0
        self.pmt_pulse_trigger = False
        self.pmt_apulse_trigger = False
        self.pmt_pulse_times = []
        self.pmt_waveform_sweep_shape = []
        self.pmt_waveform_sweep_amp = []
        self.pmt_waveform_reduced = np.array([], dtype='float')

        self.results_dict = {}
        if self.get_pmt_pulse_peak_position() < self.get_pmt_object().get_pulse_time_threshold():
            pass
        else:

            self.pmt_waveform_length = self.get_pmt_waveform().size
            self.set_pmt_baseline(np.average(self.pmt_waveform[0:self.get_pmt_object().get_setting("trigger_point")]))
            self.set_pmt_waveform_reduced()
            self.set_pmt_pulse_peak_amplitude(-1 * np.amin(self.get_pmt_waveform_reduced()))
            self.calculate_charge()

            if self.get_pmt_pulse_charge() < self.get_pmt_object().get_setting("charge_cut"):
                self.update_results_dict()
                return
            else:
                self.set_pmt_pulse_trigger(True)

                # Only sweep the waveform if there is a template
                if self.get_pmt_object().get_template_bool() and self.get_pmt_pulse_peak_position() < self.get_pmt_object().get_waveform_length() - self.get_pmt_object().get_template_pmt_pulse().size:

                    pmt_pulse = self.get_pmt_waveform()[self.get_pmt_pulse_peak_position() - np.argmin(self.get_pmt_object().get_template_pmt_pulse()): self.get_pmt_pulse_peak_position() - np.argmin(self.get_pmt_object().get_template_pmt_pulse()) + self.get_pmt_object().get_template_pmt_pulse().size] - self.get_pmt_baseline()

                    inner_product = np.dot(self.pmt_object.get_template_pmt_pulse(), pmt_pulse)
                    self.get_pmt_object().get_normalisation_factor(pmt_pulse)

                    self.set_pmt_pulse_mf_amp(inner_product)
                    self.set_pmt_pulse_mf_shape(inner_product/self.get_pmt_object().get_normalisation_factor(pmt_pulse))

                    if self.get_pmt_object().get_sweep_bool():
                        self.pmt_pulse_sweep()

        self.update_results_dict()
        # Store the pmt waveform length
        self.pmt_waveform_length = self.pmt_waveform.size
        # self.fill_pmt_hists()
        return

    def get_pulse_trigger(self):
        return self.pmt_pulse_trigger

    def get_pmt_object(self):
        return self.pmt_object

    def get_pmt_pulse_peak_amplitude(self):
        return self.pmt_pulse_peak_amplitude

    def set_pmt_pulse_peak_amplitude(self, new_amplitude: float):
        self.pmt_pulse_peak_amplitude = new_amplitude

    def get_pmt_apulse_charge(self):
        return self.pmt_apulse_charge

    def set_pmt_apulse_charge(self, new_apulse_charge: float):
        self.pmt_apulse_charge = new_apulse_charge

    def get_pmt_pulse_start(self):
        return self.pmt_pulse_start

    def get_results_dict(self):
        return self.results_dict

    def set_results_dict(self, description: str, value):
        self.get_results_dict()[description] = value

    def set_pmt_pulse_start(self, new_pmt_pulse_start: int):
        self.pmt_pulse_start = new_pmt_pulse_start

    def get_pmt_pulse_end(self):
        return self.pmt_pulse_end

    def set_pmt_pulse_end(self, new_pmt_pulse_end: int):
        self.pmt_pulse_end = new_pmt_pulse_end

    def get_pmt_pulse_trigger(self):
        return self.pmt_pulse_trigger

    def set_pmt_pulse_trigger(self, new_pmt_pulse_trigger: bool):
        self.pmt_pulse_trigger = new_pmt_pulse_trigger

    def get_pmt_waveform_reduced(self):
        return self.pmt_waveform_reduced

    def set_pmt_waveform_reduced(self):
        self.pmt_waveform_reduced = self.get_pmt_waveform() - self.get_pmt_baseline()

    def get_pmt_baseline(self):
        return self.pmt_baseline

    def set_pmt_baseline(self, new_pmt_baseline: float):
        self.pmt_baseline = new_pmt_baseline

    def get_pmt_pulse_charge(self):
        return self.pmt_pulse_charge

    def set_pmt_pulse_charge(self, new_pmt_pulse_charge: float):
        self.pmt_pulse_charge = new_pmt_pulse_charge

    def get_pmt_pulse_peak_position(self):
        return self.pmt_pulse_peak_position

    def get_pmt_waveform(self):
        return self.pmt_waveform

    def get_pmt_apulse_trigger(self):
        return self.pmt_apulse_trigger

    def set_pmt_apulse_trigger(self, new_bool: bool):
        self.pmt_apulse_trigger = new_bool

    def get_pmt_pulse(self):
        return self.pmt_waveform[self.get_pmt_pulse_start():self.get_pmt_pulse_end()] - self.pmt_baseline

    def get_pmt_waveform_length(self):
        return self.pmt_waveform_length

    def get_pmt_trace_id(self):
        return self.pmt_trace_id

    def fill_pmt_hists(self):
        #print(self.get_results_dict())
        self.pmt_object.fill_pmt_hists(self.get_results_dict())

    def update_results_dict(self):
        self.results_dict = {
            "pulse_charge": self.get_pmt_pulse_charge(),
            "pulse_amplitude": self.get_pmt_pulse_peak_amplitude(),
            "apulse_charge": self.get_pmt_apulse_charge(),
            "pulse_mf_shape": self.get_pmt_pulse_mf_shape(),
            "pulse_mf_amp": self.get_pmt_pulse_mf_amp(),
            "pulse_peak_time": self.get_pmt_pulse_peak_position(),
            "pulse_times": self.get_pmt_pulse_times(),
            "baseline": self.get_pmt_baseline()
        }

    def save_pmt_waveform_histogram(self, root_file: ROOT.TFile):
        name = self.pmt_object.get_pmt_id() + self.get_pmt_trace_id()
        pmt_waveform_hist = ROOT.TH1I(name, name, self.get_pmt_waveform_length(), 0, self.get_pmt_waveform_length())
        for timestamp in range(self.get_pmt_waveform_length()):
            pmt_waveform_hist.SetBinContent(timestamp, self.get_pmt_waveform()[timestamp])
        root_file.cd()
        pmt_waveform_hist.Write()
        del pmt_waveform_hist

    def check_cuts(self):
        # This function should be used in the init function and in the sweep function
        pass

    def pmt_pulse_sweep(self):
        sweep_start = self.get_pmt_object().get_sweep_range()[0]
        sweep_end = self.get_pmt_object().get_sweep_range()[1]
        sweep_window_length = self.pmt_object.get_template_pmt_pulse().size

        # plt.plot(self.get_pmt_waveform_reduced())
        # plt.xlabel('timestamp (ns)')
        # plt.title("Waveform {}".format(self.get_pmt_trace_id()))
        # plt.show(block=True)

        matched_filter_shape_list = []
        matched_filter_amplitude_list = []
        for i_sweep in range(sweep_start, sweep_end - sweep_window_length):
            pmt_waveform_section = self.get_pmt_waveform_reduced()[
                                   i_sweep : i_sweep + sweep_window_length]

            if pmt_waveform_section.size == self.pmt_object.get_template_pmt_pulse().size:
                pass
            else:
                print("Section size {} template size {}".format(pmt_waveform_section.size,
                                                                self.pmt_object.get_template_pmt_pulse().size))
                continue

            inner_product = np.dot(self.pmt_object.get_template_pmt_pulse(), pmt_waveform_section)

            matched_filter_amplitude_list.append(inner_product)
            matched_filter_shape_list.append(inner_product / np.sqrt(np.dot(pmt_waveform_section, pmt_waveform_section)))

            '''fig, ax1 = plt.subplots()
            color = 'tab:red'
            ax1.plot(pmt_waveform_section, color=color)
            ax1.set_xlabel('timestamp (ns)')
            ax1.set_ylabel('ADC /mV', color=color)
            ax1.tick_params(axis='y', labelcolor=color)
            color = 'tab:blue'
            ax2 = ax1.twinx()
            ax2.plot(self.pmt_object.get_template_pmt_pulse(), color=color)
            ax2.set_ylabel('Normalised ADC', color=color)  # we already handled the x-label with ax1
            ax2.tick_params(axis='y', labelcolor=color)
            fig.tight_layout()
            plt.text(0, 0, "Shape: {}".format(inner_product / np.sqrt(np.dot(pmt_waveform_section, pmt_waveform_section))))
            plt.show(block=False)
            plt.pause(0.01)
            plt.close()'''

        matched_filter_shape = np.array(matched_filter_shape_list)
        matched_filter_amplitude = np.array(matched_filter_amplitude_list)

        # TODO: check if both thresholds have been breached
        shape_peaks, _ = find_peaks(matched_filter_shape, height=self.get_pmt_object().get_setting("mf_shape_threshold"), distance=int(sweep_window_length / 2))
        amplitude_peaks, _ = find_peaks(matched_filter_amplitude, height=self.get_pmt_object().get_setting("mf_amp_threshold"), distance=int(sweep_window_length / 2))

        temp = []
        if len(shape_peaks) > 0:
            for index, value in enumerate(shape_peaks):
                if value in amplitude_peaks:
                    temp.append(value)

        if len(temp) > 0:
            self.set_pmt_apulse_trigger(True)
            self.set_pmt_pulse_times(temp)

        self.pmt_waveform_sweep_shape = matched_filter_shape
        self.pmt_waveform_sweep_amp = matched_filter_amplitude
        self.update_results_dict()

        '''fig, ax1 = plt.subplots()

        color = 'tab:red'
        ax1.set_xlabel('timestamp (ns)')
        ax1.set_ylabel('shape', color=color)
        ax1.plot(matched_filter_shape, color=color)
        ax1.plot(shape_peaks, matched_filter_shape[shape_peaks], "x", color='tab:green')
        ax1.tick_params(axis='y', labelcolor=color)
        ax1.plot(np.zeros_like(matched_filter_shape), "--", color="gray")

        ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

        color = 'tab:blue'
        ax2.set_ylabel('Waveform ADC /mV', color=color)  # we already handled the x-label with ax1
        ax2.plot(self.get_pmt_waveform_reduced(), color=color)
        ax2.tick_params(axis='y', labelcolor=color)

        ax2.set_ylabel('Amplitude', color=color)  # we already handled the x-label with ax1
        ax2.plot(matched_filter_amplitude, color=color)
        ax2.plot(amplitude_peaks, matched_filter_amplitude[amplitude_peaks], "x", color='tab:orange')
        ax2.tick_params(axis='y', labelcolor=color)

        fig.tight_layout()  # otherwise the right y-label is slightly clipped
        plt.show(block=True)'''

    def get_pmt_pulse_mf_shape(self):
        return self.pmt_pulse_mf_shape

    def set_pmt_pulse_mf_shape(self, new_shape: float):
        self.pmt_pulse_mf_shape = new_shape

    def get_pmt_pulse_mf_amp(self):
        return self.pmt_pulse_mf_amp

    def set_pmt_pulse_mf_amp(self, new_amp: float):
        self.pmt_pulse_mf_amp = new_amp

    def get_pmt_pulse_times(self):
        return self.pmt_pulse_times

    def set_pmt_pulse_times(self, new_pulse_times: list):
        self.pmt_pulse_times = new_pulse_times

    def calculate_charge(self):
        start = 0
        end = 0
        for i in range(self.get_pmt_object().get_setting("trigger_point"), int(self.get_pmt_pulse_peak_position())):
            if self.get_pmt_waveform()[i] <= self.get_pmt_waveform()[self.get_pmt_pulse_peak_position()]/self.get_pmt_object().get_setting("integration")[0]:
                start = i
                break
        for i in range(int(self.get_pmt_pulse_peak_position()), self.get_pmt_waveform_length()):
            if self.get_pmt_waveform()[i] >= self.get_pmt_waveform()[self.get_pmt_pulse_peak_position()]/self.get_pmt_object().get_setting("integration")[1]:
                end = i
                break
        self.set_pmt_pulse_start(start)
        self.set_pmt_pulse_end(end)
        self.set_pmt_pulse_charge(-1 * (np.sum(self.get_pmt_waveform()[
                                               start:end] - self.get_pmt_baseline())) / self.get_pmt_object().get_resistance())
