"""
    This is the base class for a Photomultiplier Tube
    To create an unique object:
        - PMT identifier ; PMT_id
        - Data identifier ; Data_id
    This will contain the results of waveform analysis from class PMT_Waveform
"""

import ROOT
import numpy as np


class PMT_Object:

    def __init__(self, pmt_id: str, data_id: str):
        self.pmt_id = pmt_id
        self.data_id = data_id

        self.number_of_events = 0
        self.sweep_bool = False
        self.template_bool = False

        # Default settings
        charge_cut = 6  # pC
        charge_range = [0, 100]
        nbins = 100
        amp_range = [0, 100]
        mf_shape_range = [-1, 1]
        mf_amp_range = [0, 100]
        sweep_range = [0, 500]
        pulse_time_threshold = 100
        apulse_region = 500
        resistance = 50 # Ohms
        mf_shape_threshold = 0.9
        mf_amp_threshold = 25
        baseline = 1000
        waveform_length = 8000
        trigger_point = 100
        integration = [0.3, 0.3]

        self.setting_dict = {
            "charge_cut"            : charge_cut,
            "charge_range"          : charge_range,
            "nbins"                 : nbins,
            "amp_range"             : amp_range,
            "mf_shape_range"        : mf_shape_range,
            "mf_amp_range"          : mf_amp_range,
            "sweep_range"           : sweep_range,
            "pulse_time_threshold"  : pulse_time_threshold,
            "apulse_region"         : apulse_region,
            "resistance"            : resistance,
            "mf_shape_threshold"    : mf_shape_threshold,
            "mf_amp_threshold"      : mf_amp_threshold,
            "waveform_length"       : waveform_length,
            "baseline"              : baseline,
            "trigger_point"         : trigger_point,
            "integration"           : integration
        }

        self.template_pmt_pulse = np.array([], dtype='float')

        self.histogram_names_list = ["pulse_charge_hist",
                                     "pulse_amplitude_hist",
                                     "apulse_charge_hist",
                                     "pulse_mf_shape_hist",
                                     "pulse_mf_amplitude_hist",
                                     "pulse_mf_shape_mf_amplitude_hist",
                                     "pulse_mf_shape_p_amplitude_hist",
                                     "pulse_mf_amplitude_p_amplitude",
                                     "pulse_peak_time_hist",
                                     "pulse_times_hist",
                                     "baseline_hist"
                                     ]

    def set_up_histograms(self):
        pmt_pulse_charge_hist = ROOT.TH1F(self.get_pmt_id() + "_pulse_charge_spectrum",
                                          self.get_pmt_id() + "_pulse_charge_spectrum",
                                          self.get_nbins(),
                                          self.get_charge_range()[0],
                                          self.get_charge_range()[1])

        pmt_pulse_amplitude_hist = ROOT.TH1F(self.get_pmt_id() + "_pulse_amplitude_spectrum",
                                             self.get_pmt_id() + "_pulse_amplitude_spectrum",
                                             self.get_nbins(),
                                             self.get_amp_range()[0],
                                             self.get_amp_range()[1])

        pmt_apulse_charge_hist = ROOT.TH1F(self.get_pmt_id() + "_apulse_charge_spectrum",
                                           self.get_pmt_id() + "_apulse_charge_spectrum",
                                           self.get_nbins(),
                                           self.get_charge_range()[0],
                                           self.get_charge_range()[1])

        pmt_pulse_mf_shape_hist = ROOT.TH1F(self.get_pmt_id() + "_mf_shape",
                                            self.get_pmt_id() + "_mf_shape",
                                            self.get_nbins(),
                                            self.get_mf_shape_range()[0],
                                            self.get_mf_shape_range()[1])

        pmt_pulse_mf_amplitude_hist = ROOT.TH1F(self.get_pmt_id() + "_mf_amplitude",
                                                self.get_pmt_id() + "_mf_amplitude",
                                                self.get_nbins(),
                                                self.get_mf_amp_range()[0],
                                                self.get_mf_amp_range()[1])

        pmt_pulse_mf_shape_mf_amplitude_hist = ROOT.TH2F(self.get_pmt_id() + "_mf_shape_mf_amplitude",
                                                         self.get_pmt_id() + "_mf_shape_mf_amplitude",
                                                         self.get_nbins(),
                                                         self.get_mf_shape_range()[0],
                                                         self.get_mf_shape_range()[1],
                                                         self.get_nbins(),
                                                         self.get_mf_amp_range()[0],
                                                         self.get_mf_amp_range()[1])

        pmt_pulse_mf_shape_p_amplitude_hist = ROOT.TH2F(self.get_pmt_id() + "_mf_shape_p_amplitude",
                                                        self.get_pmt_id() + "_mf_shape_p_amplitude",
                                                        self.get_nbins(),
                                                        self.get_mf_shape_range()[0],
                                                        self.get_mf_shape_range()[1],
                                                        self.get_nbins(),
                                                        self.get_amp_range()[0],
                                                        self.get_amp_range()[1])

        pmt_pulse_mf_amplitude_p_amplitude_hist = ROOT.TH2F(self.get_pmt_id() + "_mf_amplitude_p_amplitude",
                                                            self.get_pmt_id() + "_mf_amplitude_p_amplitude",
                                                            self.get_nbins(),
                                                            self.get_mf_amp_range()[0],
                                                            self.get_mf_amp_range()[1],
                                                            self.get_nbins(),
                                                            self.get_amp_range()[0],
                                                            self.get_amp_range()[1])

        pmt_pulse_peak_time_hist = ROOT.TH1I(self.get_pmt_id() + "_pulse_peak_times",
                                             self.get_pmt_id() + "_pulse_peak_times",
                                             self.get_waveform_length(),
                                             0,
                                             self.get_waveform_length())

        pmt_pulse_times_hist = ROOT.TH1I(self.get_pmt_id() + "_pulse_times",
                                         self.get_pmt_id() + "_pulse_times",
                                         self.get_waveform_length(),
                                         0,
                                         self.get_waveform_length())

        pmt_baseline_hist = ROOT.TH1F(self.get_pmt_id() + "_baseline",
                                      self.get_pmt_id() + "_baseline",
                                      self.get_nbins(),
                                      self.get_baseline_setting() - 10,
                                      self.get_baseline_setting() + 10)

        self.histogram_dict = {
            "pulse_charge_hist": pmt_pulse_charge_hist,
            "pulse_amplitude_hist": pmt_pulse_amplitude_hist,
            "apulse_charge_hist": pmt_apulse_charge_hist,
            "pulse_mf_shape_hist": pmt_pulse_mf_shape_hist,
            "pulse_mf_amplitude_hist": pmt_pulse_mf_amplitude_hist,
            "pulse_mf_shape_mf_amplitude_hist": pmt_pulse_mf_shape_mf_amplitude_hist,
            "pulse_mf_shape_p_amplitude_hist": pmt_pulse_mf_shape_p_amplitude_hist,
            "pulse_mf_amplitude_p_amplitude_hist": pmt_pulse_mf_amplitude_p_amplitude_hist,
            "pulse_peak_time_hist": pmt_pulse_peak_time_hist,
            "pulse_time_hists": pmt_pulse_times_hist,
            "baseline_hist": pmt_baseline_hist
        }

    def get_template_bool(self):
        return self.template_bool

    def set_template_bool(self,  new_bool: bool):
        self.template_bool = new_bool

    def get_normalisation_factor(self, vector: list):
        norm = 0.0
        for i in range(len(vector)):
            norm += vector[i] * vector[i]
        return np.sqrt(norm)

    def get_setting_dict(self):
        return self.setting_dict

    def get_setting(self, description: str):
        return self.get_setting_dict()[description]

    def set_setting(self, description: str, value):
        if type(value) == list:
            assert len(value) <= 2
        else:
            pass
        self.setting_dict[description] = value

    def get_resistance(self):
        return self.get_setting("resistance")

    def set_resistance(self, new_resistance: float):
        self.set_setting("resistance", new_resistance)

    def get_trigger_point(self):
        return self.get_setting("trigger_point")

    def set_trigger_point(self, new_trigger_point: int):
        self.set_setting("trigger_point", new_trigger_point)

    def get_waveform_length(self):
        return self.get_setting("waveform_length")

    def set_waveform_length(self, new_waveform_length: int):
        self.set_setting("waveform_length", new_waveform_length)

    def get_apulse_region(self):
        return self.get_setting("apulse_region")

    def set_apulse_region(self, new_apulse_pos: int):
        self.set_setting("apulse_region", new_apulse_pos)

    def get_histogram_dict(self):
        return self.histogram_dict

    def get_histogram(self, description: str):
        return self.get_histogram_dict()[description]

    def get_sweep_range(self):
        return self.get_setting("sweep_range")

    def get_sweep_range_min(self):
        return self.get_setting("sweep_range")[0]

    def get_sweep_range_max(self):
        return self.get_setting("sweep_range")[1]

    def set_sweep_range(self, new_sweep_range: list):
        assert len(new_sweep_range) == 2
        assert new_sweep_range[1] > new_sweep_range[0]
        self.set_setting("sweep_range", new_sweep_range)

    def get_pulse_time_threshold(self):
        return self.get_setting("pulse_time_threshold")

    def set_pulse_time_threshold(self, new_threshold: int):
        self.set_setting("pulse_time_threshold", new_threshold)

    def get_nbins(self):
        return self.get_setting("nbins")

    def set_nbins(self, new_nbins: int):
        self.set_setting("nbins", new_nbins)

    def get_baseline_setting(self):
        return self.get_setting("baseline")

    def set_baseline_setting(self, new_baseline: int):
        self.set_setting("baseline", new_baseline)

    def get_charge_range(self):
        return self.get_setting("charge_range")

    def set_charge_range(self, new_charge_range: list):
        assert len(new_charge_range) == 2
        assert new_charge_range[1] > new_charge_range[0]
        self.set_setting("charge_range", new_charge_range)

    def get_amp_range(self):
        return self.get_setting("amp_range")

    def set_amp_range(self, new_amp_range: list):
        assert len(new_amp_range) == 2
        assert new_amp_range[1] > new_amp_range[0]
        self.set_setting("amp_range", new_amp_range)

    def get_mf_shape_range(self):
        return self.get_setting("mf_shape_range")

    def set_mf_shape_range(self, new_mf_shape_range: list):
        assert len(new_mf_shape_range) == 2
        assert new_mf_shape_range[1] > new_mf_shape_range[0]
        self.set_setting("mf_shape_range", new_mf_shape_range)

    def get_mf_amp_range(self):
        return self.get_setting("mf_amp_range")

    def set_mf_amp_range(self, new_mf_amp_range: list):
        assert len(new_mf_amp_range) == 2
        assert new_mf_amp_range[1] > new_mf_amp_range[0]
        self.set_setting("mf_amp_range", new_mf_amp_range)

    def get_histogram_names(self):
        return self.histogram_names_list

    def get_histogram_names_dict(self):
        return self.histogram_dict

    def get_pmt_id(self):
        return self.pmt_id

    def set_pmt_id(self, new_pmt_id: str):
        self.pmt_id = new_pmt_id

    def get_data_id(self):
        return self.data_id

    def get_charge_cut(self):
        return self.get_setting("charge_cut")

    def get_event_number(self):
        return self.number_of_events

    def set_number_of_events(self, new_event_number):
        self.number_of_events = new_event_number

    def create_pmt_pulse_template(self, root_file_name: str, template_histogram_name: str):
        template_root_file = ROOT.TFile(root_file_name, "READ")
        template_histogram = template_root_file.Get(template_histogram_name)
        template_list = []
        for i_bin in range(int(template_histogram.GetEntries())):
            template_list.append(template_histogram.GetBinContent(i_bin))

        norm = self.get_normalisation_factor(template_list)

        self.set_template_pmt_pulse(np.array(template_list, dtype='float') / norm)
        # print(self.get_template_pmt_pulse())

        self.set_template_bool(True)

    def set_template_pmt_pulse(self, new_template_pmt_pulse: np.array):
        self.template_pmt_pulse = new_template_pmt_pulse

    def get_template_pmt_pulse(self):
        return self.template_pmt_pulse

    def get_pmt_pulse_charge_hist(self):
        return self.get_histogram("pulse_charge_hist")

    def fill_pmt_pulse_charge_hist(self, value: float):
        #print(value)
        self.get_histogram("pulse_charge_hist").Fill(value)

    def get_pmt_pulse_amplitude_hist(self):
        return self.get_histogram("pulse_amplitude_hist")

    def fill_pmt_pulse_amplitude_hist(self, value: float):
        self.get_histogram("pulse_amplitude_hist").Fill(value)

    def get_pmt_pulse_mf_shape_hist(self):
        return self.get_histogram("pulse_mf_shape_hist")

    def fill_pmt_pulse_mf_shape_hist(self, value: float):
        self.get_histogram("pulse_mf_shape_hist").Fill(value)

    def get_pmt_pulse_mf_amplitude_hist(self):
        return self.get_histogram("pulse_mf_amplitude_hist")

    def fill_pmt_pulse_mf_amplitude_hist(self, value: float):
        self.get_histogram("pulse_mf_amplitude_hist").Fill(value)

    def get_pmt_pulse_mf_shape_mf_amplitude(self):
        return self.get_histogram("pulse_mf_shape_mf_amplitude_hist")

    def fill_pmt_pulse_mf_shape_mf_amplitude(self, x_value: float, y_value: float):
        self.get_histogram("pulse_mf_shape_mf_amplitude_hist").Fill(x_value, y_value)

    def get_pmt_pulse_mf_shape_p_amplitude_hist(self):
        return self.get_histogram("pulse_mf_shape_p_amplitude_hist")

    def fill_pmt_pulse_mf_shape_p_amplitude_hist(self, x_value: float, y_value: float):
        self.get_histogram("pulse_mf_shape_p_amplitude_hist").Fill(x_value, y_value)

    def get_pmt_pulse_mf_amplitude_p_amplitude_hist(self):
        return self.get_histogram("pulse_mf_amplitude_p_amplitude_hist")

    def fill_pmt_pulse_mf_amplitude_p_amplitude_hist(self, x_value, y_value):
        self.get_histogram("pulse_mf_amplitude_p_amplitude_hist").Fill(x_value, y_value)

    def get_pmt_apulse_charge_hist(self):
        return self.get_histogram("apulse_charge_hist")

    def fill_pmt_apulse_charge_hist(self, value: float):
        self.get_histogram("apulse_charge_hist").Fill(value)

    def get_pmt_pulse_peak_time_hist(self):
        return self.get_histogram("pulse_peak_time_hist")

    def fill_pmt_pulse_peak_time_hist(self, value: int):
        self.get_histogram("pulse_peak_time_hist").Fill(value)

    def get_pmt_pulse_times_hist(self):
        return self.get_histogram("apulse_time_hists")

    def fill_pmt_pulse_times_hist(self, value: list):
        for i_value in range(len(value)):
            self.get_histogram("pulse_time_hists").Fill(value[i_value])

    def get_pmt_baseline_hist(self):
        return self.get_histogram("baseline_hist")

    def fill_pmt_baseline_hist(self, value: float):
        self.get_histogram("baseline_hist").Fill(value)

    def fill_pmt_hists(self, results: dict):

        #print(results)
        pulse_charge: float = results["pulse_charge"]
        pulse_amplitude: float = results["pulse_amplitude"]
        apulse_charge: float = results["apulse_charge"]
        mf_amplitude: float = results["pulse_mf_amp"]
        mf_shape: float = results["pulse_mf_shape"]
        pulse_peak_time: int = results["pulse_peak_time"]
        pulse_times: list = results["pulse_times"]
        baseline: float = results["baseline"]

        self.fill_pmt_pulse_charge_hist(pulse_charge)
        self.fill_pmt_pulse_amplitude_hist(pulse_amplitude)
        self.fill_pmt_apulse_charge_hist(apulse_charge)
        self.fill_pmt_pulse_mf_shape_hist(mf_shape)
        self.fill_pmt_pulse_mf_amplitude_hist(mf_amplitude)
        self.fill_pmt_pulse_mf_shape_mf_amplitude(mf_shape, mf_amplitude)
        self.fill_pmt_pulse_mf_shape_p_amplitude_hist(mf_shape, pulse_amplitude)
        self.fill_pmt_pulse_mf_amplitude_p_amplitude_hist(pulse_amplitude, mf_amplitude)
        self.fill_pmt_pulse_peak_time_hist(pulse_peak_time)
        self.fill_pmt_pulse_times_hist(pulse_times)
        self.fill_pmt_baseline_hist(baseline)

        self.set_number_of_events(self.get_event_number() + 1)

    def save_to_file(self, root_file_name: str):
        root_file = ROOT.TFile(root_file_name, "RECREATE")
        root_file.cd()
        for hist in self.get_histogram_dict().keys():
            self.get_histogram_dict()[hist].Write()
        root_file.Close()

    def save_histograms(self, directory: ROOT.TDirectory):
        directory.cd()
        for hist in self.get_histogram_dict().keys():
            self.get_histogram_dict()[hist].Write()

    def save_histogram(self, root_file: ROOT.TFile, hist, write_function: str):
        if write_function in ['RECREATE', "CREATE", "UPDATE"]:
            pass
        else:
            print("Invalid write function.")

    def set_sweep_bool(self, new_bool: bool):
        self.sweep_bool = new_bool

    def get_sweep_bool(self):
        return self.sweep_bool
