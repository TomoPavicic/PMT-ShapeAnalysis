import ROOT
from scr.PMT_Object import PMT_Object


class PMT_Array:

    def __init__(self, topology: list, data_id: str):
        self.pmt_topology = topology
        self.pmt_object_array = []
        self.pmt_total_number = topology[0] * topology[1]

        for i_row in range(topology[0]):
            for i_col in range(topology[1]):
                pmt_number = i_col + i_row*topology[1]
                pmt_object = PMT_Object(str(pmt_number), data_id)
                self.append_pmt_object_array(pmt_object)
                del pmt_object

    def save_to_file(self, output_root_filename: str):
        output_root_file = ROOT.TFile(output_root_filename, "RECREATE")

        for i_pmt in range(self.get_pmt_total_number()):
            if self.get_pmt_object_number(i_pmt).get_event_number() > 0:
                output_root_file.cd()
                output_root_file.mkdir(self.get_pmt_object_number(i_pmt).get_pmt_id())
                directory = output_root_file.GetDirectory(self.get_pmt_object_number(i_pmt).get_pmt_id())
                self.get_pmt_object_number(i_pmt).save_histograms(directory)

    def get_pmt_topology(self):
        return self.pmt_topology

    def get_pmt_object_array(self):
        return self.pmt_object_array

    def get_pmt_total_number(self):
        return self.pmt_total_number

    def append_pmt_object_array(self, pmt_object: PMT_Object):
        self.pmt_object_array.append(pmt_object)

    def get_pmt_object_position(self, pmt_position: list):
        assert len(pmt_position) < 3
        if len(pmt_position) == 1:
            pmt_number = pmt_position[0]
        else:
            pmt_number = pmt_position[0] * self.get_pmt_topology()[1] + pmt_position[1]
        return self.pmt_object_array[pmt_number]

    def get_pmt_object_number(self, pmt_number: int):
        return self.get_pmt_object_array()[pmt_number]

    def set_pmt_templates(self, template_root_file_name: str, template_histogram_name_list: list):
        for i_pmt in range(self.get_pmt_total_number()):
            self.get_pmt_object_number(i_pmt).create_pmt_pulse_template(template_root_file_name, template_histogram_name_list[i_pmt])

    def apply_setting(self, config_file_name: str):
        if config_file_name is not None:
            stuff_list = self.read_config_file(config_file_name)

            for i_pmt in range(self.get_pmt_total_number()):
                for i_setting in range(len(stuff_list)):
                    # print(self.get_pmt_object_number(i_pmt).get_setting_dict())
                    self.get_pmt_object_number(i_pmt).set_setting(stuff_list[i_setting][0], stuff_list[i_setting][1])
                    # print(self.get_pmt_object_number(i_pmt).get_setting_dict())
                self.get_pmt_object_number(i_pmt).set_up_histograms()
        else:
            for i_pmt in range(self.get_pmt_total_number()):
                self.get_pmt_object_number(i_pmt).set_up_histograms()

    def set_pmt_id(self, pmt_id: str, pmt_object_number: int):
        self.get_pmt_object_number(pmt_object_number).set_pmt_id(pmt_id)

    def set_sweep_bool(self, new_bool: bool):
        for i_pmt in range(self.get_pmt_total_number()):
            self.get_pmt_object_number(i_pmt).set_sweep_bool(new_bool)

    def read_config_file(self, config_file_name: str):
        try:
            config_file = open(config_file_name, 'r')
        except FileNotFoundError as fnf_error:
            print(fnf_error)
            raise Exception("Error opening config file")

        output_list = []
        for i_line, line in enumerate(config_file.readlines()):
            tokens = line.split(" ")
            if tokens[0] == '#':
                description = tokens[1].strip()
                value = tokens[3].split(",")
                value_ = []
                if len(value) == 1:
                    try:
                        value_ = int(value[0].strip())
                    except:
                        value_ = float(value[0].strip())
                else:
                    for i in range(len(value)):
                        try:
                            value_.append(int(value[i].strip()))
                        except:
                            value_.append(float(value[i].strip()))

                temp_list = [description, value_]
                output_list.append(temp_list)

        return output_list

    def fit_bismuth_function(self):
        # TODO: This function needs OOP

        for i_pmt in range(self.get_pmt_total_number()):
            if "GAO607" in self.get_pmt_object_number(i_pmt).get_pmt_id():
                if self.get_pmt_object_number(i_pmt).get_pmt_pulse_charge_hist().GetEntries() == 0:
                    pass
                else:
                    canvas_name = self.get_pmt_object_number(i_pmt).get_pmt_id()
                    canvas = ROOT.TCanvas(canvas_name, canvas_name)
                    fit = ROOT.TF1("fit",
                                   "[0]*"
                                   "(7.08*TMath::Gaus(x,[1],[2]) "
                                   " + 1.84*TMath::Gaus(x,[1]*(1 + 72.144/975.651),[2]*1.036) "
                                   " + 0.44*TMath::Gaus(x,[1]*(1 + 84.154/975.651),[2]*1.042)) "
                                   " + 0.464*(exp(0.254*x)/(1 + exp((x - 28.43)/2.14)))",
                                   33, 41)

                    fit.SetParNames("A", "mu", "sigma")

                    fit.SetParLimits(0, 0, 400)
                    fit.SetParLimits(1, 34, 37)
                    fit.SetParLimits(2, 0.8, 2)
                    fit.SetParameters(319, 36, 1.09)

                    self.get_pmt_object_number(i_pmt).get_pmt_pulse_charge_hist().Fit("fit", "", "", 33, 41)
                    self.get_pmt_object_number(i_pmt).get_pmt_pulse_charge_hist().SetXTitle("Charge /pC")
                    self.get_pmt_object_number(i_pmt).get_pmt_pulse_charge_hist().SetYTitle("Counts")
                    self.get_pmt_object_number(i_pmt).get_pmt_pulse_charge_hist().SetTitle("Bi Integrated Charge Spectrum")

                    resolution = (fit.GetParameter(2) / fit.GetParameter(1)) * 100.0
                    chi2 = fit.GetChisquare() / fit.GetNDF()
                    mu = fit.GetParameter(1)
                    mu_err = fit.GetParError(1)
                    sigma = fit.GetParameter(2)
                    sigma_err = fit.GetParError(2)

                    canvas.SetGrid()
                    canvas.Update()

                    canvas.Draw()
                    ROOT.gStyle.SetOptStat(11111111)
                    ROOT.gStyle.SetOptFit(1)
                    ROOT.gStyle.SetStatY(0.9)
                    ROOT.gStyle.SetStatX(0.9)
                    ROOT.gStyle.SetStatW(0.8)
                    ROOT.gStyle.SetStatH(0.1)
                    canvas.SaveAs(canvas_name+".pdf", "pdf")

            elif "GAO612" in self.get_pmt_object_number(i_pmt).get_pmt_id():
                if self.get_pmt_object_number(i_pmt).get_pmt_pulse_charge_hist().GetEntries() == 0:
                    pass
                else:
                    canvas_name = self.get_pmt_object_number(i_pmt).get_pmt_id()
                    canvas = ROOT.TCanvas(canvas_name, canvas_name)

                    fit = ROOT.TF1("fit",
                                   "[0]*"
                                   "(7.08*TMath::Gaus(x,[1],[2]) "
                                   " + 1.84*TMath::Gaus(x,[1]*(1 + 72.144/975.651),[2]*1.036) "
                                   " + 0.44*TMath::Gaus(x,[1]*(1 + 84.154/975.651),[2]*1.042)) "
                                   " + 0.515*(exp(0.2199*x)/(1 + exp((x - 31.68)/2.48)))",
                                   37, 44)

                    fit.SetParNames("A", "mu", "sigma")

                    fit.SetParLimits(0, 0, 400)
                    fit.SetParLimits(1, 37, 40)
                    fit.SetParLimits(2, 0.8, 2)

                    fit.SetParameters(319, 39, 1.09)

                    self.get_pmt_object_number(i_pmt).get_pmt_pulse_charge_hist().Fit("fit", "", "", 37, 44)
                    self.get_pmt_object_number(i_pmt).get_pmt_pulse_charge_hist().SetXTitle("Charge /pC")
                    self.get_pmt_object_number(i_pmt).get_pmt_pulse_charge_hist().SetYTitle("Counts")
                    self.get_pmt_object_number(i_pmt).get_pmt_pulse_charge_hist().SetTitle("Bi Integrated Charge Spectrum")
                    self.get_pmt_object_number(i_pmt).get_pmt_pulse_charge_hist().Draw()

                    resolution = (fit.GetParameter(2) / fit.GetParameter(1)) * 100.0
                    chi2 = fit.GetChisquare() / fit.GetNDF()
                    mu = fit.GetParameter(1)
                    mu_err = fit.GetParError(1)
                    sigma = fit.GetParameter(2)
                    sigma_err = fit.GetParError(2)

                    canvas.SetGrid()
                    canvas.Update()

                    canvas.Draw()
                    ROOT.gStyle.SetOptStat(11111111)
                    ROOT.gStyle.SetOptFit(1)
                    ROOT.gStyle.SetStatY(0.9)
                    ROOT.gStyle.SetStatX(0.9)
                    ROOT.gStyle.SetStatW(0.8)
                    ROOT.gStyle.SetStatH(0.1)
                    canvas.SaveAs(canvas_name+".pdf", "pdf")
