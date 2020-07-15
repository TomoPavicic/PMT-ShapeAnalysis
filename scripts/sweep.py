import sys
sys.path.insert(1, '..')

from scr.PMT_Array import PMT_Array
from functions.other_functions import sncd_parse_arguments, get_run_number, get_data_path
from functions.data_reader_functions import process_crd_file, read_filenames
import time as TIME
import ROOT


def main():
    args = sncd_parse_arguments()
    data_file_names = args.i
    config_file_name = args.c
    sweep_bool = args.sweep
    template_file_names = args.t
    output_file_name = args.o

    run_number = get_run_number(data_file_names)

    topology = [20, 20]

    template_list = read_filenames(template_file_names)
    list_of_data_file_names = read_filenames(data_file_names)

    pmt_array = PMT_Array(topology, run_number)
    if config_file_name is not None:
        pmt_array.apply_setting(config_file_name)
    if template_file_names is not None:
        pmt_array.set_pmt_templates(template_list, "template")
    if sweep_bool == "True":
        pmt_array.set_sweep_bool(True)

    print(">>> Looping of files...")

    data_path = get_data_path(data_file_names)

    waveform_output_root_file = ROOT.TFile("waveforms"+run_number+".root", "RECREATE")

    temp_start = TIME.time()
    for index, data_file_name in enumerate(list_of_data_file_names):
        process_crd_file(data_path + data_file_name, pmt_array, waveform_output_root_file)

        intermediate = TIME.time()
        time_length = intermediate - temp_start
        print(">>>\n>>>  %.3f s.\n" % (intermediate - temp_start))
        temp_start = intermediate
        print("Processed {} of the {} data files...".format(index + 1, len(list_of_data_file_names)))
        estimate = time_length * (len(list_of_data_file_names) - index - 1)
        if estimate >= 60 and estimate < 3600:
            minute = estimate/60.0
            print(">>> Estimated time till termination %.1f minutes\n\n" % minute)
        elif estimate >= 3600:
            hours = estimate / 3600.0
            print(">>> Estimated time till termination %.1f hours\n\n" % hours)
        else:
            print(">>> Estimated time till termination %.1f seconds\n\n" % estimate)

    hist_2D_shapes = ROOT.TH2F("PMT_Pulse_Shapes",
                               "PMT_Pulse_Shapes",
                               pmt_array.get_pmt_topology()[0],
                               0,
                               pmt_array.get_pmt_topology()[0],
                               pmt_array.get_pmt_topology()[1],
                               0,
                               pmt_array.get_pmt_topology()[1])

    hist_2D_shapes.SetXTitle("Column OM")
    hist_2D_shapes.SetYTitle("Row OM")

    hist_2D_shapes_SD = ROOT.TH2F("PMT_Pulse_Shapes_SD",
                                  "PMT_Pulse_Shapes_SD",
                                  pmt_array.get_pmt_topology()[0],
                                  0,
                                  pmt_array.get_pmt_topology()[0],
                                  pmt_array.get_pmt_topology()[1],
                                  0,
                                  pmt_array.get_pmt_topology()[1])

    hist_2D_mapping = ROOT.TH2F("PMT_Event_Mapping",
                                "PMT_Event_Mapping",
                                pmt_array.get_pmt_topology()[0],
                                0,
                                pmt_array.get_pmt_topology()[0],
                                pmt_array.get_pmt_topology()[1],
                                0,
                                pmt_array.get_pmt_topology()[1])

    hist_2D_num_waveforms = ROOT.TH2F("PMT_Pre_Pulse_Rate",
                                      "PMT_Pre_Pulse_Rate",
                                      pmt_array.get_pmt_topology()[0],
                                      0,
                                      pmt_array.get_pmt_topology()[0],
                                      pmt_array.get_pmt_topology()[1],
                                      0,
                                      pmt_array.get_pmt_topology()[1])

    hist_2D_shapes_SD.SetXTitle("Column OM")
    hist_2D_shapes_SD.SetYTitle("Row OM")

    for i_row in range(pmt_array.get_pmt_topology()[0]):
        for i_col in range(pmt_array.get_pmt_topology()[1]):
            print("Row: ", i_row, "Col: ", i_col, "Event Number: ", pmt_array.get_pmt_object_position([i_row, i_col]).get_event_number())
            hist_2D_shapes.Fill(i_col, i_row, pmt_array.get_pmt_object_position([i_row, i_col]).get_pmt_pulse_mf_shape_hist().GetMean())
            hist_2D_shapes_SD.Fill(i_col, i_row, pmt_array.get_pmt_object_position([i_row, i_col]).get_pmt_pulse_mf_shape_hist().GetStdDev())
            hist_2D_mapping.Fill(i_col, i_row, pmt_array.get_pmt_object_position([i_row, i_col]).get_event_number())
            #hist_2D_num_waveforms.Fill

    hist_2D_root_file = ROOT.TFile("2D_shape_hist_"+run_number+"_.root", "RECREATE")
    hist_2D_root_file.cd()
    hist_2D_shapes.Write()
    hist_2D_shapes_SD.Write()
    hist_2D_mapping.Write()
    hist_2D_root_file.Close()

    pmt_array.save_to_file(output_file_name)


def create_template_histogram_temp():
    source_file = ROOT.TFile(
        "/Users/willquinn/Documents/PhD/SuperNEMO/SNEMO_ComData_Analysis/GV_XW_ComData/Output_files/Average_Waveforms_run214.root",
        "READ")
    source_waveform = ROOT.TH1F(source_file.Get("Waveform_10_3_0"))
    # source_file.Close()

    output_file = ROOT.TFile("/Users/willquinn/Documents/PhD/SuperNEMO/SNEMO_ComData_Analysis/template.root",
                             "RECREATE")
    template_hist = ROOT.TH1F("template", "template", 80, 0, 80)
    output_file.cd()

    pmt_object = PMT_Object("template", "template")

    temp_list = []
    for i_bin in range(int(source_waveform.GetEntries())):
        temp_list.append(str(source_waveform.GetBinContent(i_bin)))
        # print(i_bin, " ", str(source_waveform.GetBinContent(i_bin + 1)))
    pmt_waveform = PMT_Waveform(temp_list, pmt_object)

    for i_bin in range(80):
        template_hist.SetBinContent(i_bin + 1,
                                    pmt_waveform.get_pmt_waveform()[260 + i_bin] - pmt_waveform.get_pmt_baseline())
    output_file.Write()
    output_file.Close()


if __name__ == '__main__':
    main()
