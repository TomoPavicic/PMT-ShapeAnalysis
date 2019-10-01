###########################################################################

## William Quinn - SuperNEMO Calorimeter Data analysis ##
## Waveform,Spectra and mapping code ##
## Version 4.0 20/3/2019 ##
## X-Wall Gamma Veto ##

###########################################################################

import time

start = time.time()

from ROOT import TFile, TH1F, TH2F, TMultiGraph, TGraph, TCanvas, TH2I
from math import log10
from numpy import sqrt


def getNORM(ADC_values, calculated_baseline):
    intermediate_normalisation_factor = 0
    for i in range(len(ADC_values) - 1):
        intermediate_normalisation_factor += (float(ADC_values[i]) - calculated_baseline) ** 2
    normalisation_factor = abs(intermediate_normalisation_factor) ** (0.5)
    return normalisation_factor


def getAmplitude(ADC_values,pre_PULSE_region,peak_cell):
    calculated_baseline = getBaseline(ADC_values,pre_PULSE_region)
    amplitude = float(ADC_values[peak_cell]) - calculated_baseline
    return amplitude


def UpdateAveragePulse(average_pulse_trace, ADC_values, pre_PULSE_region, peak_cell):
    baseline_calculated = getBaseline(ADC_values, pre_PULSE_region)
    NORM = getNORM(ADC_values, baseline_calculated)
    # NORM = 1.0
    for i in range(len(average_pulse_trace)):
        # print i
        # print (ADC_values[int(peak_cell)-60+i])
        # print int(peak_cell)-60+i
        average_pulse_trace[i] += (float(ADC_values[i]) - baseline_calculated) / NORM
    return average_pulse_trace


def getTotalArea(ADC_values, pre_PULSE_region, integration_window):
    baseline_calculated = getBaseline(ADC_values, pre_PULSE_region)
    integrated_area = 0
    for i in range(pre_PULSE_region, pre_PULSE_region + integration_window):
        # print ADC_values[i], " ", baseline_calculated
        integrated_area += (float(ADC_values[i]) - baseline_calculated)
    Total_Area = abs(integrated_area)
    return Total_Area


def getBaseline(ADC_values, pre_PULSE_region):
    calculated_baseline = 0
    pre_TRIGGER_region = pre_PULSE_region - 40
    for i in range(pre_TRIGGER_region):
        calculated_baseline += float(ADC_values[i])
    averaged_calculated_baseline = calculated_baseline / pre_TRIGGER_region
    return averaged_calculated_baseline


def getBaselineAmplitude(ADC_values, pre_PULSE_region):
    baseline_calculated = getBaseline(ADC_values, pre_PULSE_region)
    neg_amp = baseline_calculated
    for i in range(pre_PULSE_region):
        # print neg_amp," ",ADC_values[i]
        if float(ADC_values[i]) < neg_amp:
            neg_amp = float(ADC_values[i])
    return neg_amp


def getReadFile(PMT_data_filename, slot_num, channel_num, pre_PULSE_region, average_pulse_trace,
                       average_pulse_counter, waveform_length,Spectra,Spectra_cal,Peak_Cell_Hist,Rise_Time_Hist,Amplitude_Hist,Ratio_Hist):
    PMT_Data_file = open(PMT_data_filename, 'r')
    i_waveform = 0
    i_line = 0
    newWaveform = False
    event_num_HT = 0
    event_num_LTO = 0
    event = False

    for index, line in enumerate(PMT_Data_file.readlines()[10:]):
        info = line.split(" ")
        if info[0] == "=" and info[1] == "HIT":
            newWaveform = True
            i_waveform += 1
            i_line = 0

        if newWaveform and i_line == 1:
            sl = int(info[1])
            ch = int(info[3])
            event_id_HT = int(info[7])
            event_id_LTO = int(info[5])
            peak_cell = int(info[27])
            charge = float(info[29])
            rise_time = float(info[39])

            if int(event_id_HT) != 0 and int(slot_num) == int(sl) and int(channel_num) == int(ch):
                event_num_HT += 1
                event = True
            elif int(event_id_LTO) != 0 and int(slot_num) == int(sl) and int(channel_num) == int(ch):
                event_num_LTO += 1
                event = True
            else:
                event = False

        elif event and i_line == 2:
            ADC_values = []
            for i_ADC in range(waveform_length):
                ADC_values.append(float(info[i_ADC]))

            # Store some waveforms
            if average_pulse_counter < 5:
                trace = TH1F("Waveform_" + str(slot_num) + "_" + str(channel_num) + "_" + str(i_waveform),
                             "Waveform_" + str(slot_num) + "_" + str(channel_num) + "_" + str(i_waveform),
                             waveform_length, 0, waveform_length)
                trace.GetXaxis().SetTitle("Sample Number")
                trace.GetYaxis().SetTitle("ADC Counts")
                for i in range(waveform_length):
                    trace.SetBinContent(i, ADC_values[i])
                trace.Write("", TFile.kOverwrite)
                del trace

            # print ("Averaging a waveform")
            # integrated_area = getTotalArea(ADC_values,trigger_point,400)
            # if integrated_area > 1000:
            average_pulse_trace = UpdateAveragePulse(average_pulse_trace, ADC_values, pre_PULSE_region, peak_cell)
            average_pulse_counter += 1

            i_line = -1
            event = False

            Spectra.Fill(abs(charge))

            cal_charge = getTotalArea(ADC_values,pre_PULSE_region,864)
            cal_amplitude = getAmplitude(ADC_values,pre_PULSE_region,peak_cell)
            if cal_charge > 0:
                Ratio = abs(cal_amplitude/cal_charge)
            else:
                Ratio = 10000000000000

            Spectra_cal.Fill(abs(cal_charge))
            Peak_Cell_Hist.Fill(peak_cell)
            Rise_Time_Hist.Fill(rise_time)
            Amplitude_Hist.Fill(abs(cal_amplitude))
            Ratio_Hist.Fill(Ratio)

        i_line += 1

        if average_pulse_counter == 1000:
            return average_pulse_trace, average_pulse_counter, event_num_HT, event_num_LTO

    # print ("Slot ",slot_num, " Channel ", channel_num, "number  of waveforms")
    return average_pulse_trace, average_pulse_counter, event_num_HT, event_num_LTO


def getfile_again(shape_vector, PMT_data_filename, pre_PULSE_region, template_vector, waveform_length):
    PMT_Data_file = open(PMT_data_filename, 'r')

    for index, line in enumerate(PMT_Data_file.readlines()[10:]):
        if (index + 1) % 3 == 1:
            info = line.split(" ")
            sl = info[1]
            ch = info[3]
            event_id = info[7]
            peak_cell = info[27]

            if int(ch) == 13:
                event_id = 0
        if (index + 1) % 3 == 2:
            if int(event_id) != 0:
                # print event_id
                if 60 <= int(peak_cell) <= 400:
                    waveform_vector = [] * waveform_length
                    ADC_values = line.split(" ")
                    baseline_calculated = getBaseline(ADC_values, pre_PULSE_region)
                    NORM = 0
                    for i in range(waveform_length):
                        # print i
                        # print ADC_values[int(peak_cell)-60+i]
                        # print int(peak_cell)-60+i
                        waveform_vector.append(float(ADC_values[int(peak_cell) - 60 + i]) - baseline_calculated)
                        NORM += (float(ADC_values[int(peak_cell) - 60 + i]) - baseline_calculated) ** 2
                    NORM = NORM ** 0.5
                    shape_index = 0.0
                    for i in range(waveform_length):
                        waveform_vector[i] = waveform_vector[i] / NORM
                        shape_index += waveform_vector[i] * template_vector[i]
                    # print int(sl) + 20*int(ch)
                    shape_vector[int(sl) + int(ch) * 20].append(shape_index)

    return shape_vector


def getDataFiles(DataFiles_filename):
    DataFiles = open(DataFiles_filename, 'r')
    Data_files_list = []
    for line in DataFiles.readlines():
        Data_files_list.append(str(line))
    return Data_files_list


def AnalyseWaveforms(topology, PATH, PMT_data_filenames, root_filename, pre_PULSE_region, waveform_length):
    try:  # Will not create a new file if one exists already
        waveformfile = open(root_filename, "r")
        print(">>> File exists")
    except IOError:
        print(">>> File not found")
        print(">>> Creating .root file with the Average Waveforms")

        root_file = TFile(root_filename, 'update')

        map_HT = TH2I("Mapping_HT", "Mapping_HT", topology[0], 0, topology[0], topology[1], 0, topology[1])
        map_LTO = TH2I("Mapping_LTO", "Mapping_LTO", topology[0], 0, topology[0], topology[1], 0, topology[1])
        map = TH2I("Mapping", "Mapping", topology[0], 0, topology[0], topology[1], 0, topology[1])

        for slot_num in range(10, topology[0]):

            for channel_num in range(topology[1]):

                Spectra = TH1F("Spectra_" + str(slot_num) + "_" + str(channel_num),"Spectra_" + str(slot_num) + "_" + str(channel_num), 200, 0, 100000)
                Spectra.GetXaxis().SetTitle("Raw Charge")

                Spectra_cal = TH1F("Spectra_cal_" + str(slot_num) + "_" + str(channel_num),"Spectra_cal_" + str(slot_num) + "_" + str(channel_num), 200, 0, 100000)
                Spectra_cal.GetXaxis().SetTitle("Raw Charge")

                Peak_Cell = TH1F("Peak_Cell_"+str(slot_num)+"_"+str(channel_num),"Peak_Cell_"+str(slot_num)+"_"+str(channel_num), 200,200,350)
                Peak_Cell.GetXaxis().SetTitle("Peak Cell Position /Sample Number")

                Rise_Time = TH1F("Rise_Time_"+str(slot_num)+"_"+str(channel_num),"Rise_Time_"+str(slot_num)+"_"+str(channel_num),200,50,200)
                Rise_Time.GetXaxis().SetTitle("Rise Time")

                Amplitude = TH1F("Amplitudes_"+str(slot_num)+"_"+str(channel_num),"Amplitudes_"+str(slot_num)+"_"+str(channel_num),200,0,4000)
                Amplitude.GetXaxis().SetTitle("Amplitude /ADC Counts")

                Ratio = TH1F("Ratio_"+str(slot_num)+"_"+str(channel_num),"Ratio_"+str(slot_num)+"_"+str(channel_num),200,0,0.05)
                Ratio.GetXaxis().SetTitle("Amplitude/Raw Charge")

                event_num_HT = 0
                event_num_LTO = 0
                # Create a blank waveform for each OM
                average_pulse_trace = [0] * waveform_length
                average_pulse_counter = 0
                # Create a histogram to be stored
                average_pulse_trace_hist = TH1F("Waveform_" + str(slot_num) + "_" + str(channel_num) + "_average",
                                                "Waveform_" + str(slot_num) + "_" + str(channel_num) + "_average",
                                                waveform_length, 0, waveform_length)
                average_pulse_trace_hist.GetXaxis().SetTitle("Sample Number")
                average_pulse_trace_hist.GetYaxis().SetTitle("Averaged ADC Counts /AU")

                for file_num in range(len(PMT_data_filenames)):
                    average_pulse_trace, average_pulse_counter, event_number_HT, event_number_LTO = getReadFile(
                        PATH + PMT_data_filenames[file_num].rstrip(),
                        slot_num,
                        channel_num,
                        pre_PULSE_region,
                        average_pulse_trace,
                        average_pulse_counter,
                        waveform_length,
                        Spectra,
                        Spectra_cal,
                        Peak_Cell,
                        Rise_Time,
                        Amplitude,
                        Ratio)

                    event_num_HT += event_number_HT
                    event_num_LTO += event_number_LTO

                #print("Slot: ", slot_num, " Channel: ", channel_num, " Number of events: ",event_num_HT + event_num_LTO)
                map.Fill(slot_num, channel_num, event_num_HT + event_num_LTO)
                map_HT.Fill(slot_num, channel_num, event_num_HT)
                map_LTO.Fill(slot_num, channel_num, event_num_LTO)
                #print("Slot: ", slot_num, " Channel: ", channel_num, " Number of waveforms averaged: ",average_pulse_counter)

                Spectra.Write("", TFile.kOverwrite)
                Spectra_cal.Write("", TFile.kOverwrite)
                Peak_Cell.Write("", TFile.kOverwrite)
                Rise_Time.Write("", TFile.kOverwrite)
                Amplitude.Write("", TFile.kOverwrite)
                Ratio.Write("", TFile.kOverwrite)
                del Spectra
                del Spectra_cal
                del Peak_Cell
                del Rise_Time
                del Amplitude
                del Ratio

                for bin in range(len(average_pulse_trace)):
                    # Avoid dividing by zero
                    if average_pulse_counter == 0:
                        average_pulse_counter = 1
                    average_pulse_trace_hist.SetBinContent(bin,(average_pulse_trace[bin] / average_pulse_counter) * 1000.0)
                    average_pulse_trace_hist.SetBinError(bin, 1.0)

                average_pulse_trace_hist.Write()

                del average_pulse_trace

            #print(">>> Progress: ", float(slot_num + 1) / float(topology[0] - 11) * 100, "%")
            print ("Slot: ",slot_num," complete")
            print (">>>")
            print(">>>")

        map.Write("", TFile.kOverwrite)
        map_HT.Write("", TFile.kOverwrite)
        map_LTO.Write("", TFile.kOverwrite)
        root_file.Close()

    print(">>>")


def createTemplateVector(root_file, slot_num, channel_num, waveform_length, pre_PULSE_region):
    selected_average_waveform_trace = root_file.Get("Waveform_" + str(slot_num) + "_" + str(channel_num) + "_average")

    template_vector = []
    for bin in range(waveform_length):
        template_vector.append(float(selected_average_waveform_trace.GetBinContent(bin)))
    calculated_baseline = getBaseline(template_vector, pre_PULSE_region)
    # calculated_baseline = 0.0
    normalisation_factor = getNORM(template_vector, calculated_baseline)
    # print "Template normalisation factor = ", normalisation_factor
    temp = 0

    for i in range(waveform_length):
        # print template_vector[i]
        template_vector[i] = (template_vector[i] - calculated_baseline) / normalisation_factor
        temp += template_vector[i] * template_vector[i]
    # print "temp ",temp

    # normalisation_factor = 1.0
    del selected_average_waveform_trace

    return template_vector, normalisation_factor


def checkWaveforms(root_filename, topology, waveform_length, pre_PULSE_region):
    print(">>> Checking the shapes of the waveforms... ")

    try:
        test_file = open(root_filename, "r")
    except IOError:
        print(">>> No file in PATH:", root_filename)
        return

    root_file = TFile(root_filename, "update")

    shape_index_2D_hist1 = TH2F("shape_index_2D_hist1", "shape_index_2D_hist1", topology[0], 0, topology[0],
                                topology[1], 0, topology[1])
    # shape_index_2D_hist2 = TH2F("shape_index_2D_hist2","shape_index_2D_hist2",topology[0],0,topology[0],topology[1],0,topology[1])

    template_vector1, template_normalisation_factor1 = createTemplateVector(root_file, 11, 5, waveform_length,
                                                                            pre_PULSE_region)
    # template_vector2,template_normalisation_factor2 = createTemplateVector(root_file,0,2,waveform_length,pre_PULSE_region)

    # print "templates"
    # print template_vector1
    # print template_vector2

    Template_trace_hist1 = TH1F("Template_trace_hist1", "Template_trace_hist1", waveform_length, 0, waveform_length)
    # Template_trace_hist2 = TH1F("Template_trace_hist2","Template_trace_hist2",waveform_length,0,waveform_length)
    for i in range(waveform_length):
        # print "saving templates"
        # print template_vector1[i]
        Template_trace_hist1.SetBinContent(i, float(template_vector1[i]))
        # Template_trace_hist2.SetBinContent(i,float(template_vector2[i]))
    Template_trace_hist1.Write("", TFile.kOverwrite)
    # Template_trace_hist2.Write("",TFile.kOverwrite)
    del Template_trace_hist1
    # del Template_trace_hist2

    for slot_num in range(10, topology[0]):

        for channel_num in range(topology[1]):

            # Sample_Graph = TGraph()
            # Comparison_Graph = TMultiGraph("Comparison_Graph_"+str(slot_num)+"_"+str(channel_num),"Comparison_Graph_"+str(slot_num)+"_"+str(channel_num))

            average_waveform_trace = root_file.Get("Waveform_" + str(slot_num) + "_" + str(channel_num) + "_average")
            waveform_vector = []
            for bin in range(len(template_vector1)):
                # pulse_baseline = getBaseline()
                waveform_vector.append(float(average_waveform_trace.GetBinContent(bin)))
            calculated_baseline = getBaseline(waveform_vector, pre_PULSE_region)
            # calculated_baseline = 0.0
            waveform_normalisation_factor = getNORM(waveform_vector, calculated_baseline)

            if waveform_normalisation_factor == 0:
                #print(">>> blank ")
                # print waveform_vector
                waveform_normalisation_factor = 1
            # print waveform_normalisation_factor**0.5

            temp2 = 0
            for k in range(len(waveform_vector)):
                waveform_vector[k] = (waveform_vector[k] - calculated_baseline) / waveform_normalisation_factor
                # temp2 += waveform_vector[k]*waveform_vector[k]
            # print "temp2 ",temp2

            shape_index1 = 0
            # shape_index2 = 0
            for i in range(len(template_vector1)):
                '''
                shape_index1 += (float(waveform_vector[i])*float(template_vector1[i]))/(waveform_normalisation_factor*template_normalisation_factor1)
                shape_index2 += (float(waveform_vector[i])*float(template_vector2[i]))/(waveform_normalisation_factor*template_normalisation_factor2)
                '''

                shape_index1 += (float(waveform_vector[i]) * float(template_vector1[i]))
                # shape_index2 += (float(waveform_vector[i])*float(template_vector2[i]))

            '''
            if shape_index1 > 0:
                print("Slot: ", slot_num, " Channel: ", channel_num, " = ", shape_index1 * 100)
                # print ("Slot: ",slot_num," Channel: ",channel_num," = ",shape_index2*100)
            else:
                print("Slot: ", slot_num, " Channel: ", channel_num, " = ", shape_index1)
                # print ("Slot: ",slot_num," Channel: ",channel_num," = ",shape_index2)
            '''

            shape_index_2D_hist1.Fill(slot_num, channel_num, shape_index1 * 100)
            # shape_index_2D_hist2.Fill(slot_num,channel_num,shape_index2*100)

            del average_waveform_trace
    shape_index_2D_hist1.Write("", TFile.kOverwrite)
    # shape_index_2D_hist2.Write("",TFile.kOverwrite)
    root_file.Close()


def CheckWaveformsAgain(PATH, root_filename, PMT_data_filenames, topology, waveform_length, pre_PULSE_region):
    try:
        test_file = open(root_filename, "r")
    except IOError:
        print(">>> No file in PATH:", root_filename)

    root_file = TFile(root_filename, "update")
    shape_index_2D_hist3 = TH2F("shape_index_2D_hist3", "shape_index_2D_hist3", topology[0], 0, topology[0],
                                topology[1], 0, topology[1])

    shape_vector = [[] for i in range(260)]
    # print shape_vector

    template_vector, template_normalisation_factor = createTemplateVector(root_file, 0, 1, waveform_length,
                                                                          pre_PULSE_region)

    for i in range(waveform_length):
        template_vector[i] = template_vector[i] / template_normalisation_factor

    for file_num in range(len(PMT_data_filenames)):
        shape_vector = getfile_again(shape_vector, PATH + PMT_data_filenames[file_num].rstrip(), pre_PULSE_region,
                                     template_vector, waveform_length)
    average_shapes = []

    slot_num = 0
    channel_num = 0
    for i in range(260):
        average = 0

        for j in range(len(shape_vector[i])):
            average += shape_vector[i][j]
        average = average / len(shape_vector[i])
        average_shapes.append(average)
        shape_index_2D_hist3.Fill(slot_num, channel_num, average)
        slot_num += 1
        if slot_num == 20:
            slot_num = 0
            channel_num += 1

    shape_index_2D_hist3.Write()

    print(shape_vector)


def main(run,DATA_PATH,ROOT_PATH):
    DATA_PATH = DATA_PATH+"run_"+run+"/"
    #ROOT_PATH = "/Users/willquinn/Documents/PhD/SuperNEMO/GV_XW_ComData/"
    Data_files = getDataFiles(DATA_PATH +"filenames_"+run+".txt")

    root_filename = "Average_Waveforms_run"+run+".root"

    # Define the topology of the Wall
    topology = [15, 16]
    trigger_point = 160
    pulse_length = 1024

    AnalyseWaveforms(topology, DATA_PATH, Data_files, ROOT_PATH + root_filename, trigger_point, pulse_length)

    checkWaveforms(ROOT_PATH + root_filename, topology, pulse_length, trigger_point)

    # CheckWaveformsAgain(DATA_PATH,ROOT_PATH+root_filename,Data_files,topology,pulse_length,trigger_point)


if __name__ == '__main__':
    print ("Input the run number: ")
    run = str(raw_input())
    print(">>>")
    print ("Input the Data PATH: ")
    DATA_PATH = str(raw_input())
    print(">>>")
    print("Input the output file PATH: ")
    OUTPUT_PATH = str(raw_input())
    print(">>> Analysing...")
    main(run,DATA_PATH,OUTPUT_PATH)
    end = time.time()
    print(">>>\n>>> The program lasted %.3f seconds.\n" % (end - start))
