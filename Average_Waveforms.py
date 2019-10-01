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


def getAmplitude(ADC_values,calculated_baseline,peak_cell):
    amplitude = float(ADC_values[peak_cell]) - calculated_baseline
    return amplitude


def UpdateAveragePulse(average_pulse_trace, ADC_values,calculated_baseline,peak_cell):
    NORM = getNORM(ADC_values, calculated_baseline)
    # NORM = 1.0
    for i in range(len(average_pulse_trace)):
        # print i
        # print (ADC_values[int(peak_cell)-60+i])
        # print int(peak_cell)-60+i
        average_pulse_trace[i] += (float(ADC_values[i]) - calculated_baseline) / NORM
    return average_pulse_trace


def getTotalArea(ADC_values, pre_PULSE_region, calculated_baseline, waveform_length):
    integrated_area = 0
    for i in range(pre_PULSE_region, waveform_length):
        # print ADC_values[i], " ", baseline_calculated
        integrated_area += (float(ADC_values[i]) - calculated_baseline)
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


def Read_Data(PMT_data_filename,pre_PULSE_region,waveform_length,average_pulse_counter_vector,event_num_HT_vector,event_num_LTO_vector,average_pulse_trace_vector,Spectra_hist_vector,Spectra_cal_hist_vector,Rise_Time_hist_vector,Peak_Cell_hist_vector,Amplitude_hist_vector,Ratio_hist_vector,blank_num_vector):

    try:
        PMT_Data_file = open(PMT_data_filename, 'r')
    except FileNotFoundError as fnf_error:
        print(fnf_error)
        raise Exception("Error opening data file. Skip to the next file...")
        return

    newWaveform = False
    line_number = 0
    waveform_number = 0

    for index, line in enumerate(PMT_Data_file.readlines()[10:]):
        info = line.split(" ")

        if info[0] == "=" and info[1] == "HIT":
            newWaveform = True
            line_number = 0
        else:
            pass

        if newWaveform and line_number == 1:
            sl = int(info[1])
            ch = int(info[3])
            event_id_LTO = int(info[5])
            event_id_HT = int(info[7])
            peak_cell = int(info[27])
            charge = float(info[29])
            rise_time = float(info[39])
            if int(event_id_HT) != 0:
                event_num_HT_vector[sl][ch] += 1
            elif int(event_id_LTO) != 0:
                event_num_LTO_vector[sl][ch] += 1
            else:
                blank_num_vector[sl][ch] += 1 # Maybe add in a blank waveform counter

        elif newWaveform and line_number == 2:
            ADC_values = []
            for i_ADC in range(waveform_length):
                ADC_values.append(float(info[i_ADC]))

            # Store some waveforms
            if average_pulse_counter_vector[sl][ch] < 5:
                trace = TH1F("Waveform_" + str(sl) + "_" + str(ch) + "_" + str(waveform_number),
                                 "Waveform_" + str(sl) + "_" + str(ch) + "_" + str(waveform_number),
                                 waveform_length, 0, waveform_length)
                trace.GetXaxis().SetTitle("Sample Number")
                trace.GetYaxis().SetTitle("ADC Counts")
                for i in range(waveform_length):
                    trace.SetBinContent(i, ADC_values[i])
                trace.Write("", TFile.kOverwrite)
                del trace

            Calculated_Baseline = getBaseline(ADC_values,pre_PULSE_region)

            Spectra_hist_vector[sl][ch].Fill(abs(charge))

            cal_charge = getTotalArea(ADC_values, pre_PULSE_region, Calculated_Baseline,waveform_length)
            cal_amplitude = getAmplitude(ADC_values,Calculated_Baseline,peak_cell)

            if cal_charge > 0:
                Ratio = abs(cal_amplitude / cal_charge)
            else:
                Ratio = 10000000000000

            Spectra_cal_hist_vector[sl][ch].Fill(abs(cal_charge))
            Peak_Cell_hist_vector[sl][ch].Fill(peak_cell)
            Rise_Time_hist_vector[sl][ch].Fill(rise_time)
            Amplitude_hist_vector[sl][ch].Fill(abs(cal_amplitude))
            Ratio_hist_vector[sl][ch].Fill(Ratio)

            if average_pulse_counter_vector[sl][ch] < 100:
                average_pulse_trace_vector[sl][ch] = UpdateAveragePulse(average_pulse_trace_vector[sl][ch], ADC_values,Calculated_Baseline,peak_cell)
                average_pulse_counter_vector[sl][ch] += 1

        line_number += 1

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

        # First create all the conatiners for the histograms

        map_HT = TH2I("Mapping_HT", "Mapping_HT", topology[0], 0, topology[0], topology[1], 0, topology[1])
        map_LTO = TH2I("Mapping_LTO", "Mapping_LTO", topology[0], 0, topology[0], topology[1], 0, topology[1])
        map = TH2I("Mapping", "Mapping", topology[0], 0, topology[0], topology[1], 0, topology[1])

        Spectra_hist_vector = []
        Spectra_cal_hist_vector = []
        Peak_Cell_hist_vector = []
        Rise_Time_hist_vector = []
        Amplitude_hist_vector = []
        Ratio_hist_vector = []
        Average_Waveform_hist_vector = []

        # Create a blank waveform for each OM
        average_pulse_trace_vector = []
        average_pulse_counter_vector = []
        event_num_HT_vector = []
        event_num_LTO_vector = []
        blank_num_vector = []

        # Fill the containers with the objects for each OM in this Loop

        for slot_num in range(topology[0]):

            Spectra_hist_vector.append([])
            Spectra_cal_hist_vector.append([])
            Peak_Cell_hist_vector.append([])
            Rise_Time_hist_vector.append([])
            Amplitude_hist_vector.append([])
            Ratio_hist_vector.append([])
            Average_Waveform_hist_vector.append([])
            average_pulse_trace_vector.append([])
            average_pulse_counter_vector.append([])
            event_num_HT_vector.append([])
            event_num_LTO_vector.append([])
            blank_num_vector.append([])

            for channel_num in range(topology[1]):

                Spectra_hist = TH1F("Spectra_" + str(slot_num) + "_" + str(channel_num),"Spectra_" + str(slot_num) + "_" + str(channel_num), 200, 0, 200000)
                Spectra_hist.GetXaxis().SetTitle("Raw Charge")
                Spectra_hist_vector[slot_num].append(Spectra_hist)

                Spectra_cal_hist = TH1F("Spectra_cal_" + str(slot_num) + "_" + str(channel_num),"Spectra_cal_" + str(slot_num) + "_" + str(channel_num), 200, 0, 200000)
                Spectra_cal_hist.GetXaxis().SetTitle("Raw Charge")
                Spectra_cal_hist_vector[slot_num].append(Spectra_cal_hist)

                Peak_Cell_hist = TH1F("Peak_Cell_"+str(slot_num)+"_"+str(channel_num),"Peak_Cell_"+str(slot_num)+"_"+str(channel_num), waveform_length,0,waveform_length)
                Peak_Cell_hist.GetXaxis().SetTitle("Peak Cell Position /Sample Number")
                Peak_Cell_hist_vector[slot_num].append(Peak_Cell_hist)

                Rise_Time_hist = TH1F("Rise_Time_"+str(slot_num)+"_"+str(channel_num),"Rise_Time_"+str(slot_num)+"_"+str(channel_num),200,50,200)
                Rise_Time_hist.GetXaxis().SetTitle("Rise Time")
                Rise_Time_hist_vector[slot_num].append(Rise_Time_hist)

                Amplitude_hist = TH1F("Amplitudes_"+str(slot_num)+"_"+str(channel_num),"Amplitudes_"+str(slot_num)+"_"+str(channel_num),200,0,2100)
                Amplitude_hist.GetXaxis().SetTitle("Amplitude /ADC Counts")
                Amplitude_hist_vector[slot_num].append(Amplitude_hist)

                Ratio_hist = TH1F("Ratio_"+str(slot_num)+"_"+str(channel_num),"Ratio_"+str(slot_num)+"_"+str(channel_num),200,0,0.05)
                Ratio_hist.GetXaxis().SetTitle("Amplitude/Raw Charge")
                Ratio_hist_vector[slot_num].append(Ratio_hist)

                average_pulse_trace_vector[slot_num].append([0]*waveform_length)
                average_pulse_counter_vector[slot_num].append(0)
                event_num_HT_vector[slot_num].append(0)
                event_num_LTO_vector[slot_num].append(0)
                blank_num_vector[slot_num].append(0)

                average_pulse_trace_hist = TH1F("Waveform_" + str(slot_num) + "_" + str(channel_num) + "_average","Waveform_" + str(slot_num) + "_" + str(channel_num) + "_average",waveform_length, 0, waveform_length)
                average_pulse_trace_hist.GetXaxis().SetTitle("Sample Number")
                average_pulse_trace_hist.GetYaxis().SetTitle("Averaged ADC Counts /AU")
                Average_Waveform_hist_vector[slot_num].append(average_pulse_trace_hist)

                del Spectra_hist
                del Spectra_cal_hist
                del Peak_Cell_hist
                del Rise_Time_hist
                del Amplitude_hist
                del Ratio_hist
                del average_pulse_trace_hist

        # Now run ove the files
        print(">>> Reading Data Files... ")
        processing_start = time.time()
        temp_start = processing_start
        for file_num in range(len(PMT_data_filenames)):
            print (">>> File: ", PATH + PMT_data_filenames[file_num].rstrip())
            Read_Data(PATH + PMT_data_filenames[file_num].rstrip(),
                      pre_PULSE_region,
                      waveform_length,
                      average_pulse_counter_vector,
                      event_num_HT_vector,
                      event_num_LTO_vector,
                      average_pulse_trace_vector,
                      Spectra_hist_vector,
                      Spectra_cal_hist_vector,
                      Rise_Time_hist_vector,
                      Peak_Cell_hist_vector,
                      Amplitude_hist_vector,
                      Ratio_hist_vector,
                      blank_num_vector)
            intermidiate = time.time()
            time_length = intermidiate - processing_start
            print(">>>\n>>>  %.3f seconds.\n" % (intermidiate-temp_start))
            temp_start = intermidiate
            print("Processed %.2f of data..." % (float(file_num+1)/float(len(PMT_data_filenames))*100))
            Estimate = (time_length/float(file_num+1)) * (len(PMT_data_filenames) - file_num - 1)
            print(">>> Estimated time till termination %.3f seconds\n\n" % Estimate)



        for slot_num in range(topology[0]):
            root_file.mkdir("Slot"+str(slot_num))
            root_file.GetDirectory("Slot" + str(slot_num)).cd()
            Channel_Check = topology[1]
            for channel_num in range(topology[1]):
                root_file.GetDirectory("Slot" + str(slot_num)).mkdir("Channel" + str(channel_num))
                root_file.GetDirectory("Slot" + str(slot_num)).GetDirectory("Channel" + str(channel_num)).cd()

                Spectra_hist_vector[slot_num][channel_num].Write("", TFile.kOverwrite)
                Spectra_cal_hist_vector[slot_num][channel_num].Write("", TFile.kOverwrite)
                Rise_Time_hist_vector[slot_num][channel_num].Write("", TFile.kOverwrite)
                Amplitude_hist_vector[slot_num][channel_num].Write("", TFile.kOverwrite)
                Peak_Cell_hist_vector[slot_num][channel_num].Write("", TFile.kOverwrite)
                Ratio_hist_vector[slot_num][channel_num].Write("", TFile.kOverwrite)

                for bin in range(len(average_pulse_trace_vector[slot_num][channel_num])):
                    # Avoid dividing by zero
                    if average_pulse_counter_vector[slot_num][channel_num] == 0:
                        average_pulse_counter_vector[slot_num][channel_num] = 1

                    Average_Waveform_hist_vector[slot_num][channel_num].SetBinContent(bin,
                                                           (average_pulse_trace_vector[slot_num][channel_num][bin] / average_pulse_counter_vector[slot_num][channel_num]) * 1000.0)
                    Average_Waveform_hist_vector[slot_num][channel_num].SetBinError(bin, 1.0)

                Average_Waveform_hist_vector[slot_num][channel_num].Write("", TFile.kOverwrite)
                map.Fill(slot_num, channel_num, event_num_HT_vector[slot_num][channel_num] + event_num_LTO_vector[slot_num][channel_num])
                map_HT.Fill(slot_num, channel_num, event_num_HT_vector[slot_num][channel_num])
                map_LTO.Fill(slot_num, channel_num, event_num_LTO_vector[slot_num][channel_num])

                if event_num_HT_vector[slot_num][channel_num] + event_num_LTO_vector[slot_num][channel_num] + blank_num_vector[slot_num][channel_num] == 0:
                    root_file.GetDirectory("Slot"+str(slot_num)).Delete("Channel"+str(channel_num)+";1")
                    #print("\n>>> Deleted Slot: ",slot_num, " Channel: ", channel_num, "  as there were no events\n")
                    Channel_Check = Channel_Check - 1
            if Channel_Check == 0:
                root_file.Delete("Slot" + str(slot_num) + ";1")
                print('\n>>> Deleted Slot: ', slot_num, '  as there were no events\n')

        root_file.cd()
        map.Write("", TFile.kOverwrite)
        map_HT.Write("", TFile.kOverwrite)
        map_LTO.Write("", TFile.kOverwrite)
        root_file.Close()

    print(">>>")


def createTemplateVector(Directory,template_num, waveform_length, pre_PULSE_region):
    Directory.cd()
    selected_average_waveform_trace = Directory.Get("Waveform_" + str(template_num[0]) + "_" + str(template_num[1]) + "_average")

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


def checkWaveforms(root_filename, topology, waveform_length, pre_PULSE_region, template_num):
    print(">>> Checking the shapes of the waveforms... ")

    try:
        test_file = open(root_filename, "r")
    except FileNotFoundError as fnf_error:
        print(fnf_error)
        raise Exception("Error opening .root file. TERMINATING.")
        return

    root_file = TFile(root_filename, "update")

    shape_index_2D_hist = TH2F("shape_index_2D_hist1", "shape_index_2D_hist1", topology[0], 0, topology[0],
                                topology[1], 0, topology[1])

    Template_Channel =  root_file.GetDirectory("Slot"+str(template_num[0])).GetDirectory("Channel"+str(template_num[1]))

    template_vector, template_normalisation_factor = createTemplateVector(Template_Channel, template_num,waveform_length,
                                                                            pre_PULSE_region)

    Template_trace_hist = TH1F("Template_trace_hist", "Template_trace_hist", waveform_length, 0, waveform_length)

    for i in range(waveform_length):
        Template_trace_hist.SetBinContent(i, float(template_vector[i]))

    Template_trace_hist.Write("", TFile.kOverwrite)

    del Template_trace_hist

    for slot_num in range(topology[0]):
        Slot_Check = True

        Slot = root_file.GetDirectory("Slot"+str(slot_num))
        if Slot == None:
            Slot_Check = False

        if Slot_Check:
            Slot.cd()
            for channel_num in range(topology[1]):
                Channel_Check = True

                Channel = Slot.GetDirectory("Channel"+str(channel_num))
                if Channel == None:
                    Channel_Check = False

                if Channel_Check:
                    Channel.cd()
                    average_waveform_trace = Channel.Get("Waveform_" + str(slot_num) + "_" + str(channel_num) + "_average")
                    waveform_vector = []
                    for bin in range(len(template_vector)):
                        # pulse_baseline = getBaseline()
                        waveform_vector.append(float(average_waveform_trace.GetBinContent(bin)))
                    calculated_baseline = getBaseline(waveform_vector, pre_PULSE_region)
                    # calculated_baseline = 0.0
                    waveform_normalisation_factor = getNORM(waveform_vector, calculated_baseline)

                    if waveform_normalisation_factor == 0:
                        waveform_normalisation_factor = 1

                    for k in range(len(waveform_vector)):
                        waveform_vector[k] = (waveform_vector[k] - calculated_baseline) / waveform_normalisation_factor

                    shape_index = 0

                    for i in range(len(template_vector)):

                        shape_index += (float(waveform_vector[i]) * float(template_vector[i]))


                    shape_index_2D_hist.Fill(slot_num, channel_num, shape_index * 100)

                    del average_waveform_trace

    root_file.cd()
    shape_index_2D_hist.Write("", TFile.kOverwrite)
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


def main(run,DATA_PATH,ROOT_PATH,topology,trigger_point,pulse_length,template_num):
    DATA_PATH = DATA_PATH+"run_"+run+"/"
    #ROOT_PATH = "/Users/willquinn/Documents/PhD/SuperNEMO/GV_XW_ComData/"
    Data_files = getDataFiles(DATA_PATH +"filenames_"+run+".txt")

    root_filename = "Average_Waveforms_run"+run+".root"

    AnalyseWaveforms(topology, DATA_PATH, Data_files, ROOT_PATH + root_filename, trigger_point, pulse_length)

    checkWaveforms(ROOT_PATH + root_filename, topology, pulse_length, trigger_point,template_num)

    # CheckWaveformsAgain(DATA_PATH,ROOT_PATH+root_filename,Data_files,topology,pulse_length,trigger_point)


if __name__ == '__main__':
    print (">>> Input the run number: ")
    run = str(raw_input(">>> "))
    print(">>>")
    '''
    print (">>> Input the Data PATH: ")
    DATA_PATH = str(raw_input(">>> "))
    print(">>>")
    print(">>> Input the output file PATH: ")
    OUTPUT_PATH = str(raw_input(">>> "))
    print(">>> Analysing...")
    '''

    #run = "214"
    DATA_PATH = "/unix/nemo3/SN_Calo_Commissioning_Runs/GVETO_XWALL_Runs/"
    OUTPUT_PATH = "/home/wquinn/GV_XW_ComData/"

    topology = [16,16]
    trigger_point = 160
    pulse_length= 1024
    template_num = [14,15]

    main(run,DATA_PATH,OUTPUT_PATH,topology,trigger_point,pulse_length,template_num)

    end = time.time()
    print(">>>\n>>> The program lasted %.3f seconds.\n" % (end - start))
