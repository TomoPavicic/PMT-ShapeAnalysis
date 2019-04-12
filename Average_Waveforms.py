###########################################################################

            ## William Quinn - SuperNEMO Calorimeter Data analysis ##
                  ## Waveform,Spectra and mapping code ##
                        ## Version 4.0 20/3/2019 ##
                            ## Italian Shape ##

###########################################################################

import time
start = time.time()
 
from ROOT import TFile,TH1F,TH2F,TMultiGraph,TGraph
from math import log10
from numpy import sqrt

def getNORM(ADC_values,calculated_baseline):
    intermediate_normalisation_factor = 0
    for i in range(len(ADC_values)-1):
        intermediate_normalisation_factor += (float(ADC_values[i]) - calculated_baseline)**2
    normalisation_factor = abs(intermediate_normalisation_factor)**(0.5)
    return normalisation_factor

def UpdateAveragePulse(average_pulse_trace,ADC_values,pre_PULSE_region,peak_cell):
    baseline_calculated = getBaseline(ADC_values,pre_PULSE_region)
    #NORM = getNORM(ADC_values,baseline_calculated)
    NORM = 1.0
    for i in range(len(average_pulse_trace)):
        #print i
        #print ADC_values[int(peak_cell)-60+i]
        #print int(peak_cell)-60+i
        average_pulse_trace[i] += (float(ADC_values[int(peak_cell)-60+i]) - baseline_calculated)/NORM
    return average_pulse_trace

def getTotalArea(ADC_values,pre_PULSE_region,integration_window):
    baseline_calculated = getBaseline(ADC_values,pre_PULSE_region)
    integrated_area = 0
    for i in range(pre_PULSE_region,pre_PULSE_region+integration_window):
        #print ADC_values[i], " ", baseline_calculated
        integrated_area += (float(ADC_values[i]) - baseline_calculated)
    Total_Area = abs(integrated_area)
    return Total_Area

def getBaseline(ADC_values,pre_PULSE_region):
    calculated_baseline = 0
    pre_TRIGGER_region = pre_PULSE_region - 40
    for i in range(pre_TRIGGER_region):
        calculated_baseline += float(ADC_values[i])
    averaged_calculated_baseline = calculated_baseline/pre_TRIGGER_region
    return averaged_calculated_baseline

def getBaselineAmplitude(ADC_values,pre_PULSE_region):
    baseline_calculated = getBaseline(ADC_values,pre_PULSE_region)
    neg_amp = baseline_calculated
    for i in range(pre_PULSE_region):
        #print neg_amp," ",ADC_values[i]
        if float(ADC_values[i]) < neg_amp:
            neg_amp = float(ADC_values[i])
    return neg_amp

def getfile(PMT_data_filename,root_filename,file_num,slot_num,channel_num,pre_PULSE_region,average_pulse_trace,average_pulse_counter):

    PMT_Data_file = open(PMT_data_filename,'r')
    
    for index,line in enumerate(PMT_Data_file.readlines()[10:]):
        if (index + 1) % 3 ==1:
            info = line.split(" ")
            sl = info[1]
            ch = info[3]
            event_id = info[7]
            peak_cell = info[27]
        if (index + 1) % 3 ==2:
            if int(event_id) != 0:
                #print event_id
                if int(slot_num) == int(sl):
                    if int(channel_num) == int(ch):
                        if 60 <= int(peak_cell) <= 400:
                            ADC_values = line.split(" ")
                            #integrated_area = getTotalArea(ADC_values,trigger_point,400)
                            #if integrated_area > 1000:
                            average_pulse_trace = UpdateAveragePulse(average_pulse_trace,ADC_values,pre_PULSE_region,peak_cell)
                            average_pulse_counter += 1

    return average_pulse_trace,average_pulse_counter

def getfile_again(shape_vector,PMT_data_filename,pre_PULSE_region,template_vector,waveform_length):

    PMT_Data_file = open(PMT_data_filename,'r')
    
    for index,line in enumerate(PMT_Data_file.readlines()[10:]):
        if (index + 1) % 3 ==1:
            info = line.split(" ")
            sl = info[1]
            ch = info[3]
            event_id = info[7]
            peak_cell = info[27]

            if int(ch) == 13:
                event_id = 0
        if (index + 1) % 3 ==2:
            if int(event_id) != 0:
                #print event_id
                if 60 <= int(peak_cell) <= 400:
                    waveform_vector = []*waveform_length
                    ADC_values = line.split(" ")
                    baseline_calculated = getBaseline(ADC_values,pre_PULSE_region)
                    NORM = 0
                    for i in range(waveform_length):
                        #print i
                        #print ADC_values[int(peak_cell)-60+i]
                        #print int(peak_cell)-60+i
                        waveform_vector.append(float(ADC_values[int(peak_cell)-60+i]) - baseline_calculated)
                        NORM += (float(ADC_values[int(peak_cell)-60+i]) - baseline_calculated)**2
                    NORM = NORM**0.5
                    shape_index = 0.0
                    for i in range(waveform_length):
                        waveform_vector[i] = waveform_vector[i]/NORM
                        shape_index += waveform_vector[i]*template_vector[i]
                    #print int(sl) + 20*int(ch)
                    shape_vector[int(sl) + int(ch)*20].append(shape_index)


    return shape_vector
                            

    return average_pulse_trace,average_pulse_counter

def getDataFiles(DataFiles_filename):
    DataFiles = open(DataFiles_filename,'r')
    Data_files_list = []
    for line in DataFiles.readlines():
        Data_files_list.append(str(line))
    return Data_files_list

def AnalyseWaveforms(topology,PATH,PMT_data_filenames,root_filename,pre_PULSE_region,waveform_length):
    try: 
        waveformfile = open(root_filename,"r")
    except IOError:
        print ">>> File not found"
        print ">>> Creating .root file with the Average Waveforms"

        root_file = TFile(root_filename,'update')
        
        for slot_num in range(topology[0]):

            for channel_num in range(topology[1]):

                average_pulse_trace = [0]*waveform_length
                average_pulse_counter = 0
                average_pulse_trace_hist = TH1F("Waveform_"+str(slot_num)+"_"+str(channel_num)+"_average","Waveform_"+str(slot_num)+"_"+str(channel_num)+"_average",waveform_length,0,waveform_length)

                for file_num in range(len(PMT_data_filenames)):

                    average_pulse_trace,average_pulse_counter = getfile(PATH+PMT_data_filenames[file_num].rstrip(),root_filename,file_num,slot_num,channel_num,pre_PULSE_region,average_pulse_trace,average_pulse_counter)

                for bin in range(len(average_pulse_trace)):
                    if average_pulse_counter == 0:
                        average_pulse_counter = 1
                    average_pulse_trace_hist.SetBinContent(bin, (average_pulse_trace[bin]/average_pulse_counter)*1000000.0)
                    average_pulse_trace_hist.SetBinError(bin,1.0)

                average_pulse_trace_hist.Write()

                del average_pulse_trace

            print ">>> Progress: ",float(slot_num + 1)/float(topology[0]) * 100,"%"

        root_file.Close()

    print ">>>"

def createTemplateVector(root_file,slot_num,channel_num,waveform_length,pre_PULSE_region):
    selected_average_waveform_trace = root_file.Get("Waveform_"+str(slot_num)+"_"+str(channel_num)+"_average")

    template_vector = []
    for bin in range(waveform_length):
        template_vector.append(float(selected_average_waveform_trace.GetBinContent(bin)))
    #calculated_baseline = getBaseline(template_vector,pre_PULSE_region)
    calculated_baseline = 0.0
    normalisation_factor = getNORM(template_vector,calculated_baseline)
    #print "Template normalisation factor = ", normalisation_factor
    temp = 0

    for i in range(waveform_length):
        #print template_vector[i]
        template_vector[i] = template_vector[i]/normalisation_factor
        temp += template_vector[i]*template_vector[i]
    #print "temp ",temp

    normalisation_factor = 1.0

    return template_vector,normalisation_factor


def checkWaveforms(root_filename,topology,waveform_length,pre_PULSE_region):
    print ">>> Checking the shapes of the waveforms... "

    try: 
        test_file = open(root_filename,"r")
    except IOError:
        print ">>> No file in PATH:",root_filename

    root_file = TFile(root_filename,"update")

    shape_index_2D_hist1 = TH2F("shape_index_2D_hist1","shape_index_2D_hist1",topology[0],0,topology[0],topology[1],0,topology[1])
    shape_index_2D_hist2 = TH2F("shape_index_2D_hist2","shape_index_2D_hist2",topology[0],0,topology[0],topology[1],0,topology[1])

    template_vector1,template_normalisation_factor1 = createTemplateVector(root_file,0,1,waveform_length,pre_PULSE_region)
    template_vector2,template_normalisation_factor2 = createTemplateVector(root_file,0,2,waveform_length,pre_PULSE_region)

    #print "templates"
    #print template_vector1
    #print template_vector2

    Template_trace_hist1 = TH1F("Template_trace_hist1","Template_trace_hist1",waveform_length,0,waveform_length)
    Template_trace_hist2 = TH1F("Template_trace_hist2","Template_trace_hist2",waveform_length,0,waveform_length)
    for i in range(waveform_length):
        #print "saving templates"
        #print template_vector1[i]
        Template_trace_hist1.SetBinContent(i,float(template_vector1[i]))
        Template_trace_hist2.SetBinContent(i,float(template_vector2[i]))
    Template_trace_hist1.Write("",TFile.kOverwrite)
    Template_trace_hist2.Write("",TFile.kOverwrite)
    del Template_trace_hist1
    del Template_trace_hist2

    for slot_num in range(topology[0]):
         for channel_num in range(topology[1]):
            #Sample_Graph = TGraph()
            #Comparison_Graph = TMultiGraph("Comparison_Graph_"+str(slot_num)+"_"+str(channel_num),"Comparison_Graph_"+str(slot_num)+"_"+str(channel_num))
            average_waveform_trace = root_file.Get("Waveform_"+str(slot_num)+"_"+str(channel_num)+"_average")
            waveform_vector = []
            for bin in range(len(template_vector1)):
                #pulse_baseline = getBaseline()
                waveform_vector.append(float(average_waveform_trace.GetBinContent(bin)))
            #calculated_baseline = getBaseline(waveform_vector,pre_PULSE_region)
            calculated_baseline = 0.0
            waveform_normalisation_factor = getNORM(waveform_vector,calculated_baseline)

            temp = 0
            for j in range(len(waveform_vector)):
                temp += waveform_vector[j]*waveform_vector[j]
            #print "temp ",sqrt(temp)
            waveform_normalisation_factor = sqrt(temp)

            if waveform_normalisation_factor == 0:
                print ">>> blank "
                #print waveform_vector
                waveform_normalisation_factor = 1
            #print waveform_normalisation_factor**0.5

            temp2 = 0

            for k in range(len(waveform_vector)):
                waveform_vector[k] = waveform_vector[k]/waveform_normalisation_factor
                temp2 += waveform_vector[k]*waveform_vector[k]
            #print "temp2 ",temp2

            shape_index1 = 0
            shape_index2 = 0
            i = 0
            for i in range(len(template_vector1)):

                '''
                shape_index1 += (float(waveform_vector[i])*float(template_vector1[i]))/(waveform_normalisation_factor*template_normalisation_factor1)
                shape_index2 += (float(waveform_vector[i])*float(template_vector2[i]))/(waveform_normalisation_factor*template_normalisation_factor2)
                '''

                shape_index1 += (float(waveform_vector[i])*float(template_vector1[i]))
                shape_index2 += (float(waveform_vector[i])*float(template_vector2[i]))

            if shape_index1 > 0 and shape_index2 > 0:
                print "Slot: ",slot_num," Channel: ",channel_num," = ",shape_index1*100
                print "Slot: ",slot_num," Channel: ",channel_num," = ",shape_index2*100
            else:
                print "Slot: ",slot_num," Channel: ",channel_num," = ",shape_index1
                print "Slot: ",slot_num," Channel: ",channel_num," = ",shape_index2

            shape_index_2D_hist1.Fill(slot_num,channel_num,shape_index1*100)
            shape_index_2D_hist2.Fill(slot_num,channel_num,shape_index2*100)
            
            del average_waveform_trace
    shape_index_2D_hist1.Write("",TFile.kOverwrite)
    shape_index_2D_hist2.Write("",TFile.kOverwrite)
    root_file.Close()

def CheckWaveformsAgain(PATH,root_filename,PMT_data_filenames,topology,waveform_length,pre_PULSE_region):
    try: 
        test_file = open(root_filename,"r")
    except IOError:
        print ">>> No file in PATH:",root_filename

    root_file = TFile(root_filename,"update")
    shape_index_2D_hist3 = TH2F("shape_index_2D_hist3","shape_index_2D_hist3",topology[0],0,topology[0],topology[1],0,topology[1])

    shape_vector = [[] for i in range(260)]
    #print shape_vector

    template_vector,template_normalisation_factor = createTemplateVector(root_file,0,1,waveform_length,pre_PULSE_region)

    for i in range(waveform_length):
        template_vector[i] = template_vector[i]/template_normalisation_factor

    for file_num in range(len(PMT_data_filenames)):
        shape_vector = getfile_again(shape_vector,PATH+PMT_data_filenames[file_num].rstrip(),pre_PULSE_region,template_vector,waveform_length)
    average_shapes = []

    slot_num = 0
    channel_num = 0
    for i in range(260):
        average = 0
        
        for j in range(len(shape_vector[i])):
            average += shape_vector[i][j]
        average = average/len(shape_vector[i])
        average_shapes.append(average)
        shape_index_2D_hist3.Fill(slot_num,channel_num,average)
        slot_num += 1
        if slot_num == 20:
            slot_num = 0
            channel_num += 1

    shape_index_2D_hist3.Write()

    print shape_vector



def main():
    DATA_PATH = "/Users/willquinn/Documents/PhD/SNEMO_ComData_Analysis/French_CaloWall/data/"
    ROOT_PATH = "/Users/willquinn/Documents/PhD/SNEMO_ComData_Analysis/French_CaloWall/ROOT_files/"
    Data_files = getDataFiles(DATA_PATH+"filenames_104.txt")

    root_filename = "Average_Waveforms_run104_test.root"

    #Define the topology of the Wall
    topology = [20,13]
    trigger_point = 160
    pulse_length = 600

    AnalyseWaveforms(topology,DATA_PATH,Data_files,ROOT_PATH+root_filename,trigger_point,pulse_length)
    
    #checkWaveforms(ROOT_PATH+root_filename,topology,pulse_length,trigger_point)

    CheckWaveformsAgain(DATA_PATH,ROOT_PATH+root_filename,Data_files,topology,pulse_length,trigger_point)

 
if __name__ == '__main__':
    main()
    end = time.time()
    print ">>>\n>>> The program lasted %.3f seconds.\n" % (end-start)
 