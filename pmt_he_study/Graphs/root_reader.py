import ROOT
import os
import numpy as np
import matplotlib.pyplot as plt
import sys

'''
files = ["201025_A1400_B1400_t1827_output.root", "201026_A1400_B1400_t1010_output.root",  "201030_A1400_B1400_t1010_output.root"]

file_date = ["201025","201026","201030"]
file_time = ["1827","1010","1010"]

def extracting_data_example(date, time):
  #creating the file path
  file_path = "ROOT_files/1400V/"
  file_name = date + "_A1400_B1400_t" + time + "_output.root"

  #Opening the root file
  root_file  = ROOT.TFile(file_path + file_name, "READ")
  root_file.cd()
  
  #Getting the histogram from the tree
  apulse_hist = root_file.Get("201025_GAO607_apulse_num_1400V")
  
  #Diplaying relavent information
  print("Number of entries: " + str(apulse_hist.GetEntries()))
  print("Number of bins: " + str(apulse_hist.GetNbinsX()))
  print("Contents of a single bin: " + str(apulse_hist.GetBinContent(2)))
  
  #in histograms the 1st bin (index 0) is the overflow bin
  #in afterpulse the 2nd bin (index 1) is the 0-1 (so no after pulses)
  all_bins = []
  for i in range(2, apulse_hist.GetNbinsX()):
    all_bins.append(apulse_hist.GetBinContent(i))

  print("All relavent bins (at least 1 afterpulse): ")
  print(all_bins)

extracting_data_example(file_date[0],file_time[0])
'''

#ROOT_files_path = "ROOT_files/1400V/"
def directory_list(folder_path):
  files = []
  if os.path.exists(folder_path):
    for filename in os.listdir(folder_path):
      files.append(filename)
    return files
  else:
    print("Folder not found")
    sys.exit()

#directory_list(ROOT_files_path)
  
def average_afterpulse_rate(date, time):

  #creating the file path
  file_path = "ROOT_files/1400V/"
  file_name = date + "_A1400_B1400_" + time + "_output.root"

  #Opening the root file
  root_file  = ROOT.TFile(file_path + file_name, "READ")
  root_file.cd()
  
  #PMT names 612 control, 607 exposed
  #pmts = ["GAO607","GAO612"]
  
  #Getting information from the control pmt
  control_hist = root_file.Get(date + "_GAO612_apulse_num_1400V" )
  control_number = 0
  for i in range(2, control_hist.GetNbinsX()):
    control_number += control_hist.GetBinContent(i)
  
  control_rate = (control_number/control_hist.GetEntries()) * 100
  control_error = (np.sqrt(1/control_rate + 1/control_hist.GetEntries())) * control_rate

  #Getting information from the exposed pmt
  exposed_hist = root_file.Get(date + "_GAO607_apulse_num_1400V" )
  exposed_number = 0
  for i in range(2, exposed_hist.GetNbinsX()):
    exposed_number += exposed_hist.GetBinContent(i)
    
  exposed_rate = (exposed_number/exposed_hist.GetEntries()) * 100
  exposed_error = (np.sqrt(1/exposed_rate + 1/exposed_hist.GetEntries())) * exposed_rate
  
  return control_rate, control_error, exposed_rate, exposed_error

def graph(dates, cr, cr_e, er, er_e, name):
  
  start_date = np.min(dates)
  dates = dates - start_date
  
  plt.errorbar(dates, cr, yerr=cr_e, linestyle='none', color="red")
  control = plt.scatter(dates, cr, color="red")
  
  plt.errorbar(dates, er, yerr=er_e, linestyle='none', color="blue")
  exposed = plt.scatter(dates, er, color="blue")
  
  plt.legend([control,exposed],['Control rate','Exposed rate'])
  
  plt.xlim(0, np.max(dates) +1)
  
  date = str(start_date)
  plt.xlabel("Days since {}/{}/{}".format(date[4:6],date[2:4],date[0:2]))
  plt.ylabel("Average afterpulse rate /%")
  
  plt.savefig(name+'.png')
  print("Image created")
  plt.close()
  

def main(folder_path, image_name):
  files = directory_list(folder_path)
  dates = []
  control_rates = []
  control_errors = []
  exposed_rates = [] 
  exposed_errors = []
  
  for i_file in files:
    split = i_file.split("_")
    date, time = split[0], split[3]
    control_rate, control_error, exposed_rate, exposed_error = average_afterpulse_rate(date, time)
    
    dates.append(int(date))
    control_rates.append(control_rate)
    control_errors.append(control_error)
    exposed_rates.append(exposed_rate)
    exposed_errors.append(exposed_error)
  
  graph(dates, control_rates, control_errors, exposed_rates, exposed_errors, image_name)
  
  
#main("ROOT_files/1400V/")  

if __name__=="__main__":
  if len(sys.argv)!=3:
    print("Must enter folder path and output image name")
    print("Example 'python root_reader.py nemo/rootfiles/ final_image'")
  else:
    folder_path = sys.argv[1]
    image_name = sys.argv[2]
    main(folder_path, image_name)
  
  
  
  
  
  
  
  
  
  
    
