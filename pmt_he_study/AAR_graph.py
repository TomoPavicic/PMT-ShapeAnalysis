import ROOT
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import fnmatch
from datetime import date

def directory_list(folder_path):
  files = []
  if os.path.exists(folder_path):
    for filename in os.listdir(folder_path):
      if fnmatch.fnmatch(filename, '*1400*'):
        files.append(filename)
    return files
  else:
    print("Folder not found")
    sys.exit()
  
def average_afterpulse_rate(file_path, date, time):

  #creating the file path
  file_name = date + "_A1400_B1400_" + time + "_output.root"
  
  print(file_name)
  #Opening the root file
  root_file  = ROOT.TFile(file_path + file_name, "READ")
  root_file.cd()
  

  #PMT names 612 control, 607 exposed
  #pmts = ["GAO607","GAO612"]
  
  #Getting information from the control pmt
  control_hist = root_file.Get(date + "_GAO612_he_apulse_num_1400V" )
  
  #A check to see if the histogram has the valid info
  try:
     control_hist.GetNbinsX()
  except:
    return 0,0,0,0
  
  control_number = 0
  for i in range(2, control_hist.GetNbinsX()):
    control_number += control_hist.GetBinContent(i)
  
  control_rate = (control_number/control_hist.GetEntries()) * 100
  #control_error = (np.sqrt(1/control_rate + 1/control_hist.GetEntries())) * control_rate
  control_error = np.sqrt((control_rate/100) * (1-control_rate/100) / control_hist.GetEntries()) * 100

  #Getting information from the exposed pmt
  exposed_hist = root_file.Get(date + "_GAO607_he_apulse_num_1400V" )
  exposed_number = 0
  for i in range(2, exposed_hist.GetNbinsX()):
    exposed_number += exposed_hist.GetBinContent(i)
    
  exposed_rate = (exposed_number/exposed_hist.GetEntries()) * 100
  #exposed_error = (np.sqrt(1/exposed_rate + 1/exposed_hist.GetEntries())) * exposed_rate
  exposed_error = np.sqrt((exposed_rate/100) * (1-exposed_rate/100) / exposed_hist.GetEntries()) * 100
  
  return control_rate, control_error, exposed_rate, exposed_error

def graph(days, cr, cr_e, er, er_e, name, date):
  
  plt.errorbar(days, cr, yerr=cr_e, linestyle='none', color="red", capsize=1, elinewidth=0.1)
  control = plt.scatter(days, cr, color="red",s=1)
  
  plt.errorbar(days, er, yerr=er_e, linestyle='none', color="blue", capsize=1, elinewidth=0.1)
  exposed = plt.scatter(days, er, color="blue", s=1)
  
  plt.legend([control,exposed],['Control rate','Exposed rate'])
  
  plt.xlim(0, np.max(days) + 5)
  
  plt.xlabel("Days since {}/{}/{}".format(date[4:6],date[2:4],date[0:2]))
  plt.ylabel("Average afterpulse rate /%")
  
  plt.title("Average He afterpulse rate for both PMTs with binomial error")
  
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
    control_rate, control_error, exposed_rate, exposed_error = average_afterpulse_rate(folder_path,date, time)
    
    #Discarding invalid solutions
    if date != 0 or control_rate != 0 or exposed_error != 0:  
      dates.append(int(date))
      control_rates.append(control_rate)
      control_errors.append(control_error)
      exposed_rates.append(exposed_rate)
      exposed_errors.append(exposed_error)
    else:
      None
  
  start_date = str(np.min(dates))
  days = day_extractor(dates)
  
  graph(days, control_rates, control_errors, exposed_rates, exposed_errors, image_name, start_date)
  

def day_extractor(dates):
  start = str(np.min(dates))
  start_date = date(int("20"+start[0:2]), int(start[2:4]), int(start[4:6]))
  days = []
  
  for i_date in dates:
    i_date = str(i_date)
    current_date = date(int("20"+i_date[0:2]),int(i_date[2:4]),int(i_date[4:6]))
    difference = current_date - start_date
    days.append(difference.days)
  
  return days


if __name__=="__main__":
  if len(sys.argv)!=3:
    print("Must enter folder path and output image name")
    print("Example 'python root_reader.py nemo/rootfiles/ final_image'")
  else:
    folder_path = sys.argv[1]
    image_name = sys.argv[2]
    main(folder_path, image_name)

  
  
  
  
  
  
  
  
  
    
