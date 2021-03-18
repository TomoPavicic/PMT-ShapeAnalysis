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

ATCs = ["ATC25","ATC30","ATC35","ATC40"]

def average_afterpulse_rate(file_path, date, time):

  #creating the file path
  file_name = date + "_A1400_B1400_" + time + "_output.root"
  
  #Opening the root file
  root_file  = ROOT.TFile(file_path + file_name, "READ")
  root_file.cd()
  
  #Getting information from the control pmt
  control_hist = root_file.Get(date + "_GAO612_he_apulse_num_1400V" )
  
  #A check to see if the histogram has the valid info
  try:
     control_hist.GetNbinsX()
  except:
    return 0,0,0,0
  
  #Getting information from the exposed pmt
  exposed_hist = root_file.Get(date + "_GAO607_he_apulse_num_1400V" )
  exposed_number = 0
  for i in range(2, exposed_hist.GetNbinsX()):
    #exposed_number += exposed_hist.GetBinContent(i)
    exposed_number += i*exposed_hist.GetBinContent(i)
  
  #exposed_rate = (exposed_number/exposed_hist.GetEntries()) * 100
  #exposed_error = np.sqrt((exposed_rate/100) * (1-exposed_rate/100) / exposed_hist.GetEntries()) * 100
  
  #exposed_rate = exposed_number/exposed_hist.GetEntries()
  #exposed_error = np.sqrt((exposed_rate/100) * (1-exposed_rate/100) / exposed_hist.GetEntries())
  
  exposed_rate = exposed_number
  exposed_error = np.sqrt((exposed_rate/100) * (1-exposed_rate/100) / exposed_hist.GetEntries())
  
  return exposed_rate, exposed_error

def exposed_rate(ATC,files):
  
  exposed_rates = [] 
  exposed_errors = []
  
  folder_path = ATC + "/ROOT_files/1400V/"
  
  for i_file in files:
    split = i_file.split("_")
    date, time = split[0], split[3]
    
    exposed_rate, exposed_error = average_afterpulse_rate(folder_path,date, time)
    print(folder_path)
    #Discarding invalid solutions
    if date != 0 or control_rate != 0 or exposed_error != 0:  
      exposed_rates.append(exposed_rate)
      exposed_errors.append(exposed_error)
    else:
      None
  
  return exposed_rates, exposed_errors

def day_extractor(dates):
  start = "191107"
  start_date = date(int("20"+start[0:2]), int(start[2:4]), int(start[4:6]))
  days = []
  
  for i_date in dates:
    i_date = str(i_date)
    current_date = date(int("20"+i_date[0:2]),int(i_date[2:4]),int(i_date[4:6]))
    difference = current_date - start_date
    days.append(difference.days)
  
  return days


def main():
  files = directory_list("ATC25/ROOT_files/1400V/")
  dates = []
  
  for i_file in files:
    split = i_file.split("_")
    date, time = split[0], split[3]
    dates.append(int(date))
  
  date = "191107"
  days = day_extractor(dates)
  
  fig, ax1 = plt.subplots()
  plt.title("Multiplicity for different afterpulse time cuts")
  ax1.set_xlabel("Days since {}/{}/{}".format(date[4:6],date[2:4],date[0:2]))
  ax1.set_ylabel("Afterpulse multiplicity")
  
  fmts = ["r.","b.","g.","y."]
  count = 0
  
  for ATC in ATCs:
    rates, errors = exposed_rate(ATC,files)
    
    ax1.errorbar(days, rates, yerr=errors, linestyle='none', capsize=1, elinewidth=0.1,fmt=fmts[count], label=ATC)
    count += 1 
    print(rates)
    
  plt.legend(loc='upper left')
  
  #plt.xlim(380,384)
  #plt.ylim(60,64)
  
  plt.savefig("MultiplicityTimeCutsNewConf.png")
  print("Image created")
  plt.close()
  

main()

















  