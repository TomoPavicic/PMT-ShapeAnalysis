import ROOT
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import fnmatch
from datetime import date

def directory_list(folder_path):
  Ch0_files = []
  Ch1_files = []
  if os.path.exists(folder_path):
    for filename in os.listdir(folder_path):
      if fnmatch.fnmatch(filename, '*Ch0*'):
        Ch0_files.append(filename)
      elif fnmatch.fnmatch(filename, '*Ch1*'):
        Ch1_files.append(filename)
    
    del Ch0_files[0]
    del Ch1_files[0]
    
    return Ch0_files,Ch1_files
  else:
    print("Folder not found")
    sys.exit()
  

def graph(e_days, e_rates, e_errors, c_days,c_rates,c_errors,name,date):
  
  plt.errorbar(c_days, c_rates, yerr=c_errors, linestyle='none', color="red", capsize=1, elinewidth=0.1)
  control = plt.scatter(c_days, c_rates, color="red",s=1)
  
  plt.errorbar(e_days, e_rates, yerr=e_errors, linestyle='none', color="blue", capsize=1, elinewidth=0.1)
  exposed = plt.scatter(e_days, e_rates, color="blue", s=1)
  
  plt.legend([control,exposed],['Control rate','Exposed rate'])
  
  plt.xlim(left=0)
  plt.ylim(0,1)
  
  plt.xlabel("Days since {}/{}/{}".format(date[4:6],date[2:4],date[0:2]))
  plt.ylabel("Average afterpulse rate /%")
  
  plt.title("Average afterpulse rate for both PMTs at 1000V")
  
  plt.savefig(name+'.png')
  print("Image created")
  plt.close()
  

def AAR(path, channel, date, time):
  if channel == "Ch0":
    PMT = "GAO607"
  if channel == "Ch1":
    PMT = "GAO612"

  #creating the file path
  file_name = date + "_A1000_B1000_" + time + "_" + channel +"_output.root"

  #Opening the root file
  root_file  = ROOT.TFile(path + file_name, "READ")
  root_file.cd()
  
  hist = root_file.Get(date + "_" + PMT +"_apulse_num_1000V" )

  try:
    if hist.GetEntries() == 0:
      return 0,0
  except:
    return 0,0
    
  number = 0
  for i in range(2, hist.GetNbinsX()):
    number += hist.GetBinContent(i)
  
  rate = (number/hist.GetEntries()) * 100
  if rate > 2:
    print(file_name)
    print("rate: {}".format(rate))
    print("number: {}".format(number))
    print("Entries: {}".format(hist.GetEntries()))
  #error = (np.sqrt(1/rate + 1/hist.GetEntries())) * rate
  error = np.sqrt((rate/100) * (1-rate/100) / hist.GetEntries()) * 100
    
  return rate, error

def AAR_files(path, files, channel):
  dates = []
  rates = []
  errors = []
  for i_file in files:
    split = i_file.split("_")
    date, time = split[0], split[3]
    rate, error = AAR(path, channel,date,time)  
    if rate != 0:
      dates.append(int(date))
      rates.append(rate)
      errors.append(error)
    
  return dates, rates, errors

def ratio_graph(day, ratio):
  plt.scatter(day,ratio)
  plt.xlabel("Days since 08/10/19")
  plt.ylabel("Ratio of Ch1 afterpulse rate to Ch0 afterpulse rate")
  plt.title("Ratio between channels")
  plt.ylim(0.5,2)
  plt.savefig('RatioPlot2 zoomed.png')
  print("Image created")
  plt.close()

#Main running of the program
def main(folder_path, image_name):
  Ch0,Ch1 = directory_list(folder_path)
  
  Ch0_dates = []
  Ch0_rates = []
  Ch0_errors = []
  
  Ch0_dates, Ch0_rates, Ch0_errors = AAR_files(folder_path, Ch0, "Ch0")
  
  Ch1_dates, Ch1_rates, Ch1_errors = AAR_files(folder_path, Ch1, "Ch1")
  
  #print(Ch0_dates)
  #print(type(Ch0_dates[0]))
  start_date = str(np.min(Ch0_dates))
  Ch0_days = day_extractor(Ch0_dates)
  Ch1_days = day_extractor(Ch1_dates)
  
  
  day_ratio = []
  ratio = []
  for i in range(0,len(Ch0_days)):
    for j in range(0,len(Ch1_days)):
      if Ch0_days[i] == Ch1_days[j]:
        day_ratio.append(Ch0_days[i])
        ratio.append(Ch1_rates[j]/Ch0_rates[i])
  
  ratio_graph(day_ratio, ratio)
  
  
  #graph(Ch0_days,Ch0_rates, Ch0_errors,Ch1_days,Ch1_rates, Ch1_errors, image_name, start_date)
  
#Takes the array of dates and returns the amount of days since the beginning date
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


  