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
  
  #print(file_name)
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

def linear_fit(days, cr):
  #Control rate LOBF
  #calculation of gradient of line of best fit
  a= (days - np.mean(days))*cr
  b= (days-np.mean(days))**2
  m1 = np.sum(a)/np.sum(b)
  
  #calculation of y intercept of line
  y_bar= np.mean(cr)           #ybar
  x_bar= np.mean (days)              #xbar
  c1= y_bar - m1*x_bar
  
  #generating the line of best fit
  y = m1*days +c1
  
  #Calculating errors in m and c
  d= cr - m1*days - c1
  D = np.sum(b)
  
  #uncertainty in gradient 
  u_m1= round(np.sqrt((1/D) * (np.sum(d**2))*1/(len(days)-2)),5)
  
  #uncertainty in intercept
  u_c1= round(np.sqrt(((1/len(days)) + (x_bar**2/D)) * ((np.sum(d**2))/(len(days)-2))),3)
  
  return y

def second_poly_fit(days, cr):
  coef = np.polyfit(days,cr,2) 
  y = coef[0]*days**2 + coef[1]*days + coef[2]*days**0
  return y
  
  print(co)
  
def graph(days, cr, cr_e, er, er_e, name, st_date, times):
  
  #fit = "Linear"
  fit = "Second order polynomial"
  
  #error_correction = "ratio"
  error_correction = "flat"
  
  
  start_date = date(int("20"+st_date[0:2]), int(st_date[2:4]), int(st_date[4:6]))
  days_of_ten = date(2020,2,11) - start_date
  
  days = np.array(days)
  cr = np.array(cr)
  er = np.array(er)
  
  #er_e = np.array(er_e)
  #cr_e = np.array(cr_e)
  
  if fit == "Linear":
    y = linear_fit(days, cr)
  else:
    y = second_poly_fit(days,cr)

  rms = np.sqrt(np.mean(np.square(cr-y)))
  
  residuals = cr - y
  
  #error = sqrt(olderror**2 + rms**2)
  
  if error_correction == "ratio":
    control_ratio = rms/cr
    
    corrected_er_e = np.sqrt(np.square(er_e) + np.square(control_ratio*er))
    corrected_cr_e = np.sqrt(np.square(cr_e) + np.square(control_ratio*cr))

  else:
    corrected_er_e = np.sqrt(np.square(er_e) + rms**2)
    corrected_cr_e = np.sqrt(np.square(cr_e) + rms**2) 
  
  
  fig = plt.figure(figsize = (7,10))
  
  plt.subplot(3,1,1)
  plt.title("Raw data with fitting from 1% helium flow")
  plt.xlabel("Days since {}/{}/{}".format(st_date[4:6],st_date[2:4],st_date[0:2]))
  plt.ylabel("Average afterpulse rate /%")
  plt.errorbar(days, cr, yerr=cr_e, linestyle='none', color="red", capsize=1, elinewidth=0.1,fmt='r.', label="Control Rate")
  plt.errorbar(days, er, yerr=er_e, linestyle='none', color="blue", capsize=1, elinewidth=0.1,fmt='b.',label="Exposed Rate")
  plt.plot(days, y,'-', color="black",label="{} fitting".format(fit))
  plt.axvline(days_of_ten.days, label="10% Helium", color = "orange")
  plt.legend(loc="upper left")
  
  plt.subplot(3,1,2)
  plt.hist(residuals, bins=16)
  plt.title("Hisogram spread of residuals, RMS value: {}".format(round(rms,3)))
  plt.xlabel('Residual value')
  plt.ylabel("Frequency")
  
  plt.subplot(3,1,3)
  plt.title("{} adjusted error bars from {} fitting".format(error_correction,fit))
  plt.xlabel("Days since {}/{}/{}".format(st_date[4:6],st_date[2:4],st_date[0:2]))
  plt.ylabel("Average corrected afterpulse rate /%")
  plt.errorbar(days, cr, yerr=corrected_cr_e, linestyle='none', color="red", capsize=1, elinewidth=0.1,fmt='r.', label="Control Rate")
  plt.errorbar(days, er, yerr=corrected_er_e, linestyle='none', color="blue", capsize=1, elinewidth=0.1,fmt='b.',label="Exposed Rate")
  plt.plot(days, y,'-', color="black",label="{} fitting".format(fit))
  plt.axvline(days_of_ten.days, label="10% Helium", color = "orange")
  plt.legend(loc="upper left")
  
  
  plt.tight_layout()
  
  plt.savefig(name+ "_" + fit +"_" + error_correction +'.png')
  print("Image created")
  plt.close()
  
#Main running of the program
def main(folder_path, image_name):
  files = directory_list(folder_path)
  del files[0]
  dates = []
  times = []
  control_rates = []
  control_errors = []
  exposed_rates = [] 
  exposed_errors = []

  for i_file in files:
    split = i_file.split("_")
    date, time = split[0], split[3]
    
    control_rate, control_error, exposed_rate, exposed_error = average_afterpulse_rate(folder_path,date, time)
    
    #Discarding invalid solutions
    if exposed_rate > 0.1:  
      if int(date) > 191226:
        dates.append(int(date))
        time = time[1:]
        times.append(int(time))
        control_rates.append(control_rate)
        control_errors.append(control_error)
        exposed_rates.append(exposed_rate)
        exposed_errors.append(exposed_error)
      else:
        None
    else:
      None
  
  start_date = str(np.min(dates))
  days = day_extractor(dates)
  
  graph(days, control_rates, control_errors, exposed_rates, exposed_errors, image_name, start_date, times)
  
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
