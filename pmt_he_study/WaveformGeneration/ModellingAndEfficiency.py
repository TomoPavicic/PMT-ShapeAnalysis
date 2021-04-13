import ROOT
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import fnmatch
from datetime import date
from scipy.optimize import curve_fit
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
  
def get_mean(path):
  filename = path.split("/")[-1]
  date = filename.split("_")[0]
  
  root_file  = ROOT.TFile(path, "READ")
  root_file.cd()
  
  hist = root_file.Get(date + "_GAO607_he_apulse_num_1400V")
  
  """
  IN CASE I NEED FREQUENCIES FOR ERROR
  freq = []
  for i in range(1, hist.GetNbinsX()):
     freq.append(int(hist.GetBinContent(i)))
  freq = np.array(freq)
  """
  return hist.GetMean()

def mychi(O,Exp,Err):
  
  chi = 0
  parameters = 1
  bins = 0
  
  for i in range(len(O)):
    if O[i] == 0:
      None
    else:
      bins = bins + 1
      chi = chi + (O[i]-Exp[i])**2/Err[i]**2
      
  return chi/(bins-parameters)

def efficiency_factor(t,A,B,C):
  exp = np.exp(-C*t)
  return A*(1-1/(1+B*exp))
  
def model(t, p0, p1, p2):
    y = []
    for i in range(len(t)):
        temp = 0
        for n in range(1,11):
            temp += ((-1)**n/n**2)*(1 - np.exp(-(n**2)*(np.pi**2)*t[i]/p1))
        f2 = (2/np.pi**2)*p1*temp
        y.append(p0*(t[i] + f2) + p2)
    return y
'''
def model_eff(t,p0,p1,p2,sf):
    #p0 = 0.004281808926303662
    #p1 = 468.27671443406996
    #p2 = 0.03867020632284615

    y = []
    for i in range(len(t)):
        temp = 0
        for n in range(1,11):
            temp += ((-1)**n/n**2)*(1 - np.exp(-(n**2)*(np.pi**2)*t[i]/p1))
        f2 = (2/np.pi**2)*p1*temp
        y.append(p0*(t[i] + f2) + p2)
    

    eff = sf*efficiency_factor(t,0.9445,47,0.0108)
    
    return y*eff
'''

def model_eff(t,p1,p2,sf):
    #p0 = 0.004281808926303662
    #p1 = 468.27671443406996
    #p2 = 0.03867020632284615

    y = []
    for i in range(len(t)):
        temp = 0
        for n in range(1,11):
            temp += ((-1)**n/n**2)*(1 - np.exp(-(n**2)*(np.pi**2)*t[i]/p1))
        f2 = (2/np.pi**2)*p1*temp
        y.append((t[i] + f2) + p2)
    

    eff = sf*efficiency_factor(t,0.9445,47,0.0108)
    
    return y*eff



def main(path, name):
  files = directory_list(path)
  del files[0]
  dates = []
  means = []
  
  for i_file in files:
    split = i_file.split("_")
    i_date = split[0]
    if "." in i_date or i_date == "200317":
      None
    else:
      if int(i_date) < 200211:
        None
      else:
        dates.append(int(i_date))
        means.append(get_mean(path+i_file))
  
  del dates[-1]
  del means[-1]
  
  del dates[-1]
  del means[-1]
  
  
  st_date = str(np.min(dates))
  days = day_extractor(dates)
  
  days = np.array(days)
  means = np.array(means) 
  
  print(days)
  print()
  
  start_date = date(int("20"+st_date[0:2]), int(st_date[2:4]), int(st_date[4:6]))
  days_of_ten = date(2020,2,11) - start_date
  
  #err = 0.01*means
  
  popt, pcov = curve_fit(model, days, means)
  #popt, pcov = curve_fit(model, days, means)
  
  print("OLD MODEL:")
  p0 = popt[0]
  p1 = popt[1]
  p2 = popt[2]
  
  print(p0)
  print(p1)
  print(p2)
  for g in pcov:
    print(g)

    
  days = np.sort(days)
  #order_days = days
  fitting =  model(days, p0,p1,p2)
  '''
  old params
  A = 0.9669579020045902
  B = 6954.156822629187
  C = 0.012268293728577519
  '''
  
  '''
  New He Params
  A = 0.9443040849127243
  B = 46.7234792182298
  C = 0.010837913743546992
  
  eff = efficiency_factor(order_days,A,B,C)
  
  new_model = eff*fitting
  '''
  print("")
  print("")
  popt, pcov = curve_fit(model_eff,days,means)
  #popt, pcov = curve_fit(model_eff,days,means, maxfev=2000)
  p0 = popt[0]
  p1 = popt[1]
  p2 = popt[2]
  #scale = popt[3]
  
  print("New model")
  print(p0)
  print(p1)
  print(p2)  
  #print(scale)

  for g in pcov:
    print(g)

  new_days = np.arange(0,1000,1)
  
  err = means*0.01
  print(mychi(means,fitting,err)/3)
  
  new_model = model_eff(days,p0,p1,p2)
  
  print(mychi(means,new_model,err)/3)
  
  
  residuals = means - fitting
  residuals_new = means - new_model
  
  fig, (ax0, ax1) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]})
  #fig, ax0 = plt.subplots()
  
  ax0.set_title("Permeation modelling with efficiency factor")
  ax0.scatter(days,means, s=1, color='b', label="Data")
  #ax0.errorbar(days,means,yerr=err,color='b', fmt='.', elinewidth=1)
  ax0.plot(days,fitting, alpha=1, color='r',label="Fitting")
  ax0.plot(days,new_model, alpha=1, color='g',label="Fitting with efficiency")
  ax0.axvline(days_of_ten.days, label="10% Helium", color = "orange")
    
  #ax0.annotate("A: "+str(round(A,2)), xy=(0.01,0.9), xycoords="axes fraction",horizontalalignment ="left")

  
  ax0.set_ylabel("Average afterpulse number")
  ax0.legend(loc="upper left")
  #ax0.set_xlabel("Days since " + st_date[4:6] + "/" + st_date[2:4] + "/" + st_date[0:2])
  
  
  #ax1.scatter(days,residuals, color="r", s=1)
  ax1.scatter(days,residuals_new, color="g", s=1)
  ax1.set_xlabel("Days since " + st_date[4:6] + "/" + st_date[2:4] + "/" + st_date[0:2])
  ax1.set_ylabel("Residuals")
  ax1.axhline(0, color="y", linestyle="--")
  

  fig.tight_layout()
  fig.savefig(name + ".pdf")  
  
  
if __name__=="__main__":
  if len(sys.argv)!=3:
    print("Must enter folder path and output image name")
    print("Example 'python ModellingAndEfficiency.py nemo/rootfiles/ final_image'")
  else:
    folder_path = sys.argv[1]
    image_name = sys.argv[2]
    main(folder_path, image_name)
    
  
  
  