import ROOT
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import fnmatch
from datetime import date
import pandas as pd
from scipy.special import factorial
from scipy.stats import chisquare, moyal
from scipy.optimize import curve_fit
from datetime import date
from matplotlib import gridspec

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


def poisson(x,mu,A):
  factorials = factorial(x)
  return A*np.exp(-mu)*mu**x/factorials

def gaussian(x,mu,sd,A):
  exponent = -0.5*np.square((x-mu)/sd)
  return A*np.exp(exponent)

def gaus_mod(x,mu,sd,A):
  diff = x - mu
  exponent = np.exp( -0.5*(diff/sd)**2 )
  #const = 1/(sd*np.sqrt(2*np.pi))
  #dist = const*np.exp(exponent)
  dist = A*exponent
  
  
  return dist_smear(dist)

def poisson_mod(x,mu,A):
  factorials = factorial(x)
  dist = A*np.exp(-mu)*mu**x/factorials
  return dist_smear(dist)
  
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
  
def landau_mod(x,loc,sf):
  dist = moyal.pdf(x,loc,sf)
  return dist_smear(dist)


def day_plot(path,name):
  
  files = directory_list(path)
  del files[0]

  dates =[]
  chi_values = []

  mean_diff = []
  
  for i_file in files:
    split = i_file.split("_")
    date = split[0]
    if "." in date:
      None
    else:
    
      dates.append(int(date))
      
      root_file  = ROOT.TFile(path+i_file, "READ")
      root_file.cd()
      
      hist = root_file.Get(date + "_GAO607_apulse_num_1400V")
      print(date)
      mean = hist.GetMean()
      
      bins = []
      freq = []
      
      for i in range(1, hist.GetNbinsX()):
         bins.append(int(i-1))
         freq.append(int(hist.GetBinContent(i)))
  
    
      bins = np.array(bins)
      freq = np.array(freq)
      normed = freq/sum(freq)
  
      chivalue = np.round(mychi(normed,poisson(bins,mean),np.sqrt(freq)/sum(freq)),5)
  
      chi_values.append(chivalue)
      
      popt, pcov = curve_fit(poisson,bins,normed)
      fitted_mean = popt[0]
  
      mean = hist.GetMean()
      
      mean_diff.append(abs(fitted_mean-mean))
      
  st_date = str(np.min(dates))
  days = day_extractor(dates)
  
  #start_date = date(int("20"+st_date[0:2]), int(st_date[2:4]), int(st_date[4:6]))
  
  plt.scatter(days,mean_diff, s=1)
  #plt.ylabel(r"$\chi^2_{R}$")
  plt.ylabel("Difference in fitted mean and ROOT mean")
  plt.xlabel("Days since " + st_date)
  
  plt.title("Difference in fitted Poisson mean and ROOT mean over time")
  plt.savefig(name+".png")


def single_hist(path, name):
  filename = path.split("/")[-1]
  date = filename.split("_")[0]
  
  root_file  = ROOT.TFile(path, "READ")
  root_file.cd()
  
  hist = root_file.Get(date + "_GAO607_apulse_num_1400V")
  mean = hist.GetMean()

  bins = []
  freq = []
  values = []
  
  for i in range(1, hist.GetNbinsX()):
     bins.append(int(i-1))
     freq.append(int(hist.GetBinContent(i)))

  bins = np.array(bins)
  freq = np.array(freq)
  normed = freq/sum(freq)
  
  popt, pcov = curve_fit(poisson,bins,normed)

  x = np.arange(0,20,0.1)
  
  fitted_mean = popt[0]
  fitted_y = poisson(x,fitted_mean)
  
  chivalue = np.round(mychi(normed,poisson(bins,mean),np.sqrt(freq)/sum(freq)),5)
  
  plt.errorbar(bins+0.5,normed, xerr=0.5, yerr=np.sqrt(freq)/sum(freq), fmt='.', color="b")
  plt.scatter(bins+0.5,normed, alpha =1, color='b',label="Real data", s=10)
  
  placement=0.0000001
  
  plt.text(x=13,y=placement,s=r"$\chi^2_{R}$ :"  + str(chivalue))
  plt.text(x=13,y=0.1*placement, s="Mean from ROOT: " + str(round(mean,2)))
  plt.text(x=13,y=0.01*placement, s="Mean from Curve_fit: " + str(round(fitted_mean,2)))
  
  plt.bar(bins,poisson(bins,mean),align='edge',width=1,alpha=0.5,color='r', label="Poisson distribution")
  plt.plot(x+0.5,poisson(x,mean), alpha=1, color='r')
  plt.plot(x+0.5,fitted_y, alpha=1, color='g', label="Curve fit Poisson")
  
  plt.yscale('log')
  plt.legend()
  plt.savefig(name + ".png")
  

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

def get_pulses(path,filename):
  #creating the file path
  full_path = path + filename
  
  #Opening the root file
  root_file  = ROOT.TFile(full_path, "READ")
  root_file.cd()
  
  hist = root_file.Get("110011_GAO607_apulse_times_1400V")
  total = hist.GetEntries()
  
  return total

def get_bins_values(path):
  filename = path.split("/")[-1]
  date = filename.split("_")[0]
  
  root_file  = ROOT.TFile(path, "READ")
  root_file.cd()
  
  hist = root_file.Get(date + "_GAO607_apulse_num_1400V")
  mean = hist.GetMean()

  bins = []
  freq = []
  values = []
  
  for i in range(1, hist.GetNbinsX()):
     bins.append(int(i-1))
     freq.append(int(hist.GetBinContent(i)))

  bins = np.array(bins)
  freq = np.array(freq)
  normed = freq/sum(freq)
  
  return bins, normed, mean, freq

def dist_smear(start_values):
  df = pd.read_csv("Matrix1.txt", index_col = 0)
  df = df/500
  df.loc["i0","0"] = 1
  M = df.values
  transformed_values = M.T @ start_values
  return transformed_values

def mean_calc(bins,dist):
  dot_prod = np.dot(bins,dist)
  return dot_prod/sum(dist)


def main(path, name):
  graph = "transform"
  fit = "gaus"
  
  bins, real_values, real_mean, freq = get_bins_values(path)
  
  #real_mean = 2.6
  #sd = 2.3
  

  freq = [0.000001 if x==0 else x for x in freq]
  print(freq)
  #freq = freq.astype(float)
  #freq[13:] = 0.0000001

  err = np.sqrt(freq)/sum(freq)
  
  x = np.arange(0,20,0.1)
  
  #err = 0.1 * real_values
  #err[13:] = 0.0001
  
  if fit == "gaus":
    #popt, pcov = curve_fit(gaussian,bins,real_values)

        
    popt,pcov = curve_fit(gaus_mod,bins,real_values ,sigma=err)
    real_mean = popt[0]
    sd = popt[1]
    A = popt[2]
    
    start_values = gaussian(bins,real_mean,sd,A)
    continuous_plot = gaussian(x,real_mean,sd,A)
    
    title = "Starting Gaussian mean: {}, sd: {}, and Constant: {}".format(round(real_mean,2),round(sd,2),round(A,2))
    
  if fit == "poisson":
    #popt, pcov = curve_fit(poisson,bins,real_values)
    popt,pcov = curve_fit(poisson_mod,bins,real_values, sigma = err)

  
    real_mean = popt[0]
    A = popt[1]
    start_values = poisson(bins,real_mean,A)
    continuous_plot = poisson(x,real_mean,A)
    
    title = "Starting Poisson mean: {}, A: {}".format(round(real_mean,2),round(A,2))
    
  if fit == "landau":
    popt,pcov = curve_fit(landau_mod,bins,real_values,sigma=err)
    
    loc = popt[0]
    sf = popt[1]
    
    start_values = moyal.pdf(bins,loc,sf)
    continuous_plot = moyal.pdf(x,loc,sf)
    title = "Approximated Landau Distribution"

  
  
  transformed_values = dist_smear(start_values)
  
  trans_residuals = real_values - transformed_values
  start_residuals = real_values - start_values
  
  chivalue = np.round(mychi(real_values,transformed_values,err),5)
  
  data_mean = round(mean_calc(bins,real_values),3)
  dist_mean = round(mean_calc(bins,start_values),3)
  smear_mean = round(mean_calc(bins,transformed_values),3)
  
  print(data_mean)
  print(dist_mean)
  print(smear_mean)
  
  fig, (ax0, ax1) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]})
  
  ax0.set_title(title)
  ax0.errorbar(bins+0.5,real_values, xerr=0.5, yerr=err, fmt='.', color="b")
  ax0.scatter(bins+0.5,real_values, alpha =1, color='b',label="Real data", s=10)
  ax0.annotate(r"$\chi^2_{R}$ :"  + str(chivalue), xy=(1,0.6), xycoords="axes fraction",
  horizontalalignment ="right")
  ax0.annotate("Data mean: {}".format(data_mean), xy=(1,0.5), xycoords="axes fraction", horizontalalignment ="right")
  ax0.annotate("Distribution mean: {}".format(dist_mean), xy=(1,0.4), xycoords="axes fraction", horizontalalignment ="right")
  ax0.annotate("Smeared mean: {}".format(smear_mean), xy=(1,0.3), xycoords="axes fraction", horizontalalignment ="right")
  
  
  ax0.bar(bins,transformed_values,align='edge',width=1,alpha=0.5,color='r', label="Transformed Distribution")
  ax0.plot(x+0.5,continuous_plot, alpha=1, color='g',label="Starting Distribution")
  ax0.set_yscale('log')
  ax0.set_ylabel("Normalised frequency")
  #ax0.set_xlabel("Number of afterpulses")
  ax0.legend()
  
  ax1.scatter(bins+0.5,start_residuals, label="Starting Distribution", color="g", s=5)
  ax1.errorbar(bins+0.5,trans_residuals, xerr=0.5, yerr=err, fmt='.', color="r",alpha=0.5)
  ax1.scatter(bins+0.5,trans_residuals, label="Transformed", color ="r", alpha =0.5,s=5)
  ax1.set_ylabel("Residuals")
  ax1.set_xlabel("Number of afterpulses")
  ax1.hlines(0,xmin=0,xmax=18, color="y", linestyles="dashed")
  ax1.legend(loc="lower right")
  #ax1.yscale("log")
    
  fig.tight_layout()
  fig.savefig(name + ".png")  


  if graph == "single":
    single_hist(path,name)
  
  if graph == "all_days":
    day_plot(path,name)


if __name__=="__main__":
  if len(sys.argv)!=3:
    print("Must enter folder path and output image name")
    print("Example 'python root_reader.py nemo/rootfiles/ final_image'")
  else:
    folder_path = sys.argv[1]
    image_name = sys.argv[2]
    main(folder_path, image_name)