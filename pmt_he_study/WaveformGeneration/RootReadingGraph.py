import ROOT
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import fnmatch
from datetime import date
import pandas as pd
import matplotlib as mpl

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

#path is RandomApPlacementSim/RootFiles/ROOT_files/1400V

def get_pulses(path,filename):
  #creating the file path
  full_path = path + filename
  
  #Opening the root file
  root_file  = ROOT.TFile(full_path, "READ")
  root_file.cd()
  
  hist = root_file.Get("110011_GAO607_apulse_times_1400V")
  total = hist.GetEntries()
  
  return total

import plotly.express as px

def read_matrix(path):
  df = pd.read_csv(path, index_col = 0)
  return df

def main(path,output):
  
  '''
  files = directory_list(path)
  n_pulses = []
  pulses_found = []
  dummy = np.arange(0,19,1)
  indices = ["i0","i1","i2","i3","i4","i5","i6","i7","i8","i9","i10","i11","i12","i13","i14","i15","i16","i17","i18"]
  columns = dummy.astype(str)
  df = pd.DataFrame(index = indices, columns = columns)
  df = df.fillna(0)
  
  
  offset = []
  amps = []
    

  for ifile in files:
    broken = ifile.split("_")
    inserted = broken[4]
    n_pulses.append(int(inserted))
    
    #offset.append(int(broken[4]))
    #amps.append(int(broken[5]))
    
    
    ap_num = int(get_pulses(path,ifile))
    pulses_found.append(ap_num)
    
    #if ap_num == 2 or ap_num == 0:
    #  print(ap_num)
    #  print(ifile)
    
    df.loc["i" + inserted, str(ap_num)] = df.loc["i"+inserted, str(ap_num)] + 1
  
  
  #print(df)
  
  
  print("")
  #print(df.sum(numeric_only=True))
  #df.sum(axis=1) 
  df["total"] = df.sum(axis=1)
  #df["total"] = df.sum(axis=0)
  df.loc["total"] = df.sum()
  #diagonals = np.append(0,np.diag(df,k=0))
  diagonals = np.diag(df)
  df.loc["correct"] = diagonals
  df.loc["accuracy"] = ((df.loc["correct"] / df.loc["total"])*100).astype(int)
  print(df)
  '''
  
  
  #Opening the root file
  root_file  = ROOT.TFile(path, "READ")
  root_file.cd()
  
  hist = root_file.Get("200803_GAO612_apulse_times_1400V")
  p_times = []
  print(hist)
  for i in range(0, hist.GetNbinsX()):
    #print(hist[i])
    #p_times.append(hist.GetBinContent(i))
    p_times.append(hist[i])
  
  p_times = np.array(p_times)
  p_times = p_times/len(p_times)
  
  graph = "hist"
  if graph == "hist":
    bins = np.arange(0,7000,1)
    #hist = plt.hist(p_times,density=True)
    plt.bar(bins,p_times,align='edge',width=1,alpha=1,color='b')
    plt.ylabel("Frequency")
    plt.xlabel("Afterpulse time /ns")
    plt.savefig(output + '.pdf')
  if graph == "yesno":
    fig = px.bar(df, x=df.index, y=columns,barmode='group',
    labels = {"value":"Number of waveforms with given pulses found", "index": "Number inserted", "variable":"Number of   pulses per waveform found"}, width=1000, height=500)
    
    fig.write_image(output + ".png")
  if graph == "normal":
    #plt.plot(n_pulses,n_pulses, linewidth=1, color='red')
    plt.scatter(n_pulses, pulses_found, s=1, color='blue')
    plt.xlabel("Number of afterpulses inserted into each waveform")
    plt.ylabel("Number of afterpulses detected")
    
    plt.savefig(output + ".png")
  
  if graph == "offset":
    fig, ax = plt.subplots()
    ax.set_prop_cycle(color=['red', 'green', 'blue'])
    plot = plt.scatter(offset,amps, c=pulses_found, s=3)
    cbar = fig.colorbar(plot)
    cbar.set_label('Number of pulses found', rotation=90)
    plt.xlabel("Offset /ns")
    plt.ylabel("Amplitude /mV")
    
    plt.savefig(output + ".png", dpi=500)


if __name__=="__main__":
  if len(sys.argv)!=3:
    print("Must enter folder path and output image name")
    print("Example 'python root_reader.py nemo/rootfiles/ final_image'")
  else:
    folder_path = sys.argv[1]
    image_name = sys.argv[2]
    main(folder_path, image_name)
