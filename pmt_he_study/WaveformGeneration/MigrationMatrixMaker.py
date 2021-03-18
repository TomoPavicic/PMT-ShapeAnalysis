import sys
import ROOT
import numpy as np
import pandas as pd
import os
import fnmatch

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

def get_pulses(path,filename):
  #creating the file path
  full_path = path + filename
  
  #Opening the root file
  root_file  = ROOT.TFile(full_path, "READ")
  root_file.cd()
  
  hist = root_file.Get("110011_GAO607_apulse_times_1400V")
  total = hist.GetEntries()
  
  return total

def read_matrix(path):
  df = pd.read_csv(path, index_col = 0)
  return df
  
def main(path, output):
  
  
  files = directory_list(path)
  n_pulses = []
  pulses_found = []
  dummy = np.arange(0,19,1)
  indices = ["i0","i1","i2","i3","i4","i5","i6","i7","i8","i9","i10","i11","i12","i13","i14","i15","i16","i17","i18"]
  columns = dummy.astype(str)
  df = pd.DataFrame(index = indices, columns = columns)
  df = df.fillna(0)
  
    

  for ifile in files:
    broken = ifile.split("_")
    inserted = broken[4]
    n_pulses.append(int(inserted))

    ap_num = int(get_pulses(path,ifile))
    pulses_found.append(ap_num)

    df.loc["i" + inserted, str(ap_num)] = df.loc["i"+inserted, str(ap_num)] + 1

  df.to_csv(output + ".txt", sep=",")


if __name__=="__main__":
  if len(sys.argv)!=3:
    print("Must enter input folder path and output file name")
    print("Example 'python MigrationMatrixMaker.py WaveformGeneration/VaryingAp/ROOT_files table1'")
  else:
    folder_path = sys.argv[1]
    image_name = sys.argv[2]
    main(folder_path, image_name)