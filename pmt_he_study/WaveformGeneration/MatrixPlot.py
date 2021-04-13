import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

def main(file_path, image_name):
  df = pd.read_csv(file_path + ".txt", index_col = 0)
  df = df/500
  df.loc["i0","0"] = 1
  #M = df.values
  #M = M.T
  df = df.transpose()
  
  #fig = plt.subplots()
  p = plt.pcolor(df, cmap='coolwarm')
  cbar = plt.colorbar(p)
  cbar.set_label("Proportion found",rotation=90)
  plt.xlabel("Number of afterpulses inserted in each waveform")
  plt.ylabel("Number of afterpulses found in each waveform")
  plt.yticks(np.arange(0.5, len(df.index), 1), df.index)
  plt.xticks(np.arange(0.5, len(df.index), 1), df.index)
  plt.savefig(image_name + ".pdf")
  
if __name__=="__main__":
  if len(sys.argv)!=3:
    print("Must enter input matrix.txt and output file name")
    print("Example 'python MatrixPlot.py Matrix1 Plot1'")
  else:
    file_path = sys.argv[1]
    image_name = sys.argv[2]
    main(file_path, image_name)