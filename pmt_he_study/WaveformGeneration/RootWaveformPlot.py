import numpy as np
import ROOT
import matplotlib.pyplot as plt
import sys

#Opening the root file
file_path = "RootFiles/"
file_name = "200101WaveModified2"

def main(file_path,file_name, output_file):
  
  file = ROOT.TFile(file_path + file_name,"READ")
  tree = file.T
  file.cd()
  for event in tree:
    wave_array = []
    mf_shapes = []
    mf_amps = []
    
    for value in event.waveform:
      wave_array.append(value)
    for value in event.mf_shapes:
      mf_shapes.append(value)
    for value in event.mf_amplitudes:
      mf_amps.append(value)
    
    time_holder = np.arange(len(wave_array))
    
    filler = np.zeros(39)
    mf_shapes = np.append(mf_shapes, filler)
    mf_amps = np.append(mf_amps,filler)
    
    print(len(time_holder))
    print(len(wave_array))
    print(len(mf_shapes))
    print(len(mf_amps))
    
    if len(time_holder) == 0:
      print("Empty set")
      break
    
    #plt.plot(time_holder,wave_array,linewidth=0.5)
    fig = plt.figure(figsize = (7,10))
    
    plt.subplot(2,1,1)
    plt.title("")
    plt.xlabel("Time placeholder")
    plt.ylabel("Voltage /mV", color="blue")
    plt.plot(time_holder,wave_array,linewidth=0.5, color="blue")
    plt.ylim(925,1000)
    plt.xlim(1900,2100)
    
    ax2 = plt.twinx()
    plt.hlines(0.85,xmin=0,xmax=7000, color="green")
    ax2.plot(time_holder,mf_shapes,linewidth=0.5, color="red")
    ax2.set_ylabel("Shape index", color="red")
    
    
    plt.subplot(2,1,2)
    plt.plot(time_holder,wave_array,linewidth=0.5, color="blue")
    plt.title("")
    plt.xlabel("Time placeholder")
    plt.ylabel("Voltage /mV", color="blue")
    plt.ylim(925,1000)
    plt.xlim(1900,2100)
    
    ax2 = plt.twinx()
    plt.hlines(1,xmin=0,xmax=7000,color="green")
    ax2.plot(time_holder,mf_amps,linewidth=0.5, color="red")
    ax2.set_ylabel("Amplitude index", color="red")
    ax2.set_ylim(0,75)
    
  
    
    plt.tight_layout()
  
    plt.savefig(output_file + '.png')
  
    break
    
if __name__=="__main__":
  if len(sys.argv)!=3:
    print("Must enter config path then output path")
    print("Example 'python RootWaveformPlot.py /RootFiles/test1.root FirstPlot'")
  else:
    input_path = sys.argv[1]
    output_path = sys.argv[2]

    inp_array = input_path.split("/")
    
    file_name = inp_array[len(inp_array)-1]
    
    del inp_array[len(inp_array)-1]
    file_path = "/".join(inp_array) + "/" 
    
    main(file_path, file_name, output_path)