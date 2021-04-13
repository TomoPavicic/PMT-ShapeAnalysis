'''
This is a script made to generate waveforms for XML files in the same sytle we get from the CAENScope as 
raw data files. It needs a config file and output path to run. Config file contains parameter instructions.

'''


import configparser
import sys
from xml.dom import minidom
import xml.etree.ElementTree as ET 
import numpy as np
from random import randint
import matplotlib.pyplot as plt
import ROOT

def extract_ap_times():
  #Wills date: 15/05/2020 : 1311
  #path = "/unix/nemo4/PMT_He_Study_nemo4/data/ROOT_files/200823_A1400_B1400_t1832.root"
  output_path = "ApTimeDistribution"
  image = "n"
  
  path = "/unix/nemo4/PMT_He_Study_nemo4/data/ROOT_files/200515_A1400_B1400_t1311.root"
  file = ROOT.TFile(path,"READ")
  tree = file.T
  file.cd()
  
  ap_times = []
  for event in tree:
    for value in event.apulse_times:
      if value > 1400 and value < 2000:
        ap_times.append(value)
  
  
  if image == "y":
    bins = np.arange(0,8000,40)
    hist = plt.hist(ap_times,bins=bins,density=True)
    plt.ylabel("Frequency")
    plt.xlabel("Afterpulse time /ns")
    plt.savefig(output_path + '.png')
  
  if image == "n":
    return ap_times

#Function to generate an array 7168 elements long that resemble afterpulses from a template.
def ap_template(config):
  number = int(config['afterpulses']['number'])
  spacing = config['afterpulses']['spacing_type']
  min_pos = int(config['afterpulses']['min_pos'])
  max_pos = int(config['afterpulses']['max_pos'])
  
  temp_path = config['ap_template']['file_path']
  scale_factor = float(config['ap_template']['scale_factor'])
  
  file = ROOT.TFile(temp_path,"READ")
  file.cd()
  
  #if spacing == "sampled":
  #  ap_distribution = extract_ap_times()
  
  temp = file.Get("Template_Ch0")
  
  all_aps = np.zeros(7168)
  
  #first = True
  
  for i in range(number):
    
    '''
    if first == False:
      scale_factor = 0.00002*int(config['afterpulses']['counter'])
    first = False
    '''
    
    randomN = randint(0,20)
    scale_factor = 0.0002 + 0.00002*randomN
    #config['ap_template']['scale_factor'] = str(scale)
    
    single_ap = np.zeros(7168)
    
    if spacing == 'random':
      position = randint(min_pos,max_pos)
    
    if spacing == 'fixed':
      position = int(config['afterpulses']['position'])
    
    if spacing == 'sampled':
      position = int(ap_distribution_glob[randint(0,len(ap_distribution_glob)-1)])
      
    if spacing == 'spaced':
      position = 2000 + i*int(config['afterpulses']['space_dist'])
    
    for i in range(temp.GetNbinsX()):
      single_ap[position-1+i] = int(temp.GetBinContent(i))
    
    single_ap = single_ap*scale_factor
    all_aps = all_aps + single_ap
  
  return all_aps.astype(int)

#Function to generate an array 7168 elements long that resemble afterpulses from specified properties.
def ap_gen(config):

  number = int(config['afterpulses']['number'])
  spacing = config['afterpulses']['spacing_type']
  min_pos = int(config['afterpulses']['min_pos'])
  max_pos = int(config['afterpulses']['max_pos'])
  
  width = int(config['ap_generation']['width'])
  amp = int(config['ap_generation']['amplitude'])
  offset = int(config['ap_generation']['offset'])
  
  all_aps = np.zeros(7168)
  
  for i in range(number):
    single_ap = np.zeros(7168)
    if spacing == "random":
      position = randint(min_pos,max_pos)
    
    #refer to diagram if confused
    rise_step = amp/(width*0.5 + offset)
    fall_step = amp/(width*0.5 - offset)
    peak_point = width*0.5 + offset
    
    for j in range(1,width+1):
      if j <= peak_point:
        single_ap[position-2+j] = rise_step*j
      else:
        single_ap[position-2+j] = fall_step*(width+1-j)
        
    all_aps = all_aps + single_ap
  
  all_aps = all_aps*-1

  return all_aps.astype(int)

#Function to generate an array 7168 elements long that resemble a baseline waveform from a template.
def base_template(filename):
  f = open(filename,"r")
  line_string = f.readline()
  
  parse = line_string.split(" ")
  final = []
  for i in parse:
    final.append(int(i))
  
  return final

#Function to generate an array 7168 elements long that resemble a baseline waveform from specified properties.
def base_gen(config):
  
  position = int(config['waveform']['peak_pos'])
  amp = int(config['waveform']['peak_amp'])
  width = int(config['waveform']['width'])
  offset = int(config['waveform']['offset'])

  base = np.zeros(7168)
  upper = True
  for i in range(7168):
    random_n = randint(1,100)
    if upper:
      if random_n == 1:
        base[i-1] = 983
      elif random_n <= 52:
        base[i-1] = 982
      else:
        base[i-1] = 981
      upper = False
    else:
      if random_n == 1:
        base[i-1] = 978
      elif random_n <= 42:
        base[i-1] = 979
      else:
        base[i-1] = 980
      upper = True
  
  main_peak = np.zeros(7168)
  
  rise_step = amp/(width*0.5 + offset)
  fall_step = amp/(width*0.5 - offset)
  peak_point = width*0.5 + offset
  
  for i in range(1,width+1):
    if i <= peak_point:
      main_peak[position-2+i] = rise_step*i
    else:
      main_peak[position-2+i] = fall_step*(width+1-i)
  main_peak = main_peak*-1
  
  wave = base + main_peak
  
  return wave.astype(int)

#Main function to create an array with 7168 elements which resembles a baseline with inserted afterpulses
def create_waveform(config):

  base_template_use = config['waveform']['template']
  if base_template_use == 'n':
    wave = base_gen(config)
  if base_template_use == 'y':
    wave = base_template(config['waveform']['file'])
  
  template_use = config['afterpulses']['template']
  if template_use == 'n':
    afterpulses = ap_gen(config)  
  if template_use == 'y':
    afterpulses = ap_template(config)
  
  wave = wave + afterpulses
  
  return wave

#Main program which calls the create_waveform function and converts the array into the output type specified
def main(config, output_path):
  file_type = config['file']['output_type']
  
  waveform = create_waveform(config)
  
  if file_type == "XML" or file_type == "xml":
  
    wave_str = ''
    for i in waveform:
      wave_str = wave_str + str(i) + ' '

    event = ET.Element('event')
    event.set('id','1')
    event.set('settings','22')
    event.set('digitizer','AAA')
    event.set('timestamp','1010')
       
    triggershift = ET.SubElement(event, 'triggershift')
    triggershift.set('samples','0')
    
    trace = ET.SubElement(event,'trace')
    trace.set('channel','0')
    trace.text = wave_str
    
    str_xml = ET.tostring(event)
    pretty_xml = minidom.parseString(str_xml).toprettyxml(indent=" ")
  
    with open(output_path+".xml", "w") as f: 
        f.write(pretty_xml)
        
  if file_type == "image":
    
    time_holder = np.arange(len(waveform))
  
    fig = plt.figure(figsize = (7,10))

    plt.title("")
    plt.xlabel("Time placeholder")
    plt.ylabel("Voltage /mV", color="blue")
    plt.plot(time_holder,waveform,linewidth=0.5, color="blue")
    plt.ylim(925,1000)
    #plt.ylim(977,983)
    plt.savefig(output_path + '.png')
  
  if file_type == 'noise_hist':
    temp_base = waveform[1500:7000]
    print("RMS Value: {}".format(np.sqrt(np.mean(np.square(temp_base)))))
    bin_array = [976,977, 978,979,980,981,982,983,984,985]
    hist = plt.hist(temp_base, bins=bin_array)
    print(hist[1])
    print(hist[0])
    plt.savefig(output_path + '.png')
  
  else:
    None

 
if __name__=="__main__":
  if len(sys.argv)!=3:
    print("Must enter config path then output path")
    print("Example 'python XML_waveform_gen.py /waveforms/setup.conf /waveforms/'")
  else:
    conf_path = sys.argv[1]
    output_path = sys.argv[2]
    
    config = configparser.ConfigParser()
    config.read(conf_path)
    
    loop = config['file']['loop']
    
    n_files = int(config['file']['n_files'])
    
    global ap_distribution_glob
    ap_distribution_glob = extract_ap_times()
    
    #sf = 0.00002  ##~1mv

    
    if loop == 'y':
      '''
      for i in range(n_files):
        #position = ap_distribution[randint(0,len(ap_distribution)-1)]
        #config['afterpulses']['position'] = str(position)
        config['afterpulses']['counter'] = str(i)
        for j in range(-50,51):
          #config['afterpulses']['number'] = str(j)
          #j is the space between
          config['afterpulses']['space_dist'] = str(j)
          
          base = "0000"
          base = base + str(i)
          time_stamp = base[len(str(i)):]
  
          file_name = output_path + "A1400_B1400_t" + time_stamp + "_" + str(j) + "_" + str(i)
          main(config,file_name)
        '''
        
      for i in range(n_files):
        base = "0000"
        base = base + str(i)
        time_stamp = base[len(str(i)):]
        for j in range(1,19):
          config['afterpulses']['number'] = str(j)
          file_name = output_path + "A1400_B1400_t" + time_stamp + "_" + str(j)
          main(config,file_name)
        
        
        '''
        base = "0000"
        base = base + str(i)
        time_stamp = base[len(str(i)):]
  
        file_name = output_path + "A1400_B1400_t" + time_stamp
        main(config,file_name)
        '''
        
  
        print("{}/{} completed".format(i+1,n_files))
        
    else:
      main(config, output_path)
    
    
    
    
    
    
    
    
    
    