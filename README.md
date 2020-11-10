# PMT-ShapeAnalysis

Set Up
------

This repo contains a Class file for reading any form of PMT waveform.

Requirements to use this code:
  - Must have python3
  - Must have pyroot
  
I include here, as it may be useful to others a quick tutorial in how to set up these two things. However, I will only do so for Mac OS as I have only figured it out for such.

First get Homebrew installed - https://brew.sh

Then install python3 and pip3:

  	$ brew install python3
  	$ brew install pip3
  
Now we have to install root. This is not simple as if you run brew install root6, it will work, but not for python3. Brew uses as default python2 and when you build root in this manner, importing it into python3 doesn't work. I you want to use python2 then you have no issue.

We have to build root from source. Just incase, you'll need:

  	$ brew install cmake
  	$ brew install make
  
Now go to root https://root.cern.ch/releases and get the latest release. Download this file to /usr/local and unzip. you willneed root priveleges.

Type the following:

  	$ sudo mkdir root6
  	$ cd root6
  
The next one is important you know which version of python you have downloaded (type python3 --version). Also change $Directory$ in the next line to the unzipped file you got from root-cern website:

  	$ sudo cmake ../$Directory$ -DPYTHON_EXECUTABLE=/usr/local/bin/python3 -DPYTHON_INCLUDE_DIR=/usr/local/Cellar/python/3.7.6_1/Frameworks/Python.framework/Versions/3.7/lib/libpython3.7.dylib -DPYTHON_LIBRARY=/usr/local/Cellar/python/3.7.6_1/Frameworks/Python.framework/Versions/3.7/lib/libpython3.7.dylib
  	$ sudo make
  
The last line should take a while and you will see a percentage number as the first character in each line.

Now in your .bash_profile add:

  	$ cd /usr/local/root6
  	$ source bin/thisroot.sh
  	$ cd $HOME

Once that is done you can now run python and import ROOT. And you shouldn't get any error.

Usage
-----
Now you should be able to run the test.py file. This will read a sample set of waveforms I have provided and plot some of these data to the screen and save a few files to a test root file. If this works then everything is running fine. The output file will be in the directory that you clone this repo and I advise you delete it.

The main file to focus on in this repo is the PMT_Classes.py file. The other files are examples of how I am using this code. The code is designed to be universal so that it can read any form of waveform.

With python3, pyroot installed properly and a good IDE you should be able to familiarise yourself with all the functionality of this trivial OOP for PMTs. I have tried to write a function for everything and have made the naming almost too clear so there should be no problem.
  
For those with access to the HEP cluster at UCL you can use pc204 or pc202 as both have python3 and root6 set up ready to use - as well as many other modules that this repo uses such as minidom.

I would recommend familiarising oneself with test.py as this file give a crash course in how to use some features of the class structure.

TMany of the projects in this repository revolve around understanding the  shape of PMT pulses. I have created, though not necessarily utilised fully, three object classes to clean up some of the code. Also helps with reusability.

PMT_Object.py
-------------
Class for each PMT that will hold settings and some preliminary result histograms. To set up you need to give it a unique name. Only add to these histograms or settings if you know you will utilise them regularly.

This has a set of default settings:

        - charge_cut = 6  # pC
        - charge_range = [0, 100]
        - nbins = 100
        - amp_range = [0, 100]
        - mf_shape_range = [-1, 1]
        - mf_amp_range = [0, 100]
        - sweep_range = [0, 500]
        - pulse_time_threshold = 100
        - apulse_region = 500
        - resistance = 50 # Ohms
        - mf_shape_threshold = 0.9
        - mf_amp_threshold = 25
        - baseline = 1000
        - waveform_length = 8000
        - trigger_point = 100
        - integration = [0.3, 0.3]
        - he_region = [1400, 2000]

The names of these settings are very important. use the get_setting("name") or set_setting("name", value) to access or change. 

The default histograms that are created are:

    - pulse_charge_hist
    - pulse_amplitude_hist
    - apulse_charge_hist
    - pulse_mf_shape_hist
    - pulse_mf_amplitude_hist
    - pulse_mf_shape_mf_amplitude_hist
    - pulse_mf_shape_p_amplitude_hist
    - pulse_mf_amplitude_p_amplitude
    - pulse_peak_time_hist
    - apulse_times_hist
    
These names are also important. They all have the ranges and number of bins defined in the setting above. See the set_up_histograms function in PMT_Object.py for more details.

These histograms can be access by member function or by get_histogram("name") function. To fill you can do either as well.

PMT_Waveform.py
---------------
Feed a list of strings, and a PMT_Object, and it will convert that list of strings into a PMT_Waveform object and you can fill the PMT_Object with the analysis from the PMT_Waveform.

Requires you to pass it a PMT_Object.

First it checks the charge of the pmt pulse which it finds automatically using the integration settings. It then checks the PMT_Object settings for a charge cut. By default, even after applying a setting file [see later in PMT_Array], the charge cut is the only trigger cut. If the charge of the pulse is bigger than the cut value the pmt waveform will trigger to True. This boolean can be used as a descriminator to fill the results of the analysis into the PMT_Object. This avoids null waveforms for example/ empty waveforms.

The results it calculates by default are:

            "pulse_charge"      : 0.0,
            "pulse_amplitude"   : 0.0,
            "apulse_charge"     : 0.0,
            "mf_pulse_shape"    : 0.0,
            "mf_pulse_amp"      : 0.0,
            "pulse_peak_time"   : self.get_pmt_pulse_peak_position(),
            "apulse_times"      : []
            
These names are again very important. Use get_results("name" )to access these results. I have showed the type of value above so you can be sure to use them correctly.

PMT_Object.py has a function called fill_pmt_hists() which you call after each triggered waveform. See test.py as an example of this.

PMT_Array.py
------------
A class to contain multiple PMT_Objects.

Requires parsing a PMT_Object to the class. Depending on the structure you can access these objects by the number they are associated to or the position (if they are in a 2D array). 

There are several functions in this class that make fill the PMT_Objects easier if you have lots of PMTs as well as accessing their results.

If one is not happy with the default settings you can simply call the function apply_setting("config_file_name") which will read all the settings you want from the file. See example_cofig_file.txt to see the format. The format is quite specific but the example makes  it clear. 
  
As it stands this applies the same setting to all PMTs. This will become a problem when we want to apply different settings to each PMT, so for now you will have to be more specific for special cases.

Extras
------
I have also included some example data reading functions but these I will not explain in detail.

Cluster Running
---------------
For those running on the UCL HEP computer cluster, a combination of numpy, scipy, ROOT and other modules may not be properly set up. I recommend pc202 or pc204 as they seem to work. However, the only requirement you really need is that you pc has ROOT for python3. So the first test is to do the following:

    $ python3
    >>> import ROOT
    
If this works with no error then all is good. If not find a PC that this works on. Type the following in the terminal:

    $ virtualenv --python=python3 --system-site-packages myenv
    $ source myenv/bin/activate
    $ pip install --upgrade pip
    $ pip3 install —-upgrade numpy scipy
    
You are now in a python3 virtual environment. All the code should now work and you can now download any modules that you may want or need. To deactivate the environment type:

    $ deactivate
