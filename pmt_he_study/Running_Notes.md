This TXT file describes the method and commands for the PMT permeation data aquisition.

The first thing to do is set up your remote connection to the two PC's in F9. There is the DAQ computer aka the windows PC and the labsupernemo1 PC which controls the HV. Then you will be guided in setting the experiment up, then eventually to the taking of data.

Step 0
------
Note: skip this step if you have done it already. You should only have to do it once.

You need to set up in your own dedicated area on the HEP cluster a working directory containing a clone of the supernemopmts repository. Enter the following commands:

	$ ssh -XYC $USERNAME$@plus2.hep.ucl.ac.uk
		your password
	$ mkdir PMT_Helium_WD
	$ cd PMT_Helium_WD
	$ git clone https://USERNAME:PASSWORD@github.com/willquinn/supernemopmts.git
	
Note that you will need permission to be able to do this so email william.quinn.14@ucl.ac.uk . Run all analysis code form this directory. if you want to make changes to the code, make your own branch and request a merge.

Step 1
------
This step involves using a dedicated gmail account hpge.hep.ucl@gmail.com. The password for this account is vm6sqg5pu4. In a chrome browser go to the website https://remotedesktop.google.com/access/ and log in with the hpge gmail account. There should be a list of the remote devices. The computer which you will have access to  is dw-lab-hep. The pin for this is 607612. You will then connect to the PC and eventually the desktop will appear in your browser. Log into the hep_lab user (if not already logged on) with the password 1st_floor. 

Step 2
------
CAENScope is a GUI oscilloscope. It isn't the best but once the settings are set it will work fine. The icon should be on the dock though can be easily found. When open there will be a screen with a moving crosshair. Go to 'file' and select 'connect'. A pop up box will appear. Check that is says usb and then click ok. This will then connect you to the digitizer which is hooked up to both PMTs. The digitizer should always be on and if disconnected will mean that CAENScope will not see it.

Step 3
------
Execute the following commands from your own computer terminal:

	$ ssh -XYC $USERNAME$@plus2.hep.ucl.ac.uk 
		your own password
	$ ssh -XYC labgate
		your own password
	$ ssh -XYC supernemo@labsupernemo1.hep.ucl.ac.uk
		Password: crapL0g0!
	$ ssh -X admin@192.168.0.1
		Password: admin

Once in you should see a very minimal display with three words (which are actually drop down menus) 'Main' 'Utility' and 'Setup'. Using the arrow keys:

	- Navigate to 'setup' and hit enter. 

It will then show a drop down menu with 'Transparent Mode'

	- Hit enter again. 

A pop up box within your terminal will display "Slot to Communicate"

	- Type 11 then hit enter twice. 

If nothing happens hit enter again. Again if noting happens, try naviagting to Setup again with the arrow keys and retype 11. If it works you should see a screen with A1535 at the top (within the terminal). Where the cursor is highlighted:

	- type D. 

This will bring you to the HV control "centre". The channels we are focusing on are convieniently named PMT_HE_1 and PMT_HE_2. As summary of the labels and which PMT refer to is at the bottom of this doc. Use the arrow keys to navigate across columns and rows within the terminal.

Step 4
------
Note: If the PMTs are both set to 1000V, check the current is around 124 microAmps and then skip this step. If not then do the following for both PMT's. 

	- Set the voltage in the column V0set to 100V and hit enter. 
	- In the column labelled Pw hit the spacebar to ramp up the voltage.

Check the current and voltage readings shown in the Vmon and Imon columns and wait for the voltage to stabalise about the value you have set. The resistance across the PMT base is roughly 8MOhm so the monitored values should be consistent with this so around 12 microAmps. If there is a current going through the board, that means the PMTs are set up properly and that the circuit is continuous. 

	- Now ramp the voltage up to 500V. 

We now need to check the current to see if the resistance is consistant. At this high a voltage if there are any light leaks, this should be easily detected by an excess in current. The current should be around 62.5 microAmps. If the numbers do not match to 8MOhm then ramp the PMT voltage down as there may be something wrong. 

	- Increase the voltage slowly (in 100V steps) to 1000V
	
(For subsequent steps, and once it is established that the PMTs are stable, it is not necessary to always change the voltages in 100V increments). Once the PMT is at 1000V and the currents are consistent (around 124 microAmps) then the PMTs are fine and ready for data-taking.

Step 5
------
Go back to CAENScope. The settings are a bit tricky to set up but fortunately we can use the RunSettings found in the CAENScopeData folder on the Desktop. We will start with channel 0 and then (when told to later on) we will repeat for channel 1. Please amend the file names for the corresponding channel you are recording.

	- On the GUI click the "Restore Settings" button and navigate to the CAENScopeData folder in the I: drive. 
	- Select the RunSettings_Ch0.xml file and click open. (Or RunSettings_Ch1.xml for channel 1)

For the resolution study we will be collecting data for each channel individually. Which settings file you select will depend on the channel you are running i.e RunSettings_Ch0 is for channel 0.

	- Check that the correct channel has been selected under trigger settings and channel enables. 
	- The threshold for both PMT's should be 960. 
	- The DC offset should be 10920 for channel 0 and 11790 for channel 1. 

After doing these checks, we are now ready to record some data:

	- Go to file -> "Record (text)". 

There will be a pop up window for where to save the file. 

	- First create (if it doesn't exist) and then navigate to a directory with the name format YYMMDD (this should be in "I:CAENScopeData").
	- Save the file as A1000_B1000_tHHMM_Ch0.xml 

where HHMM is the time which you start the run. **Note the unusual date format for the directory (YYMMDD). This is needed for subsequent data processing.**

	- Tick the box that says 'Stop After' in the bottom right hand corner.
	- Enter the number 90000 in the box.
	- Click on the up and down arrow keys a few times but make sure to end up on 90000  (this is the only way I have found to set the stop correctly).
	- Press the Start button in the bottom right of CAENScope. This starts the run.
	- Make sure you can see some wavforms. Rotate the 'Position' knob under the horizontal tab in the bottom left of the panel to the position '1' to make the waveforms visible.
	
You should now see waveforms bouncing around on the screen. If not then continue to adjust the horizontal position knob until you find the waveforms. The peak of the waveform should come around 600ns. If this is not the case then there is an issue and you should stop the run and seek advice.

Step 6
------
After an hour the run should be complete. Redo step 5 for the Channel 1.

Step 7
------
Go back to the terminal with the HV control centre on it. 
	
	- Set the voltage on both PMTs to be 1400V.

Step 8
------
On CAENScope:

	- Select the restore settings and choose the RunSettings.xml file in the CAENScopeData directory.

This will set up CAENScope for both channels for the afterpulse DAQ. Again, check to be sure the settings are correct as in step 5. 

	- Select file - record (text) and save in the same dated folder A1400_B1400_tHHMM.xml. 
	- Tick the box that says 'Stop After' int he bottom right hand corner
	- Enter the number 50000 in the box
	- Click on the up and down arrow keys a few times but make sure to end up on 200000
	- Select start on the GUI.
	- Rotate the the 'Position' knob under the horizontal tab in the bottom left to the position '1' to see the waveforms.

Step 9
------
The run should be finished after 4 hours.  

Step 10
-------
With the run now done. Go back to the HV control terminal

	- Set the voltages on both channels int he terminal to 1000.
	
To turn off the PMTs which should only be done if neccessary (at the end of the full experiment or in emergency) then

	- In the Pw column press the spacebar to ramp down the voltage.
	- Set the voltages on both channels int he terminal to 0 as well.

Step 11
-------
In a new terminal, type the following commands
	
	$ ssh -XYC $USERNAME$@plus2.hep.ucl.ac.uk
		your password
	$ ssh -XYC labgate
		your password
	$ cd /unix/nemo3/PMT_He_Study/data/raw_xml_files 
	$ scp -r hep_lab@DW-HEP-LAB:J:\/CAENScopeData\/YYMMDD .

Use password 1st_floor. Note that this copies the whole directory and its contents; you do not need to create an empty directory first. Be sure that the permissions on your directories are such that anyone can read the data files.

Step 12
-------
Run the area analysis code. The permisions are now set within the relevent directories (/home/wquinn/PMT_Helium/) so that no matter where you run the code (in your own cloned repository is preferrable) you can create and ammend the output files. My repository can be found in /home/wquinn/PMT_Helium/repo/supernemopmts. All output .root files will be placed in a dedicated ROOT_files directory with the naming as follows:

	Areas:
		YYMMDD_A1000_B1000_tHHMM_ChX_areas.root
		...
	Apulses:
		YYMMDD_A1400_B1400_tHHMM_analysis.root

Make sure you have an up to date repository and navigate to it on the HEP cluster. You can use pc202 as a dedicated pc for running the code. Run the following commands:
	
	$ git pull
	$ root.exe
	$ gSystem->CompileMacro("Area_Spectrum.cpp","O")
	$ Area_Spectrum("HHMM","X","YYMMDD")
	
Where you have to input the correct values of HHMM (Hour and minute), X (Channel) and YYMMDD that are in the name of the .xml file you have created.
	
After 5 mins or so two TCanvas will pop up with a histogram and a fitted red curve. Then repeat this step for the other channel (i.e 0 or 1). If the fit is poor or the spectrum looks bizare then we have a problem. Type the following commands:

	$ cat /home/wquinn/PMT_Helium/Results/resolution_vs_date_Ch0.dat
	$ cat /home/wquinn/PMT_Helium/Results/resolution_vs_date_Ch1.dat

Check that there is an entry for todays date. If the results are bizare then you can access the TCanvas .pdf in the directory /unix/nemo3/PMT_He_Study/ROOT_files/Area_Spectrum_Files/PDFs/

Step 13
-------
Run the afterpulse analysis code. Again, navigate to the repository (either your own or the one in /home/wquinn/PMT_Helium/repo/supernemopmts. Run the following command:

	$ root.exe
	$ gSystem->CompileMacro("sweep_shape_index.cpp","O")
	$ sweep_shape_index(false,"HHMM","YYMMDD")
	
Then:

	$ cat /home/wquinn/PMT_Helium/Results/apulseNUM_vs_date_Ch0.dat
	$ cat /home/wquinn/PMT_Helium/Results/apulseNUM_vs_date_Ch1.dat

Check that there is an entry for todays date.

Step 14
-------

The .root files you have created are in /unix/nemo3/PMT_He_Study/ROOT_Files/ . You need to set the permissions for other people to be able to use the files. Type something like:

	$ chmod a+rwx /Area_files/*
