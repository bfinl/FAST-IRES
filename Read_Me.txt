Read Me
This readme document describes software codes and testing human data, which have been developed under support from the National Institutes of Health via grants NS096761, EB021027, MH114233, AT009263 awarded to Dr. Bin He. The testing human data were collected as part of NIH funded research at Mayo Clinic, Rochester, overseen by Dr. Greg Worrell. 
The source codes and testing human data are provided as service to the scientific community and may be used for any non-commercial purpose.  Users should use the codes or data at their own risks.
If anyone uses the whole or part of the codes or the human data, provided here, we ask you to cite the following publication (https://doi.org/10.1038/s41467-020-15781-0) in any of your publications or presentations:

Sohrabpour et al, “Noninvasive Electromagnetic Source Imaging of Spatio-temporally Distributed Epileptogenic Brain Sources,” Nature Communications, PAGE NUMBER, 2020.
------------------------------------------------------------------------------------------------------------------------------------------

This folder contains the codes and sample data used to analyze EEG data within the FAST-IRES framework. De-identified data of one patient containing inter-ictal spikes recorded in EEG, as well as one seizure recording is provided.
The same codes provided here, were run and the results presented here are derived with this exact code, to ensure accuracy.
Solutions to the spike and seizure imaging analysis, are saved in folders named Results and Figures. Spike results/figures are stored in a sub-folder named “Results\Spike_I” and “Figures\Spike_I”.
To run these codes, you need to have MATLAB installed on your computer (the current code has been tested on MATLAB 2018a and 2018b - most codes developed on the 2013b version, so most probably this code is compatible with prior versions). The current version of the code was run on MATLAB 2018a/b and topo plots are also depicted using EEGLAB's codes. You can download and deploy EEGLAB toolbox from https://sccn.ucsd.edu/eeglab/index.php. Please use latest version of EEGLAB (version 14.1.1b recommended - also tested on version 14.1.2b).
A function called "norms" is used in these codes which is borrowed from CVX package (http://cvxr.com/cvx/). We included this function and (its nested functions) in Codes folder, for your convenience. CVX is a MATLAB toolbox that can be installed easily and freely. For the purposes of running our codes, its installation is not needed.
For connectivity analysis, the latest version of eConnectome must be installed, which is also a freely available MATLAB toolbox (https://www.nitrc.org/projects/econnectome/) as well as arfit (http://climate-dynamics.org/software/#arfit). It is safe to add EEGLAB, eConnectome, and ARFit to the MATLAB path, to avoid any errors.

Now we will quickly go through the folders and codes and how to run the provided example.

Folders:
- Codes: Contains the main codes needed to solve the inverse problem for seizure imaging, spike imaging and also codes necessary to perform the connectivity analysis, i.e. connectivity imaging.
- Connectivity : A folder to save Connectivity results, i.e. the outputs after running code "Patch_Time_Course_Extractor.m".
- Denoised Seizures : Contains the denoised seizure, seizure time-basis functions (TBFs), etc. This is the denoised seizure data which is the output after running this code "TBF_Selection_Seizure.m". If you don't want to run the pre-processing again, you can directly load these files into the main code for seizure analysis "FIRES_Sz_Rot.m". We recommend, loading the pre-loaded parameter files for best performance.
- Extracting EZ from FAST-IRES Estimates: Contains codes that extract the estimated EZ from spike, seizure, connectivity imaging analysis, to ultimately assess how good we match the clinical data. Not critical for running the imaging codes. Clinical findings may contain sensitive personal data and cannot be shared publicly. This patient became seizure-free after surgery (ILAE 1). The surgery was a left anterior temporal lobectomy and amygdalohippocampectomy. SOZ defined from iEEG was limited to the left anterior portion of the temporal lobe and some depth electrodes on the mesial side near hippocampus.
- Figures: Figures of the estimated sources (sum of all components, i.e. columns of J) - Spike results are in the sub-folder called Spike_I and the seizure images are in the root, figures are automatically saved when running "FIRES_Sz_Rot.m" and "FIRES_Spk_Rot.m" for seizures and spikes, respectively.
- Grid Location and Defined Parameters: Folder containing lead-field, mesh data, etc. for this patient. Use as is. Data was partially analyzed with a commercial package called Curry 8 (Compumedics, NC, USA).
- Raw 273: contains the raw spike and seizure data (in case you want to run the pre-processing yourself) and a generic electrode location file. 
- Results: Saved results/solutions of the estimated sources, i.e. columns of J, etc. - Spike results are in the sub-folder called Spike_I and the seizure results are in the root, results are automatically saved when running "FIRES_Sz_Rot.m" and "FIRES_Spk_Rot.m" for seizures and spikes, respectively.
- Spikes: Contains the denoised spikes and parameters, spike TBFs, etc. This is the denoised spike data which is the output after running this code "TBF_Spike_Extraction.m". If you don't want to run the pre-processing again, you can directly load these files into the main code for spike analysis "FIRES_Spk_Rot.m". We recommend, loading the pre-loaded parameter files for best performance.

Simplest way to run the code.
Seizure Imaging/Connectivity Imaging:
- First run the main seizure imaging code "FIRES_Sz_Rot.m" in Codes folder,
- Then either run "Patch_Time_Course_Extractor.m" for connectivity analysis or just go to next step to skip connectivity analysis for simplicity,
- Then look at results, figures, etc. or go to the "Extracting EZ from FAST-IRES Estimates" folder and run "Frequency_Peak_Solution.m" or "Connectivity_Segment_Solution.m" code to see the final estimations of the EZ and performance, for seizure and connectivity imaging results, respectively.

Spike Imaging:
- First run the main seizure imaging code "FIRES_Spk_Rot.m" in Codes folder,
- Then look at results, figures, etc. or go to the "Extracting EZ from FAST-IRES Estimates" folder and run "Spike_Peak_Solution.m" code to see the final estimations of the iterative zone and performance, for spike imaging results. 

A bit about the codes.
Codes Folder:
- Main Seizure Imaging Code is "FIRES_Sz_Rot.m" - You can open and see all the parameters, etc. This solves the seizure imaging problem and saves the figures and all the spatio-temporal results.
- Main Spike Imaging Code is "FIRES_Spk_Rot.m" - You can open and see all the parameters, etc. This solves the spike imaging problem and saves the figures and all the spatio-temporal results. 
Both, the main spike and seizure imaging codes, call the solvers FISTA_ADMM_IRES.m and Proj_Flat_Hyper_Ellips.m for some optimization procedures.
- Main Connectivity code is "Patch_Time_Course_Extractor.m", performs the spatio-temporal correlation analysis to find the hyper-nodes and performs Granger causality. Calls upon the DTF... files (from eConnectome) as well as the code Find_Patch, to find the piecewise homogeneous steps of the FAST-IRES solution.
- Pre-processing codes are "TBF_Selection_Seizure.m" and "TBF_Spike_Extraction.m" for seizures and spikes, respectively. You can skip these as the results are already saved and some modules such as “runica” may change from run to run or computer to computer, needing re-tuning and calibration.
- "norms.m", "cvx_check_dimension.m" and "cvx_default_dimension.m" are the codes from CVX, included here, for your convenience. Just perform vector norm operation over columns of a matrix. Nothing too complicated.

This is FAST-IRES in a nutshell.
Abbas Sohrabpour
March 20th 2020

