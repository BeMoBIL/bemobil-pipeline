# BeMoBIL BIDS tools  

## To convert .xdf to BIDS...  

The .xdf to BIDS conversion scripts introduced here are intended to be used internally within the Berlin Mobile Brain Body Imaging Lab. As much as we want to make it more general, it is inevitable becasue it relies on how channels and streams in the .xdf files are named and this may vary from setup to setup. The diversity in naming convention for lsl outlets can even be an issue within the group, so one should think this through at the time of implementation and keep each other informed.  The hope is still that these scripts can be used for other setups that use .xdf with some adjustments. 

Not that, although the EEG part is meant to pass the BIDS validator already, the motion part is not included in the current BIDS version and will change along with the specs for motion data. 


As of now the components in BeMoBIL BIDS tool are as follows

- **bemobil_xdf2bids.m**  
  this is the main function that calls other configuration scripts below  
- **bemobil_bidsconfig_general.m**   
  contains configurations that apply to all participants and all modalities  
- **bemobil_bidsconfig_participant.m**  
  contains participant or file-specific configuration   
- **bemobil_bidsconfig_motion.m**   
  contains motion specific configuration    
- **bemobil_bidsconfig_eeg.m**  
  contains eeg specific configuration  
- **bemobil_bidsmotion.m**
  deals with quat2eul conversion which is necessary for BIDS conformity


Dependencies

- **FieldTrip**     
currently data2bids and xdf2fieldtrip versions used here are not in the standard release  
(now they are [here]( https://github.com/sjeung/fieldtrip/tree/motion2bids)).  

- **natsortorder**  
a tool to aid in sorting files according to the natural order 
Stephen Cobeldick (2021). [Natural-Order Filename Sort](https://www.mathworks.com/matlabcentral/fileexchange/47434-natural-order-filename-sort), MATLAB Central File Exchange. Retrieved March 1, 2021.
  

There is no dependency on the rest of the BeMoBIL pipeline as long as the configuration is correct. 
For any questions/suggestions regarding BeMoBIL BIDS tools pleaase don't hesitate to contact <seinjeung@gmail.com>.

# How to use 

## A very standard usage 

You can use the example scripts by minimal changes using the input "bemobil_config" struct. 
Below is the list of fields in bemobil_config that are used by the batch of scripts 

       bemobil_config.study_folder             = 'E:\Project_BIDS\example_dataset_MWM\';
       bemobil_config.filename_prefix          = 'sub_';
       bemobil_config.raw_data_folder          = '0_raw-data\';
       **bemobil_config.bids_data_folder       = '1_BIDS-data\';     
       bemobil_config.raw_EEGLAB_data_folder   = '2_basic-EEGLAB\';
       bemobil_config.channel_locations_filename = 'VN_E1_eloc.elc'; 
       bemobil_config.filenames                = {'VR' 'desktop'}; 
       bemobil_config.rigidbody_streams        = {'playerTransform','rightHand', 'leftHand', 'Torso'};
       **bemobil_config.bids_rbsessions        = [1,1,1,1 ; 1,0,0,0]; 
       **bemobil_config.eeg_streamkeyword      = {'BrainVision'};
                                                  a unique keyword used to identify the eeg stream in the .xdf file             
       **bemobil_config.bids_tasklabel         = 'VNE1';

Here the starred fields are used speficially for bids processing.
So, you can leave them unspecified if you go the direct path from .xdf to .set!

Now some more detailed description of BIDS specific fields below...


       bemobil_config.bids_data_folder     = '1_BIDS-data\';  

type : STRING  
default value : '1_BIDS-data\'   
This is the folder in which all files are going to be written.
                                          
                                          
       bemobil_config.bids_rbsessions       = [1,1,1,1 ; 1,0,0,0];  

type : LOGICALS of size numel(filenames) X numel(rigidbody_streams)      
default value : []       
This indicates which streams are included in repective recording sessions. For instance, in the example above, playerTransform might only be present in session 'VR', so the first row '1,1,1,1' means all rigidbody streams are present in 'VR' session but the '1,0,0,0' in the second row means only the first type of rigidbody is in session 'desktop'.
       
      
      bemobil_config.eeg_streamkeyword     = {'BrainVision'}; 
       
type : Cell  
default value : {'BrainVision'}  
A cell containing the keyword to be used to identify EEG stream


       bemobil_config.bids_tasklabel        = 'VNE1';

type : STRING  
default value : 'taskname'  
Task label to be used in constructing bids filenames with no '_' or '-' character
  
  

Once all the config fields are filled out, you can simply call function **bemobil_xdf2bids.m** with only one additional input, namely the **numericalIDs**. 


        numericalIDS            = [1,2,3,5,10]; 
        bemobil_xdf2bids(bemobil_config, numericalIDs)

The IDs are assumed to be numerical, as it will make things easier when we convert BIDS to .set to use in the pipeline.
But this is of course a redundant restriction if one is only interested in converting .xdf to BIDS. 

## Multi-session, multi-run cases

Some caution is required when the source .xdf files are split into multiple runs.  
At the moment, we consider these runs to be split parts of a continuous recording.  
So those run files in a single session will be first saved as different run files in BIDS but will be merged when you import from BIDS to .set using BeMoBIL tools. 


Entries in **bemobil_config.filenames** will search through the raw data directory of the participant and group together .xdf files with matching keyword in the name into one session. If there are multiple files in one session, they will be given separate 'run' numbers in the file name

         
The order of runs rely on incremental name sorting using the ["Nature Order File Sorting Tool"](https://www.mathworks.com/matlabcentral/fileexchange/47434-natural-order-filename-sort)
for example,

                   sub-1\sub-1_VNE1_VR_rec1.xdf
                   sub-1\sub-1_VNE1_VR_rec2.xdf
                   sub-1\sub-1_VNE1_VR_rec11.xdf
                   sub-1\sub-1_VNE1_desktop.xdf
                   
will be organized into

                   sub-001\ses-VR\sub-001_ses-VR_task-VNE1_run-1_eeg.bdf
                   sub-001\ses-VR\sub-001_ses-VR_task-VNE1_run-2_eeg.bdf
                   sub-001\ses-VR\sub-001_ses-VR_task-VNE1_run-11_eeg.bdf
                   sub-001\ses-desktop\sub-001_ses-desktop_task-VNE1_eeg.bdf
               
and then be merged into one .set file later in that order (*'run-1, run-2, run-11'*), instead of (*'run-1, run-11, run-2'*).

This is why it is important to name the .xdf files correctly. You would not want the merged file to have random order.
This below is a bad example as the order is 1) not clear to human readers and 2) can be shuffled in the merged file

                   sub-1\sub-1_VNE1_VR_first.xdf  % don't do this!      
                   sub-1\sub-1_VNE1_VR_old1.xdf   % don't do this!
                   sub-1\sub-1_VNE1_VR_test.xdf   % don't do this!

## Selecting streams by other methods

To use these methods you would have to modify the scripts on your own. They are just briefly described here. 

- effective sampling rate 
You may also want to select streams using the range of effective sampling rate, if your stream names are tricky to work with.

- numerical indices of the stream
It is also possible to use numerical indices of the streams. However, this seems to differ from session to session. 

- automatic 
If you do NOT specify any stream, it will automatically assume that the stream of the highest sampling rate is eeg and will also try to import all other streams to be put in a single file. This might be fine if you only have EEG data in your .xdf file but complicates things downstream when you have more modalities. 

## For other modalities than EEG and motion 

There is a plan to also support some other modalities that are often used in MoBI research. 
