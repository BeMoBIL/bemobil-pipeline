# BeMoBIL BIDS tools  

## Before you try to convert .xdf to BIDS...  

The .xdf to BIDS conversion scripts introduced here are intended to be used internally within the Berlin Mobile Brain Body Imaging Lab.
As much as I want to make it more general, it is inevitable becasue it relies on how channels and streams in the .xdf files are named and this may vary from setup to setup.
The diversity in naming convention for lsl outlets can even be an issue within the group, so one should think this through at the time of implementation and keep each other informed. 
The hope is still that these scripts can be used for other setups that use .xdf with some adjustments. 

Some caution regarding dependencies 

-FieldTrip     - currently data2bids and xdf2fieldtrip versions used here are not in the standard release (now they are [here]( https://github.com/sjeung/fieldtrip/tree/motion2bids)).  
-natsortorder  - a tool to aid in sorting files according to the natural order 
  
As of now the components in BeMoBIL BIDS tool are as follows

-bemobil_xdf2bids.m  
  this is the main function that calls other configuration scripts below
-bemobil_bidsconfig_general.m
  contains configurations that apply to all participants and all modalities
-bemobil_bidsconfig_participant.m
  contains participant or file-specific configuration 
-bemobil_bidsconfig_motion.m
  contains motion specific configuration  
-bemobil_bidsconfig_eeg.m
  contains eeg specific configuration

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
       bemobil_config.rigidbody_streams        = {'playerTransform','playerTransfom','rightHand', 'leftHand', 'Torso'};
       **bemobil_config.bids_rbsessions        = [1,1; 1,1; 0,1; 0,1; 0,1]; 
       **bemobil_config.eeg_stream             = {'BrainVision'};
                                                  a unique keyword used to identify the eeg stream in the .xdf file             
       **bemobil_config.bids_tasklabel         = 'VNE1';

Here the starred fields are used speficially for bids processing.
So, no need to specify them if you go the direct path from .xdf to .set!

## Multi-session, multi-run cases

Some caution is requireed when the source .xdf files are split into multiple runs. 
At the moment, we consider these runs to be split parts of a continuous recording. 
So those run files in a single session will be merged downstream. 

## Selecting streams by other methods

You may also want to select streams using the range of effective sampling rate, if your stream names are tricky to work with.

It is also possible to use numerical indices of the streams. 
However, this seems to differ from session to session. 

If you do NOT specify any stream, it will automatically assume that the stream of the highest sampling rate is eeg and will also try to import all other streams to be put in a single file.
This might be fine if you only have EEG data in your .xdf file but complicates things downstream when you have more modalities. 

## For other modalities than EEG and motion 

