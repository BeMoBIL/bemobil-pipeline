# BeMoBIL BIDS tools 1. Source data in .xdf to BIDS

## To convert .xdf to BIDS...  

The .xdf to BIDS conversion scripts introduced here are intended to be used *internally* within the **Berlin Mobile Brain Body Imaging Lab**. On the user side, it will be less work compared to directly using standard conversion scripts (provicded by [FieldTrip](https://github.com/fieldtrip) or [EEGlab](https://github.com/sccn/bids-matlab-tools)), as the tool deals with all the idiosyncrasies of BeMoBIL recording setup. On the other hand this means that it makes a lot of assumptions about how the source data is formatted and collected. Especially processing the motion data relies on how channels and streams in the .xdf files are named. The scripts were written for data set that contains streams collected using HTC Vive for instance. The diversity in naming convention for lsl outlets can even be an issue within the group, so one should think this through *at the time of implementation* and keep each other informed. It is a trade-off between generalizability and ease. Now all we try to put most weight on the "ease" side. The hope is still that these scripts can be used for other setups that use .xdf with some adjustment and become more flexible over time. 

Note that, although the EEG part is meant to pass the BIDS validator already, the motion part is not included in the current BIDS version and will change along with the specs for motion data ([BEP029](https://bids.neuroimaging.io/get_involved.html#extending-the-bids-specification)). 


As of now the components in BeMoBIL BIDS tool are as follows

- **bemobil_xdf2bids.m**  
  this is the main function that calls other configuration scripts below  
- **bemobil_bids_motioncfg.m**   
  default motion specific configuration    
- **bemobil_bids_eegcfg.m**  
  default eeg specific configuration  
- **bemobil_bids_motionconvert.m**  
  deafault function for processing the motion data (e.g., quat2eul conversion)
- **bemobil_bids_parsemarkers.m**  
  deafault function for processing the events 


Dependencies

- **FieldTrip**     
currently data2bids and xdf2fieldtrip versions used here are not in the standard release. 
(now they are [here]( https://github.com/sjeung/fieldtrip/tree/motion2bids))

- **natsortorder**  
a tool to aid in sorting files according to the natural order 
Stephen Cobeldick (2021). [Natural-Order Filename Sort](https://www.mathworks.com/matlabcentral/fileexchange/47434-natural-order-filename-sort), MATLAB Central File Exchange. Retrieved March 1, 2021.
  

It also uses the homemade quat2eul conversion function utils_quat2eul.m.
There is no dependency on the rest of the BeMoBIL pipeline as long as the configuration is correct. 

We do hope to make your life easier, not harder.  
So if you are confused by anything, and have any suggestions regarding BeMoBIL BIDS tools please don't hesitate to contact <seinjeung@gmail.com>.

# How to use 

## Naming and organizing your files

Assuming that you have some experience with BeMoBIL pipelne and know how to organize your data so that the standard mobilab import can deal with it, you are good to go. 
But there is notable difference between BeMoBIL and BIDS conventions in how file names are broken down to different components.

If you have been using the standard import in the pipeline already, your source data would look something like this :   


<img src="/resources/bidstools/unises1.png" width="300">


The components come from bemobil config fields 'filename_prefix' and 'filenames'. 
Note that we are now assuming that there is the entry in 'filenames' is a cell that contains only one name.   
For multi-session and multi-run cases, please scroll down this page and read the corresponding subsection.
This also talks about what to do when one continuous recording has been broken down into multiple files. 


But now moving on with a single-subject, single-file case. 
BIDS has this field called 'task label' that will look like this in the BIDS-formatted file below. 



<img src="/resources/bidstools/unises2.png" width="500">

Now the name is taken from the BIDS-specific config field called 'bids_tasklabel'. This label can differ from the 'filenames' entry, and the 'filenames' entry is not even visible at this step, which can be confusing.  
Entries in field 'filenames' can be interpreted both as a tasklabel in a unisession case and session labels in a multisession case (for example, when there is one VR and one desktop session per participant, filenames are equivalent to session labels 'VR' and 'Destkop' but the task label is 'NavigationTask'. On the other hand, in a unisession case with only one 'filenames' element, this will be used like a task label). This is again, detailed later in the document. Note that, in a single session case, the value of 'bids_tasklabel' will not appear again within the pipeline. However, when you are converting the data set to BIDS in order to share it, we recommend you make it informative enough for people who do not know the data set. For example, using string "VirtualWalkOnStraightLines" instead of "VNE1" is more desirable. 'bids_tasklabel' field exists also because, unlike the pipeline, BIDS does not allow usage of certain characters in tasklabels. For instance, 'VN_E1' is OK to use in the BeMoBIL pipeline but BIDS does not allow 'task-VN_E1'. 


And then, if you later convert the unisession BIDS data to BeMoBIL compatible .set files, it will look like this. 



<img src="/resources/bidstools/unises3.png" width="350">


To sum up, we recommend you to do the following 

> single-session data set : only the single entry in 'filenames' is important, but also try to find a sensible name for 'bids_tasklabel'  
> multi-session data set : use multiple 'filenames' to indicate keywords that represent each session and 'bids_tasklabel' as the name of the task common to all sessions

If you still have questions like "So what about broken recording sessions?" or "How are sessions and runs different and how should that be reflected in my file names?"
Please do scroll down and check out the bit on multi-session and multi-run handing 

## A very standard usage 
### 1. Fill out the relevant fields in bemobil_config struct

Following is a detailed description of the example script 'example_bemobil_bids_xdf2bids.m' 
The first section writes the struct "bemobil_config" with some bids-specific fields.

    bemobil_config.study_folder             = 'M:\8_Conferences\MoBI Workshop\data\';
    bemobil_config.filename_prefix          = 'vp_';
    bemobil_config.raw_data_folder          = '0_raw-data\';
    bemobil_config.filenames                = {'_walk'}; 
    bemobil_config.resample_freq            = 250;
    bemobil_config.rigidbody_streams        = {'PhaseSpace_Rigid1','PhaseSpace_Rigid2', 'PhaseSpace_Rigid3', ...
                                            'PhaseSpace_Rigid4','PhaseSpace_Rigid5', 'PhaseSpace_Rigid6', ...
                                            'PhaseSpace_Rigid7','PhaseSpace_Rigid8', 'PhaseSpace_Rigid9'};
    bemobil_config.channel_locations_filename = [];                            

    % bids fields 
    bemobil_config.bids_data_folder         = '1_BIDS-data\';
    bemobil_config.bids_rbsessions          = true(1,numel(bemobil_config.rigidbody_streams));
    bemobil_config.bids_eegkeyword          = {'BrainVision'}; 
    bemobil_config.bids_tasklabel           = 'walking';
    
    % custom function names - customization highly recommeded, especially when not using HTC-Vive
    bemobil_config.bids_motioncfg_custom        = 'bids_motioncfg_mobiworkshop';
    bemobil_config.bids_motionconvert_custom    = 'bids_motionconvert_mobiworkshop';
    bemobil_config.bids_parsemarkers_custom     = 'bids_parsemarkers_mobiworkshop';

Here the last seven fields are used speficially for bids processing. So, you can leave them unspecified if you go the direct path from .xdf to .set and don't need to use BeMoBIL BIDS tools.

The last three fields are needed when you want to operate on the motion or event data more closely and will be described later in section **Custom functions used for .xdf to BIDS conversion**. Custom functions are essential for motion setups other than HTC Vive and also for parsing events properly. Even when it works without custom functions, it is highly recommened that one looks are each default files (bemobil_bids_motioncfg.m, bemobil_bids_motionconvert.m, bemobil_bids_parsemarkers.m) to double check what happens. 
Now some more detailed description of the four BIDS specific fields below...


       bemobil_config.bids_data_folder     = '1_BIDS-data\';  

type : STRING  
default value : '1_BIDS-data\'   
This is the folder in which all files are going to be written.
                                          
                                          
       bemobil_config.bids_rbsessions       = [1,1,1,1 ; 1,0,0,0];  

type : LOGICALS of size numel(filenames) X numel(rigidbody_streams)      
default value : LOGICAL ones of size numel(filenames) X numel(rigidbody_streams)       
This indicates which streams are included in repective recording sessions. For instance, in the example above, playerTransform might only be present in session 'VR', so the first row '1,1,1,1' means all rigidbody streams are present in 'VR' session but the '1,0,0,0' in the second row means only the first type of rigidbody is in session 'desktop'.
       
      
      bemobil_config.bids_eegkeyword      = {'EEG'}; 
       
type : Cell  
default value : {'EEG'}  
A cell containing the keyword to be used to identify EEG stream.


       bemobil_config.bids_tasklabel      = 'walking';

type : STRING  
default value : 'defaulttask'  
Task label to be used in constructing bids filenames with no '_' or '-' character
  

### 2. Specify general information about the data set

The next step is to specify general metadata about the data set. Every data set in BIDS **must** include this information about the authors, institution, and the task. For details please refer to [BIDS page for modality agnostic files](https://bids-specification.readthedocs.io/en/stable/03-modality-agnostic-files.html).

    generalinfo = [];

    % root directory (where you want your bids data to be saved)
    generalinfo.bidsroot                                = fullfile(bemobil_config.study_folder, bemobil_config.bids_data_folder); 

    % required for dataset_description.json
    generalinfo.dataset_description.Name                = 'Walking task in the young and old';
    generalinfo.dataset_description.BIDSVersion         = 'unofficial extension';

    % optional for dataset_description.json
    generalinfo.dataset_description.License             = 'n/a';
    generalinfo.dataset_description.Authors             = 'JP & KG';
    generalinfo.dataset_description.Acknowledgements    = 'Acknowledgements here';
    generalinfo.dataset_description.Funding             = 'n/a';
    generalinfo.dataset_description.ReferencesAndLinks  = 'n/a';
    generalinfo.dataset_description.DatasetDOI          = 'n/a';

    % general information shared across modality specific json files 
    generalinfo.InstitutionName                         = 'Technische Universitaet zu Berlin';
    generalinfo.InstitutionalDepartmentName             = 'Biological Psychology and Neuroergonomics';
    generalinfo.InstitutionAddress                      = 'Strasse des 17. Juni 135, 10623, Berlin, Germany';
    generalinfo.TaskDescription                         = 'Participants walked repeatedly on a straight path.';
    generalinfo.task                                    = bemobil_config.bids_tasklabel;  


### 3. Specify participant information  

Finally, participant information is **recommended** to be specified also following the description in [BIDS page for modality agnostic files](https://bids-specification.readthedocs.io/en/stable/03-modality-agnostic-files.html).

    % here describe the fields in the participant file
    % for numerical values  : 
    %       subjectData.fields.[insert your field name here].Description    = 'describe what the field contains';
    %       subjectData.fields.[insert your field name here].Unit           = 'write the unit of the quantity';
    % for values with discrete levels :
    %       subjectData.fields.[insert your field name here].Description    = 'describe what the field contains';
    %       subjectData.fields.[insert your field name here].Levels.[insert the name of the first level] = 'describe  what the level means';
    %       subjectData.fields.[insert your field name here].Levels.[insert the name of the Nth level]   = 'describe what the level means';
    %--------------------------------------------------------------------------
    subjectData.fields.age.Description          = 'age of the participant'; 
    subjectData.fields.age.Unit                 = 'years'; 
    subjectData.fields.sex.Description          = 'sex of the participant'; 
    subjectData.fields.sex.Levels.M             = 'male'; 
    subjectData.fields.sex.Levels.F             = 'female'; 
    subjectData.fields.group.Description        = 'experiment group';
    subjectData.fields.group.Levels.young       = 'younger participants under 40';
    subjectData.fields.group.Levels.old         = 'older participants over 65';
    subjectData.fields.handedness.Description    = 'handedness of the participant';
    subjectData.fields.handedness.Levels.R       = 'right-handed';
    subjectData.fields.handedness.Levels.L       = 'left-handed';

This part above is going to be used in writing a Json file for participant data. The role of the Json file is to describe the names and values of the variables that are considered to be "participant information". For generic fields such as age and sex, this kind of description may appear redundant, but the importance becomes obvious when we consider fields that are specific to the experiment, such as "group". Once you have identified and described all variables to be saved as participant information, you can start filling out the table below. 

    % names of the columns - 'nr' column is just the numerical IDs of subjects
    %                         do not change the name of this column 
    subjectData.cols = {'nr',   'age',  'sex',  'group',    'handedness'};
    subjectData.data = {64,     71,     'F',    'old',      'R' ; ...
                       66,     67,     'M',    'old',      'R' ; ...
                       76,     34,     'M',    'young',    'R' ; ...
                       78,     33,     'M',    'young',    'R' };

Note that the first column of this table is reserved for numerical IDs of the participants. It is then followed by the names of variables you have described in **subjecData.fields**. 

Finally, you can now call function **bemobil_xdf2bids.m** with all three of the structs above as inputs. 

        bemobil_xdf2bids(bemobil_config, generalinfo, subjectData)


## Custom functions used for .xdf to BIDS conversion 

### 1. Motion configuration and conversion script 
The two files below (one script and one function) can be replaced by custom mfiles. What they do is to figure out which streams are associated with which object and which channels are associated with which types of motion data (for instance, is that channel containing position x coordinates or y anlges in euler rotations). Motion data collected in BeMoBIL usually come in quaternions, which have to be converted to Euler angles.

- **bemobil_bids_motioncfg.m**   
  default motion specific configuration    
  
      % copy general fields 
      motioncfg       = cfg; 

      % data type and acquisition label
      motioncfg.datatype                                = 'motion';
      motioncfg.acq                                     = 'ViveTrackers';

      % motion specific fields in json
      motioncfg.motion.Manufacturer                     = 'HTC-VIVE';
      motioncfg.motion.ManufacturersModelName           = 'VIVE Pro';
      motioncfg.motion.RecordingType                    = 'continuous'; 

      % coordinate system
      motioncfg.coordsystem.MotionCoordinateSystem      = 'RUF';
      motioncfg.coordsystem.MotionRotationRule          = 'left-hand'; 
      motioncfg.coordsystem.MotionRotationOrder         = 'ZXY'; 

The default configuration file above is used to specify metadata about the motion recording setup and the coordinate system that your motion data uses. If any of the information above deviates from your data set, this file has to be replaced with a custom file. This is best done by copy pasting the content of bemobil_bids_motioncfg.m, modifying the fileds, and renaming it. Then you can specify the name of the custom file in bemobil_config.bids_motioncfg. The script also takes the header information to construct channel information for motion data. The header is first read from the .xdf data and then modified by the 'bemobil_bids_motionconvert.m' below, which can also be replaced with a custom function. 

- **bemobil_bids_motionconvert.m**  
  deafault function for processing the motion data (e.g., quat2eul conversion)
  
      motionOut = bemobil_bids_motionconvert(motionIn, objects, pi, si, di)
 
 This function takes inputs 'motionIn' (motion data directly read from the .xdf file), 'objects' (names of objects in 'bemobil_config.rigidbody_streams'), and then the three indices for participant, session, and datafile. The latter three indices are unused by default. They are there as inputs so one can take care of file-specific problems (e.g., misassigned trackers for some participants) here. The output 'motionOut' will now have correct entries in hdr.label, hdr.chantype, hdr.chanunit that will be used for constructing columns in channels.tsv.
 
     
        quaternionIndices = NaN(1,4); 
    
        for qi = 1:4
            quaternionIndices(qi) = find(strcmp(labelsPre, [objects{ni} '_quat_' quaternionComponents{qi}]));        
        end
    
        cartIndices = NaN(1,3); 
    
        for ci = 1:3
            cartIndices(ci) = find(strcmp(labelsPre, [objects{ni} '_rigid_' cartCoordinates{ci}]));
        end
        
Above is part of the default motionconvert function, demonstrating why it sometimes has to be customized. here the channels containing quaternion components are found by the labels of the channels in .xdf data. It assumes that the label has the format 'objectname_quat_x'. This may of course not be the case for all data sets which is why you might work on that part so that the right channel indices can be found. 
    
  
### 2. Marker processing 
The default function below deals with events. It does a bare minum processing on marker stream included in the .xdf file, namely replacing undefined durations with zeros and generate a generic events.json file. If you are not going to share the BIDS data, it is fine to leave it like this. But if you want to enable people to make sense out of your events and organize differen kinds of information as columns in the events.tsv file, use a custom function for parsing them as you convert data set to bids. 

- **bemobil_bids_parsemarkers.m**  
  deafault function for processing the events 

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
               
and then be merged into one .set file later (when you used bids2set) in that order (*'run-1, run-2, run-11'*), instead of (*'run-1, run-11, run-2'*).

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

## For other modalities than EEG and motion and other things still missing

There is a plan to also support some other modalities that are often used in MoBI research. 
