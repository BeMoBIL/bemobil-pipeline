The BeMoBIL Pipeline is an open-source MATLAB toolbox for fully automatic, transparent, and replicable processing and visualization of mobile brain/body imaging and other EEG data. It includes wrappers for EEGLAB functions, the use of additional EEGLAB plugins, and comes with additional new functionalities. 

All parameters are configurable in a central script and everything is stored in the EEG.etc struct. Additionally, analytics plots are generated for each step.

***

A comprehensive guide to installing, using, and understanding the pipeline can be found [in our wiki on github](https://github.com/BeMoBIL/bemobil-pipeline/wiki)!

***

For a quick start, we recommend you to have a look at the following scripts in the root directory of the pipeline. 

    example_bemobil_process_all_EEG_data.m
    example_bemobil_process_all_motion_data.m
    example_bemobil_config.m
    example_bids_metadata.m


The two scripts `example_bemobil_process_all_EEG_data.m` and `example_bemobil_process_all_motion_data.m` will both load `example_bemobil_config.m`. So the configuration file serves as the summary of settings used throughout the pipeline. Comments within the file explain the parameters. These scripts are the only scripts that need to be run. They contain all steps from the raw to the processed and cleaned data and allow batch processing of all subjects.

The `example_bids_metadata.m` script contains the metadata struct that one can enter as input to `bemobil_xdf2bids.m`. If you already have your data in BIDS or do not need to add metadata during xdf to bids conversion, you may skip this step by commenting out corresponding lines in `example_bemobil_process_all_EEG_data.m` which contains codes that import data. Specific instruction is given as comments.

Here is an example of the final BeMoBIL pipeline output :   

![folder structure of the pipeline output](https://raw.githubusercontent.com/BeMoBIL/bemobil-pipeline/master/wiki_images/mainwiki/output-folders.png)

***

If you have any comments, questions, or suggestions, please [open issues on git](https://github.com/BeMoBIL/bemobil-pipeline/issues), join our [discord server](https://discord.gg/xJMru7XVXY), or write an email to marius.s.klug@gmail.com!
