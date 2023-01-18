# The BeMoBIL Pipeline

The BeMoBIL Pipeline is an open-source MATLAB toolbox for fully synchronized, automatic, transparent, and replicable import, processing and visualization of mobile brain/body imaging and other EEG data. It includes wrappers for EEGLAB functions, the use of various EEGLAB plugins, and comes with additional new functionalities. 

All parameters are configurable in central scripts and everything is stored in the EEG.etc struct. Additionally, analytics plots are generated for each step.

A comprehensive guide to using and understanding the pipeline can be found [in our wiki on github](https://github.com/BeMoBIL/bemobil-pipeline/wiki)!

The BeMoBIL Pipeline can be cited as follows: 

Klug, M., Jeung, S., Wunderlich, A., Gehrke, L., Protzak, J., Djebbara, Z., Argubi-Wollesen, A., Wollesen, B., & Gramann, K. (2022). The BeMoBIL Pipeline for automated analyses of multimodal mobile brain and body imaging data. In bioRxiv (p. 2022.09.29.510051). https://doi.org/10.1101/2022.09.29.510051

Latest release: [![DOI](https://zenodo.org/badge/141275506.svg)](https://zenodo.org/badge/latestdoi/141275506)

# Quick start

For a quick start, install the pipeline and all its requirements as laid out [here](https://github.com/BeMoBIL/bemobil-pipeline/wiki/Installation-and-Requirements). Also download the example data and scripts [from OSF](https://doi.org/10.17605/OSF.IO/8UPNW). The example data repository contains a scripts folder, in which there is a script `run_all_examples.m`, which runs through the entire pipeline. This also calls the following examples within the `examples` folder of the pipeline itself:

    example_bemobil_config_script.m
    example_bemobil_import.m
    example_EEG_pipeline.m
    example_MOTION_pipeline.m

You will have to edit the filepath in the `example_bemobil_import` and `example_bemobil_config_script` scripts and then remove the error message. The import script is complete on its own. The scripts `example_EEG_pipeline.m` and `example_MOTION_pipeline.m` will all load `example_bemobil_config_script.m`. So the configuration file serves as the summary of settings used throughout the processing. Comments within the files explain the parameters. These scripts are the only scripts that need to be run. They contain all steps from the source xdf data over the raw imported data to the processed and cleaned data, and allow batch processing of all subjects.

The `example_bemobil_import.m` script contains an exemplary import process from xdf over BIDS to EEGLAB structure. If you already have your data in EEGLAB set files you may skip this step entirely, if you have your data already in BIDS, you can just run the bottom part that loads from BIDS to EEGLAB set. Specific instruction is given as comments.

Here is an example of the final BeMoBIL pipeline output :   

![folder structure of the pipeline output](https://raw.githubusercontent.com/BeMoBIL/bemobil-pipeline/master/wiki_images/mainwiki/output-folders.png)

In the single subject EEG analysis folder there are two final data sets after the complete processing is done (`xxx_preprocessed_and_ICA.set` and `xxx_cleaned_with_ICA.set`). Both sets have basic preprocessing done (line noise removal, channel locations, channel interpolation, removal of very slow trends), and contain ICA information. The `xxx_cleaned_with_ICA.set` file additionally has ICs removed as determined by the settings for ICLabel in the pipeline config. If ICA was only meant to be used for cleaning, any kind of sensor-level analysis (like ERPs at Pz electrode) can now be performed on the cleaned data. If the end goal is source analysis, and potentially analysis of muscle or eye activity in conjunction with other MoBI modalities, this data is still available in `xxx_preprocessed_and_ICA.set`. Consider [our repeated clustering approach](https://github.com/BeMoBIL/bemobil-pipeline/wiki/repeated-clustering) in that case.

***

The pipeline has more functionalities than this, so don't forget to check out [our wiki on github](https://github.com/BeMoBIL/bemobil-pipeline/wiki) with details on the parameters and usage, as well as more info on event creation and other additional options the pipeline has to offer!

***
If you have any comments, questions, or suggestions, please [open issues on git](https://github.com/BeMoBIL/bemobil-pipeline/issues), join our [discord server](https://discord.gg/xJMru7XVXY), or write an email to marius.s.klug@gmail.com!
