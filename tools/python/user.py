import os

if os.environ['USER'] in ['dspitzbart', 'dspitzba']:
    results_directory       = "/afs/hephy.at/data/dspitzbart01/nanoMET/results/"
    plot_directory          = "/afs/hephy.at/user/d/dspitzbart/www/nanoMET/"
    dpm_directory           = '/dpm/oeaw.ac.at/home/cms/store/user/dspitzba/'
    data_directory          = "/afs/hephy.at/data/dspitzbart03/nanoSamples/"
    postprocessing_output_directory = "/afs/hephy.at/data/dspitzbart03/nanoSamples/"
    analysis_results        = results_directory
    
if os.environ['USER'] in ['llechner']:
    results_directory       = "/afs/hephy.at/data/llechner03/nanoMET/results/"
    plot_directory          = "/afs/hephy.at/user/l/llechner/www/nanoMET/"
    dpm_directory           = '/dpm/oeaw.ac.at/home/cms/store/user/llechner/'
    data_directory          = "/afs/hephy.at/data/llechner03/nanoSamples/"
    postprocessing_output_directory = "/afs/hephy.at/data/llechner03/nanoSamples/"
    analysis_results        = results_directory
    
