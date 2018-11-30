
# coding: utf-8

# In[1]:


import imctools.scripts.pipeline.configuration as configuration
import imctools.scripts.pipeline.steps as steps


# In[2]:


import os
import re
import zipfile
import sys


# In[3]:


from ruamel.yaml import YAML


# In[4]:


# This should be removed in the final version
get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '1')
get_ipython().run_line_magic('aimport', 'imctools.scripts.pipeline.steps')
get_ipython().run_line_magic('aimport', 'imctools.scripts.pipeline.configuration')
get_ipython().run_line_magic('aimport', '')


# 
# # The IMC preprocessing pipeline for multiplexed image analysis
# 

# This is a pipeline to segment IMC data using Ilastik pixel classification as well as CellProfiler.
# 
# It requires:
# - CellProfiler 3.1.5: http://cellprofiler.org/releases/
# - Ilastik: http://ilastik.org/
# - The following Github repositories:
#     => Either Clone (e.g. using command line git or the github deskop client: https://desktop.github.com/) OR download and unzip these repositories to a local folder:
#     - ImcPluginsCP plugins of the development branch: develop-cp3, https://github.com/BodenmillerGroup/ImcPluginsCP/tree/develop-cp3
#           - This repository contains additional CellProfiller modules.
#           - **In the preferences of CellProfiller point the CellProfiler Plugins folder to `ImcPluginsCP/plugins`**
#       
#     - The ImcSegmentationPipeline repository: https://github.com/BodenmillerGroup/
#         - contains the conda environment that you will use bellow in the `setup` folder
#         - contains the CellProfiler pipelines to be used
# 
# 
# - Install the `imctools` conda environment
#     - If you have already an older version of this installed, please uninstall it first
#  
#     - To install the `imctools` envrionment:
#         -> Install conda: https://www.anaconda.com/download/#linux
#         -> The conda environment file is found in  `ImcSegmentationPipeline/Setup/conda_imctools.yml`
#         -> On a conda console type: `conda env create -f PATHTO/conda_imctools.yml` OR use the Anaconda GUI -> Environments -> Import -> choose `setup/conda_imctools.yml`
#         -> Start a Jupyter notebook instance *in this conda environment* and open this script: 
#             -`conda activate imctools`
#             -`conda jupyter notebook`
#             - OR in the GUI: choose the `imctools` environment, start `Jupyter Notebook`
#             - open the `ImcSegmentationPipeline/scripts/imc_preprocessing.ipynb`
#             - Execute the script cell by cell using `shift-enter`
# 
# - If compensation should be done (https://www.cell.com/cell-systems/abstract/S2405-4712(18)30063-2), the following additional requirements are needed:
#     - A spillover matrix specific to the isotope lots used for antibody conjugation
#         - Use the experimental protocol of:https://docs.google.com/document/d/195eViUqHoYRKrkoy_NkIdJPmyx1-OuDaSjiWQBy4weA/edit
#             -> Will result in a spillovermatrix in `.csv` format.
#         - Clone the Github repository: 
#         - R > 3.5 (https://cran.r-project.org/bin/windows/base/)
#         - Rstudio: https://www.rstudio.com/products/rstudio/download/
#         - CATALYST >= 1.4.2: https://bioconductor.org/packages/release/bioc/html/CATALYST.html
#         - 'tiff' R library: run `install.packages('tiff')`
# 
# 
# - Data requirements:
#     - This scripts assume that *each `.mcd`  acquisition and all `.txt` files corresponding to this '.mcd' acquisition* are saved in one a seperate `.zip` folder.
#         -> This is my recomended data format as it preserves and contains all original metadata and enforces a consistent naming scheme.
#     - see the example files that are downloaded bellow for an example
# 
# Note that the `description` image name can be found in the `..._Acquisition_meta.csv` generated together with the ome tiffs as well as in the `cpout` folder later in the script.
# After analysis the `Image.csv` metadata file generated in Cellprofiller will also contain the `Description` as well as other important metadata for each 
# image, such as acquisition frequency, time, location etc.
# 
# For working with `.txt` files only, please look at the older examples.
# 
# For any feedback please contact: Vito, vito.zanotelli@uzh.ch or even better raise an issue on this Github page!

# ### Input folders (Needs to be adapted for use)

# In[5]:


fn_config = '../config/config_pipeline.yml'


# In[6]:


get_ipython().run_cell_magic('file', '$fn_config', "# %load configuration.conf_str\n\n# This section specifies the input files.\ninput_files:\n    # This defines a list of folders  where to look for input files which must be zip files\n    # containing each 1 .mcd file together with (optionally)  all the .txt files of this\n    # acquisition session.\n    folders: ['../example_data']\n    # This regular expression can be used to further select .zip files for processing\n    file_regexp: '.zip'\n\n# This section defines the output folders\noutput_folders:\n    # This indicates the base folder for all the output subfolders:\n    base: '/home/vitoz/Data/Analysis/201811_icp_segmentation_example4'\n    subfolders:\n        # The ome subfolder contains all the acquisitions as ometiffs together with the\n        # metadata in .csv format.\n        ome: 'ometiff'\n        # The cpout subfolder will contain the final cellprofiler output.\n        cpout: 'cpout'\n        # The analysis subfolder will contain all the .zip files that are used for analysis.\n        analysis: 'analysis'\n        # The ilastik folder will contain all the .hd5 files used for Ilastik Pixel classification.\n        ilastik: 'ilastik'\n        # The uncertainty folder can contain the uncertainties that can be generated from the\n        # probabiprobabilities from the Ilastik pixel classification.\n        uncertainty: 'uncertainty'\n        # The Histocat folder will contain the images in the HistoCat folder structure.\n        histocat: 'histocat'\n\n# This configuration saves different substacks of the acquired images in a format compatible for\n# analanalysis with CellProfiler and Ilastik.\nanalysis_stack_generation:\n    # The pannel is a CSV file with a column representing the metal tag and a column containing\n    # boolean (0/1) values which channels to select for the different stacks.\n    pannel:\n        # The filename of the pannel .csv file\n        csv_filename: '../config/example_pannel.csv'\n        # The name of the column containing the metal tag\n        metal_column: 'Metal Tag'\n    # An arbitrary number of stacks can be defined. Classically a 'full' stack is defined,\n    # containing all the channels(=metals) that should be measured using CellProfiler as well as\n    # an 'ilastik' stack, containing only a subset of channels that contain information for the\n    # pixelpixel classification used for segmentation.\n    stacks:\n        # The stack name indicated the name of the stack to be generated\n        - stack_name: 'full'\n          # The stack column is the column that contains the boolean selection for the channels to\n          # be used for this stack\n          csv_column: 'full'\n          # The suffix is added to the filename before the file ending and is used to identify the\n          # stack in the later analysis.\n          suffix: '_full'\n\n        - stack_name: 'ilastik'\n          csv_column: 'ilastik'\n          suffix: '_ilastik'\n          # With aditional stack parameters, keyword arguments used by the ometiff_2_analysis script\n          additional_parameters:\n            # A commonly used one parameter to set is 'addsum', which will add the sum of all channels of the stack as\n            # first channel of the stack, which is convenient for forground/background\n            # discrimination in pixel classification.\n            addsum: True")


# In[7]:


conf = configuration.load_config(fn_config)


# In[8]:


steps.initialize_folderstructure(config=conf)


# ### Convert zipped IMC acquisitions to input format
# 
# This script works with zipped IMC acquisitions:
# Each acquisition session = (1 mcd file) should be zipped in a folder containing:
# - The `.mcd` file
# - All associated `.txt` file generated during the acquisition of this `.mcd` file -> Don't change any of the filenames!!

# In[9]:


steps.convert_rawzip_to_ome(conf)


# Generate a csv with all the acquisition metadata

# In[10]:


steps.exportacquisitionmetadata(conf)


# Generate the analysis stacks

# In[11]:


steps.generate_analysisstacks(conf)


# Convert ome.tiffs to a HistoCAT compatible format, e.g. to do some visualization and channel checking.
%%time
if not(os.path.exists(folder_histocat)):
    os.makedirs(folder_histocat)
for fol in os.listdir(folder_ome):
    ome2micat.omefolder2micatfolder(os.path.join(folder_ome,fol), folder_histocat, dtype='uint16')

# # Next steps
# 
# This concludes the conversion of the IMC rawdata into usable TIFFs.
# 
# The pipelines can be found in the `cp3_pipeline` folder in this repository. They were tested in `cellprofiler 3.1.5`.
# 
# The next steps are:
# 
# ### A) Cellprofiler: 1_prepare_ilastik
# 
# In this module we prepare the data for Ilastik pixel classification, by first removing strong outlier pixels, then scaling the images 2x and then taking random 500x500 crops to do the train the pixel classifier.
# 
# The following parts of this module need to be adapted:
# 
# 1) File list: choose all files in the `tiff` subfolder
# 
# 2) Default Output Folder: Choose the `ilastik` subfolder
# 
# No further parts need to be adapted.
# In our 16 core computer this step takes ca 5 min for the example dataset.
# 
# 
# ### B) Ilatik: Train a pixel classifier
# 
# This uses the random crops generated in the last step.
# 
# 1) Make a new `pixel classification project`. Save the project file in the `ilastik` subfolder.
# 
# 2) Add the `.h5` random crops: Raw data -> Add Seperate Images -> Select all `.h5` images in the `ilastik` subfolder.
# 
# 3) Proceed to `Feature Selection`
# 
# 4) Select suitable features (or just everythin > 1 pixels)
# 
# 5) Proceed to the classification:
#     - Add 3 labels:
#         - 1: Nuclei
#         - 2: Cytoplasma/membrane
#         - 3: Background
#         - -> For large datasets adding the labels can take a while
#     - Start labeling: 
#         - The box next to `Input Data` can change the channels. What each channel corresponds to can be seen when looking in any of the `..._ilastik.csv` files in the `tiff` folder. The 0 channel correspond to the sum of all channels, very usefull to label the background.
#         - Use window leveling change the contrast. Right click on the `Input Data` -> `Adjust Thresholds` is also very usefull
#         - Label opiniated: If you see in the nucleus channel that two nuclei are stuck together but have a faint dip in intensity in between, label this as 2: Cytoplasma. Encyrcle nuclei with Cytoplasma
#         - Diseable `Live Update` for performance
#         - Frequently check the `Uncertainties`: This indicates which pixels the classifier profits most if they are labeled. A well trained classifier has low uncertainty within class regions (e.g. Nuclei) and high uncertainty at class borders (e.g. between nuclei and cytoplasma).
#         
# 6) If you think the classifier is well trained, export the probabilities:
#         - Export Settings -> Source: Probabilities -> Choose Export Image Settings:
#             - Convert to datatype: Unsigned Integer 16 bit
#             - Renormalize: check
#             - Format: Tiff
#             - File: leave default
#         - Export all: This generates `_Probabilities.tiff` in the `ilastik` folder. They can be checked using any image viewer
#             - To generate uncertainty maps (good to identify regions that need training),
#             run the `Convert probabilities to uncertainties` section `#For training` below. This will put uncertainties in the uncertainty folder.
#             -> Well trained classifiers have low uncertainty (transparent) everywhere but at class borders which should be white.
#             
#         - Optional: Train again regions with high uncertainty, then proceed.
#         
# 7) If you think that you are finished with classification, you need to apply the classifier to your whole dataset usingBatch processing:
#         - Make sure the Export Settings are still the same as in the step 6), then go to the 'Batch Processing' step in ilastik.
#         - Select raw data files -> select all `_s2.h5` files in the `tiff` folder. (sort by filetype, select all `H5` files).
#         => This step takes a while and is computationally intensive!
#         => Ca 15 min on 10 cores on the example data
#             
#         - Optional: use the below probability to uncertainty `#For the data` to convert all proabilities to uncertainties, check if there are any regions of high uncertainty and optionally crop the corresponding image part in imagej and add it to the training data.
#         - Note: store the `ilastik` folder with all the random crops and the trained classifier for reproducibility reasons.
#         
# ### C) Cellprofiler: 2_segment_ilastik
# 
# This step will segment the probabilities into masks.
# 
# Things to adapt:
# 
# 1) File list: choose again all files from the `tiffs` folder
# 
# 2) It is important to check the `IdentifyPrimaryObjects` step, if the segmentation settings are suitable!
#     This might vary strongly between cell/tissue/training and needs attention! Use the test mode and try various settings.
#     Also note the `smooth` step immediately before: This can be also removed, I just happen get good results with this additional step.
#     
# 3) Also the `MeasureObjectSizeShape` combined with `FilterObjects` is just some personal preference of mine, feel free to change
# 
# 4) `IdentifySecondaryObjects`: Here th mask is expanded to the full cell.
# 
# 5) `Rescale objects`: note that our segmentation was done on 2x upscaled images, this scales the masks down again. Note that potentially also the nuclei mask could be scaled down and further exported and used.
# 
# 6) The `Default Output Path` does not need to be adapted for this module.
# 
# 
# Note1: Seperating mask generation from mask measurement adds modularity and is thus highly recommended, as generating masks is one of the most resource intensive steps.
# 
# 
# ### D) Cellprofiler: 3_measure_mask
# 
# This step is not necessary for `HistoCat` only analysis. If `HistoCat` should be used, use the `Generate the histocat folder with masks` section below.
# 
# #### 3_measure_mask_basic
# 
# This module measures without considering spillover correction.
# 
# 1) File list: choose again all files from the `tiffs` folder
# 
# 2) View Output settings: set the `Default output folder` to the `cpout` folder
# 
# 3) Metadata: update - this will automatically merge the mcd metadata .csv generated earlier in the script with your images.
# 
# 4) Names and types: click update
# 
# 5) `Measure Object Intensity Multichannel`: Adapt the channel numbers. Check the `_full.csv` files in the `tiffs` folder to see how many channels the stack have and adapt accordingly.
# 
# 6) `Measure Image Intensity Multichannel`: Adapt the channel numbers. Check the `_full.csv` files in the `tiffs` folder to see how many channels the stack have and adapt accordingly.
# 
# Notes:
# - In this pipeline all the intesities are scaled by `1/(2**16)`
# - The mapping between channel number c1, c2, c3 corresponds to the position in the `_full.csv`s found in the `tiffs` folder.
# - The original acquisition description, acquisition frequencies etc can be found in the `Image.csv` output as `Metdata_...` columns.
# - This outputs a lot of measurements that are acutally of little interest - usually we only look at `meanintensity` per channel and cell.
#     To reduce the outputs, select in `Export To Spreadsheet` -> `Select Measurements to Export` -> Only the measurements you want (usually all Image measurements and only the `MeanIntensity` fullstack measurements).
# - The `FullStack` can also be not measured, as it is almost identical to the `FullStackFiltered`.
# 
# #### 3_measure_mask_compensation
# This will do measurements and also single cell data compensation
# 0) Run the script: https://github.com/BodenmillerGroup/cyTOFcompensation/blob/master/scripts/imc_adaptsm.Rmd in R:
#     - Adapt the path to the spillover matrix `fn_sm='.../path/to/sm/spillmat.csv'`. In this example data it can be found at:
#         `fn_sm = 'PATHTO/ImcSegmentationPipeline/config/20170707_example_spillmat.csv'`
#     - Choose any `_full.csv` file, generated during the `Generate analysis stacks` step in the output folder, for the `fn_imc_metals = '/path/to/anyfile_full.csv' `.
#         In this example this could be: `fn_imc_metals = 'PATHTO/tiffs/20170905_Fluidigmworkshopfinal_SEAJa_s0_p0_r0_a0_ac_full.csv'`
#     - Run the script and this will produce an `PATHTO/tiffs/imc_full_sm.tiff` file
# 
# 1) File list: choose again all files from the `tiffs` folder
# 
# 2) View Output settings: set the `Default output folder` to the `cpout` folder
# 
# 3) Metadata: update - this will automatically merge the mcd metadata .csv generated earlier in the script with your images.
# 
# 4) Names and types:  Make sure that in `NamesAndTypes` the `PATHTO/tiffs/imc_full_sm.tiff` file is selected, click update
# 
# 5) `Measure Object Intensity Multichannel`: Adapt the channel numbers. Check the `_full.csv` files in the `tiffs` folder to see how many channels the stack have and adapt accordingly.
# 
# 6) `Measure Image Intensity Multichannel`: Adapt the channel numbers. Check the `_full.csv` files in the `tiffs` folder to see how many channels the stack have and adapt accordingly.
# 
# 7) `CorrectSpilloverApply`: This will generate a corrected image stack, this can be used e.g. to do measurements of intensity distribution. For measurements of intensity it is however better to correct the measurement afterward using the `CorrectSpilloverMeasurement`.
# 
# 8) `CorrectSpilloverMeasurement`: Here the intensity measurement can be spillover corrected. Note that this makes only sense for linear combinations of intensity measurements such as `MeanIntensity` or `TotalIntensity`. For these it is more accurate to do this after measurement than doing it on the pixel level beforehand. Note that for things with non linear transformations as `MedianIntensity`, this will not result in valid results and these measurements should be done on beforehand corrected images from `CorrectSpilloverApply`.

# ## Convert probabilities to uncertainties

# In[ ]:


# For training
for fn in os.listdir(folder_ilastik):
    if fn.endswith(suffix_probablities+'.tiff'):
        print(fn)
        probablity2uncertainty.probability2uncertainty(os.path.join(folder_ilastik,fn), folder_uncertainty)


# In[ ]:


# For the data
for fn in os.listdir(folder_analysis):
    if fn.endswith(suffix_probablities+'.tiff'):
        print(fn)
        probablity2uncertainty.probability2uncertainty(os.path.join(folder_analysis,fn), folder_uncertainty)


# ## Generate the histocat folder with masks

# In[ ]:


get_ipython().run_cell_magic('time', '', "for fol in os.listdir(folder_ome):\n    ome2micat.omefolder2micatfolder(os.path.join(folder_ome,fol), folder_histocat, \n                                         fol_masks=folder_analysis, mask_suffix=suffix_mask, dtype='uint16')")

