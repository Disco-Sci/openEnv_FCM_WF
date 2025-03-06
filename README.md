# openEnv_FCM_WF
Development space for flow cytometry (FCM) data (.fcs) processing workflows with open source software

05MAR2025

An introduction

The R script "flowCore_wf_mvp_script.R" contains many wrapper functions, yet possibly not enough, that are called so that all FCS files within a directory and it's sub-directories are processed. 

My primary goal for this workflow procedure is to auto-gate CD3+/CD4/CD8 populations using only the flowCore package. This includes compensation using the defined spill matrix from the device and transforming using the "logicle" method (AKA hyperbolic sine transformation). 

The method I employ is using k-means clustering which I think works well enough with distinct signaling populations. With this method I am able to automatically select my lymphocyte population based on size (FSC) then select the CD3 positive population and finally enumerate the number of CD4-/CD8-, CD4+/CD8-, CD4-/CD8+, and CD4+/CD8+ stained objects. How this is performed is detailed in the Process section and I will do my best to describe the mapping there. 

The secondary goals are to organize and structure the data in an automated fashion. Device inputs are stored in FCS files which flowCore can access. Therefore, if the setup of the flow cytometry data acquisition at the bench is defined at some level, then at the backend the data can be arranged and structured to data analysis. 

I demonstrate here in this R file a procedure that produces a list of data.frames of CD3 positive CD4 and CD8 populations. They are grouped by experiment names that are defined in the FCS file. Further, I generate basic descriptive statistics such as the five-number summary and the % populations. 

My goal is continue developing my workflows and implement better gating solutions as well as provide graphical outputs at some point in time. Hopefully I continue to find fulfillment in this work and so can advance my understading of computational methods and computer science as a whole. 

Process: 

A listing of files with the .FCS extension is generated by searching through the assigned directory path. Each file is then read in as a flowCore::flowFrame. 

Specific device data written into the FCS file, which I refer to as the metadata, is assigned into a list. This data identifies specific device configuration, device serial, experimental inputs, etc. that are relevant to understanding the flow cytometry experimentation and some quality control interests. Importantly, the marker and channel can be found in a key-value pairing (named vector) that I did not employ, but regret not doing as I am writing this. 

Despite that simple feature missing, the next steps are compensation and transformation. Compensation works by taking the spill matrix that is already part of the FCS file. Defining your own spill matrix can be added; however, this was not needed for the data set I used to develop this procedures. The transformation likewise does not need to be defined as the functions are readily available, but if a bespoke transformation is needed then the flowCore allows you to define that as an object and apply to your data. Unfortunately the transformation is not always successful and generates errors. This is wrapped by a Try call and my workaround is to create NA values. Out of 453 files, I only generated 10 or so errors. 

Next I use the pre-determined channel parameters that I happen to know are consistent throughout my file set. Though I checked programmatically, like I mentioned before I could have coded that into my workflow. 

Individually the FSC-A, Am-Cyan-A, APC-Cy7-A, and Pacific-Blue-A measures are fed into the kmeans clustering function with defined expected groupings. The unwanted objects are gated  out from the data until the I have only the CD3+ cells and I then count the number of CD4 and CD8 positive and negative counts.

The counts are then placed into a data frame along with metadata identifiers such as the file name, experimental group and so on. 

That dataframe is the further processed to group the data by experiment which is structured into a list. 

The next set of functions are applied to that list to statistical values. 

Performance: 

I process 453 FCS files. This takes about 7 minutes to complete on my computer.   
