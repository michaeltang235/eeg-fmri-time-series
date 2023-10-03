# eeg-fmri-time-series
This repository highlights some of the scripts I wrote for the project of using fMRI time series features to detect brain areas
generating seizure-related acitivity. Below is a brief summary of the project. 

**Background:**\
The current clinical gold standard to identify epileptic brain areas involves implanting electrodes in the brain 
during presurgical workup, yet only offers limited brain coverage with risks of developing a host of complications.
Resting-state fMRI (rs-fMRI), on the other hand, measures spontaneous activity at high spatial resolution 
across the entire brain with no invasive procedures required during image acquisition. Previous rs-fMRI studies employed
time series features including regional homogeneity (ReHo) and amplitude of low-frequency fluctuations (ALFF) to measure
epileptic activity in local brain areas, and function connectivity (FC) to measure correlation in rs-fMRI signal 
fluctuations between distinct areas of interest.

Despite indications that the time series features mentioned above could provide information about the locations of 
epileptic areas, te number of studies directly compared the intensity of seizure activity with magnitude of time series features 
for higher accuracy was rather limited. Moreover, the dynamic nature of brain acitivty, which fluctuates over time 
scales ranging from seconds to minutes, was not taken into account in majority of the studies published thus far. Hence, 
this project investigated static and dynamic features (ReHo, ALFF, fractional fALFF, and FC) of rs-fMRI time series 
in focal epilepsy patients with simultaneous intracranial EEG recordings to confirm the pathological relevance 
of the observed rs-fMRI abnormalities. Specifically, the influence of interictal epileptic discharges (IEDs) 
on static and dynamic rs-fMRI properties was determined.

**Analyses:**\
In static analysis, rs-fMRI features (ReHo, ALFF, and fractional ALFF, FC) were computed over the entire duration
of fMRI scan and were correlated to the rate of IEDs occurring at the same locations.

*related scripts: filenames with 'static_' included

In dynamic analysis, the trasient nature of epileptic activity, which fluctuated over the duration of fMRI scan, was considered.
Accordingly, a sliding-window approach was implemented in which time series was decomposed into segments of
120 s with 50% overlap. In every window, rs-fMRI features were computed and compared with the instantaneous IED rate.

*related scripts: filenames with 'dynamic_' included

**Folder structure:**\
Under 'subject_code', there are scripts for subject-specific analyses, such as preprocessing of fMRI images, 
matching IED onset times to the start of fMRI scans. <br>

Under 'code/', there are scripts for analyzing EEG-fMRI data across all subjects in database all at once. <br>

Under 'plots_code/', there are scripts that were used to make plots. <br>

**Conclusions**\
This work provides evidence that rs-fMRI time series features could locate epileptic brain areas for presurgical evaluation, 
guiding (offering complementary information for) physicians on implantation of intracranial electrodes and resection surgery. 
Among the features studied, the magnitude of ReHo, fALFF and FC exhibited robust correlation with IED rates in 
both static and dynamic analyses. Regions with high IED rates exhibited reduced ReHo and fALFF, and high FC 
when compared to regions with no IEDs registered. ALFF, on the other hand, displayed significant correlation with 
IED rates only in dynamic analysis, where the magnitude of ALFF increased with temporal fluctuations in 
instantaneous IED rate, illustrating that additional connectivity would have been left unnoticed if temporal 
variations in rsfMRI features and IED rates were neglected. 
Overall, with no information about epileptogenicity in any voxel of the brain, areas that are epileptic could 
be estimated based solely on the magnitude of ReHo and fALFF, and on temporal variations in ALFF. 
With initial knowledge of some voxels having IEDs recorded, other voxels that are likely to possess the same 
type of IEDs could be identified with FC.
