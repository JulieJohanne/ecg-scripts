# ECG scripts

Use igbextract to extract phie at specific nodes. All phie files should be ascii/text files, not binary files.

## calculate_12_lead_ecg.py
* Calculates the 12 leads of 12 lead ECG 
* Input: phie at 10 electrode locations for 12 lead ECG. 
For the maternal 12 lead ECG, should be node numbers: 
125 
4297 
1505
4725
565
996083
1075031
6113
6817
4198

## calculate_fetal_ecg_features.py
* Calculates, prints and plots QT_dur, RR, QTc, QRS_dur  for fetal ECG.
* Input: (To reconstruct lead 2 for fetus)

    File 1: Phie extracted at fetal left leg

    File 2: Phie extracted at fetal  right arm


## calculate_maternal_ecg_features.py
* Calculates, prints and plots QT_dur, RR, QTc, QRS_dur for maternal ECG.
* Input: phie at 10 electrode locations for 12 lead ECG as input, then the script extracts lead 2. 

## calculate_ecg_amplitudes.py
* Calculates T wave and QRS amplitude
* Input: ECG trace (so for instance one trace like lead 2 or entire 12 lead ecg should also work) 