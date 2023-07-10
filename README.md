# ECG scripts

## calculate_12_lead_ecg.py
Calculates the 12 leads of 12 lead ECG 
Input: phie at 10 electrode locations for 12 lead ECG. 
For the maternal 12 lead ECG, should be node numbers (but double check perhaps in meshalyzer that locations look consistent with 12 lead ECG locations): 
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

calculate_fetal_ecg_features.py
Calculates, prints and plots QT_dur, RR, QTc, QRS_dur  for fetal ECG.
Input: 1. Phie extracted at fetal left leg (node number 1395799) . 2. Phie extracted at fetal  right arm (node number 1264521). 
(To reconstruct lead 2 for fetus)

calculate_maternal_ecg_features.py
Calculates, prints and plots QT_dur, RR, QTc, QRS_dur for maternal ECG.
Input: phie at 10 electrode locations for 12 lead ECG as input, then the script extracts lead 2. 

calculate_ecg_amplitudes.py
Calculates T wave and QRS amplitude
Input: ECG trace (so for instance one trace like lead 2 or entire 12 lead ecg should also work) 