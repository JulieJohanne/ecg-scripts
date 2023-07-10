import sys
import argparse
import numpy as np
from scipy.signal import find_peaks


def get_arguments(input_arguments):
    parser = argparse.ArgumentParser(description="")
    parser.add_argument(
        "paths", help="Enter path to input files")
    return parser.parse_args()

def get_baseline(ecg, t=450):
    baseline = ecg[:, t]
    return baseline

def get_qrs_interval(p_R, interval=50):
    qrs_region_min, qrs_region_max = p_R - interval, p_R + interval
    if qrs_region_min < 0:
        qrs_region_max = qrs_region_max - qrs_region_min
        qrs_region_min = 0
    return qrs_region_min, qrs_region_max


def get_t_interval(p_R, interval=50):
    return p_R + 260 - interval, p_R + 260 + interval  # 270


def calculate_qrs_amplitude(ecg, p_R, interval=50):
    qrs_region_min, qrs_region_max = get_qrs_interval(p_R, interval=interval)

    qrs_amplitude = np.max(ecg[:, qrs_region_min: qrs_region_max], axis=1) - \
        np.min(
        ecg[:, qrs_region_min:qrs_region_max], axis=1)
    qrs_amplitude_indexes = [
        np.argmax(ecg[:, qrs_region_min:qrs_region_max],
                  axis=1) + qrs_region_min,  # +1
        np.argmin(ecg[:, qrs_region_min:qrs_region_max],
                  axis=1) + qrs_region_min  # + 1
    ]
    return qrs_amplitude, qrs_amplitude_indexes


def calculate_t_amplitude(ecg, p_R, baseline=None, interval=50):
    t_region_min, t_region_max = get_t_interval(p_R, interval=interval)

    if baseline is None:
        baseline = get_baseline(ecg)

    t_amplitude, t_index = [], []
    for i, row in enumerate(ecg):
        cycle = row[t_region_min: t_region_max]

        t_amplitude_1 = np.max(cycle) - baseline[i]
        t_amplitude_2 = np.min(cycle) - baseline[i]

        if np.abs(t_amplitude_1) < np.abs(t_amplitude_2):
            t_amplitude_index = np.argmin(cycle) + t_region_min
            t_amplitude.append(np.abs(t_amplitude_2))
            t_index.append(t_amplitude_index+1)
        else:
            t_amplitude_index = np.argmax(cycle) + t_region_min
            t_amplitude.append(t_amplitude_1)
            t_index.append(t_amplitude_index+1)
    return t_amplitude, t_index


def find_p_R(ecg, no_qrs=3):
    r_peak, _ = find_peaks(np.abs(ecg))  # -ecg
    r_peak_index = ecg[r_peak].argsort()[: no_qrs]  # [::-1]
    p_R = r_peak[r_peak_index]
    return sorted(p_R)


def find_cycle(ecg):
    start_index, stop_index = -1, 0
    p_R_index = 0
    p_Rs = find_p_R(ecg[1, :], no_qrs=int(ecg.shape[1]/450))
    p_R = p_Rs[p_R_index]

    while start_index < 0 or stop_index > ecg.shape[1]:
        start_index, stop_index = p_R-70, p_R+430
        if start_index < 0:
            if len(p_Rs) <= p_R_index+1:
                p_R = p_Rs[p_R_index]
                break
            p_R = p_Rs[p_R_index+1]
            p_R_index += 1
            # stop_index = stop_index - start_index
            # start_index = 0

        if stop_index > ecg.shape[1]:
            p_R = p_Rs[p_R_index-1]
            p_R_index -= 1

        start_index, stop_index = p_R-70, p_R+430

    if start_index < 0:
        stop_index = stop_index - start_index
        start_index = 0

    if stop_index > ecg.shape[1]:
        stop_index = ecg.shape[1]
    # p_R_index += 1

    return p_R, p_R_index, start_index, stop_index


if __name__=="__main__":
    args = get_arguments(sys.argv)
    path = args.paths

    ecg = np.loadtxt(path)
    
    p_R, p_R_index, start_index, stop_index = find_cycle(ecg)
    qrs_region_min, qrs_region_max = get_qrs_interval(p_R, interval=50)
    t_region_min, t_region_max = get_t_interval(p_R, interval=50)
    cycle = ecg[:, start_index:stop_index]

    baseline = get_baseline(cycle)
    cycle = cycle - baseline[:, np.newaxis]
    cycle = cycle * 1e3

    qrs_amplitude, qrs_index = calculate_qrs_amplitude(cycle, p_R - start_index)
    t_amplitude, t_index = calculate_t_amplitude(cycle, p_R - start_index)

    tqrs_ratio = t_amplitude/qrs_amplitude

    print(qrs_amplitude)
    print(t_amplitude)
    print(tqrs_ratio)
