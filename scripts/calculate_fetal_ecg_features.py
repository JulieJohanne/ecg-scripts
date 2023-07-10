import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from ecgdetectors import Detectors
from compute_12_lead_ecg import calculate_ecg



def get_arguments(input_args):
    parser = argparse.ArgumentParser(description="")
    parser.add_argument(
        "filename", help="Enter path to input file")
    parser.add_argument(
        "filename_ref", help="Enter path to input file")
    return parser.parse_args()


def detect_qrs(ecg_signal, fs=360):
    detectors = Detectors(fs)
    detector = detectors.hamilton_detector
    return np.array(detector(ecg_signal))

def calculate_consecutive_angles(arr):
    v = np.diff(arr)  # check correct axis
    u = v / np.linalg.norm(v, axis=0)
    w = np.zeros(u.shape[1])
    for i in range(len(u)-1):
        w[i] = np.dot(u[:, i], u[:, i+1])
    #w = np.dot(u[:, :-1], u[:, 1:])
    angles = np.arccos(w)
    return angles


def find_perpendicular_angle(angles, tol=1e-3):
    indexes_bool = (angles < (np.pi/2 + tol)) & (angles > (np.pi/2 - tol))
    indexes = np.arange(len(angles))[indexes_bool]
    return indexes[-1]  


def calculate_points(ecg):
    # Detect number of QRS, but not very accurate so neeed to do some refinement
    qrs = detect_qrs(ecg)
    no_qrs = qrs.shape[0]

    # Use find_peaks, but finds too many peaks, so use number of actual QRS found above to sort out the correct ones
    p_R, _ = find_peaks(ecg)
    p_R_index = ecg[p_R].argsort()[-no_qrs:][::-1]
    p_R = p_R[p_R_index]
    p_R = np.sort(p_R)

    # can now find the other features
    p_bQ, p_aS, p_T, p_aT = [], [], [], []

    for i in range(no_qrs):
        # 1) Define search area
        if i == (no_qrs-1):
            cycle = ecg[p_R[i]:]
        else:
            cycle = ecg[p_R[i]:(p_R[i+1]-1)]

        # 1) Find p_aS
        # p_aS.append(np.argsort(np.abs(cycle))[0] + qrs[i])
        peaks, _ = find_peaks(np.gradient(cycle))
        peaks = peaks + p_R[i]
        p_aS.append(peaks[1])

        # 2) Find p_bQ
        # Define search area
        if i == 0:
            search_bQ = ecg[:(p_R[i]-1)]
            angles = calculate_consecutive_angles(
                np.stack((np.arange(len(search_bQ)), search_bQ)))
            p_bQ.append(find_perpendicular_angle(angles))

        else:
            search_bQ = ecg[p_R[i-1]:p_R[i]]
            angles = calculate_consecutive_angles(
                np.stack((np.arange(len(search_bQ)), search_bQ)))
            p_bQ.append(find_perpendicular_angle(angles) + p_R[i-1])

        # 3) Find T-wave
        # Define search area
        search_T = cycle[(p_aS[i]-p_R[i]):(p_bQ[i]-p_R[i])]
        p_T.append(np.argmax(search_T) + (p_aS[i]))

    for i in range(no_qrs):

        # 1) Define search area
        if i == (no_qrs-1):
            cycle = ecg[p_R[i]:]
        else:
            cycle = ecg[p_R[i]:(p_R[i+1]-1)]

        # 4) Find p_bT
        # Define search area
        if i == no_qrs-1:
            search_aT = cycle[(p_T[i]-p_R[i]):]

        else:
            search_aT = cycle[(p_T[i]-p_R[i]):(p_bQ[i+1]-p_R[i])]

        p_aT.append(np.argsort(np.abs(search_aT))[0] + p_T[i])

    return np.stack((p_R, p_bQ, p_aS, p_T, p_aT))  # p_bT,


def calculate_features(points, ecg):
    p_R, p_bQ, p_aS, p_T, p_aT = points

    QT_dur = p_aT - p_bQ
    RR = np.diff(p_R)
    QTc = QT_dur/np.sqrt(np.mean(RR)) * np.sqrt(10**3)
    QRS_dur = p_aS - p_bQ

    return QT_dur, RR, QTc, QRS_dur


def print_features(name, feature, unit, range=[None, None]):
    from colorama import Fore
    from colorama import Style

    if range[0] == None and range[1] == None:
        print(
            f"{name}:\n\t{feature.mean():.2f} +- {feature.std():.2f} {unit}\n")

    elif range[0] == None:
        if feature.mean() < range[1]:
            print(
                f"{name}:\n\t{Fore.GREEN}{feature.mean():.2f}{Style.RESET_ALL} +- {feature.std():.2f} {unit}\n\trange: {range[0]}, {range[1]}")
        else:
            print(
                f"{name}:\n\t{Fore.RED}{feature.mean():.2f}{Style.RESET_ALL} +- {feature.std():.2f} {unit}\n\trange: {range[0]}, {range[1]}")

    elif range[1] == None:
        if feature.mean() > range[1]:
            print(
                f"{name}:\n\t{Fore.GREEN}{feature.mean():.2f}{Style.RESET_ALL} +- {feature.std():.2f} {unit}\n\trange: {range[0]}, {range[1]}")
        else:
            print(
                f"{name}:\n\t{Fore.RED}{feature.mean():.2f}{Style.RESET_ALL} +- {feature.std():.2f} {unit}\n\trange: {range[0]}, {range[1]}")

    else:
        if (feature.mean() > range[0]) and (feature.mean() < range[1]):
            print(
                f"{name}:\n\t{Fore.GREEN}{feature.mean():.2f}{Style.RESET_ALL} +- {feature.std():.2f} {unit}\n\trange: {range[0]}, {range[1]}")
        else:
            print(
                f"{name}:\n\t{Fore.RED}{feature.mean():.2f}{Style.RESET_ALL} +- {feature.std():.2f} {unit}\n\trange: {range[0]}, {range[1]}")


def print_results(points, ecg):
    QT_dur, RR, QTc, QRS_dur = calculate_features(
        points, ecg)

    # values from table 1, week 29-35 in (https://iopscience.iop.org/article/10.1088/1361-6579/ab0a2c/pdf)
    QT_min, QT_max = 232, 295  # or less than half of R-R
    QTc_min, QTc_max = 352.9, 448.1
    QRS_dur_min, QRS_dur_max = 43.2, 55.2
    # amplitude greater than 0.5 mV in at least one standard lead, and greater than 1.0 mV in at least one precordial lead. Upper limit of normal amplitude is 2.5 - 3.0 mV (https://elentra.healthsci.queensu.ca/assets/modules/ECG/normal_ecg.html)

    print_features("QT_dur", QT_dur, "ms", [QT_min, QT_max])
    print_features("QTc", QTc, "ms", [QTc_min, QTc_max])
    print_features("RR", RR, "ms")
    print_features("QRS_dur", QRS_dur, "ms", [QRS_dur_min, QRS_dur_max])


if __name__ == "__main__":
    args = get_arguments(sys.argv)
    filename = args.filename
    filename_ref = args.filename_ref

    ecg = np.loadtxt(filename).T
    ecg_ref = np.loadtxt(filename_ref).T

    ecg = ecg - ecg_ref
    # ecg = ecg - ecg.mean(axis=1)[:, np.newaxis]

    ecg = ecg - ecg.mean()

    points = calculate_points(ecg)

    p_R, p_bQ, p_aS, p_T, p_aT = points

    plt.plot(ecg)

    plt.plot(p_R, ecg[p_R], "*", label="p_R")
    plt.plot(p_aS, ecg[p_aS], "*", label="p_aS")
    plt.plot(p_bQ, ecg[p_bQ], "*", label="p_bQ")
    plt.plot(p_T, ecg[p_T], "*", label="p_T")
    plt.plot(p_aT, ecg[p_aT], "*", label="p_aT")

    plt.legend()
    plt.show()
    plt.close()

    print_results(points, ecg)
