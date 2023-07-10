import numpy as np


def calculate_ecg(ecg_raw):
    la, ra, ll, rl, v1, v2, v3, v4, v5, v6 = ecg_raw
    # la, ra, ll, v1, v2, v3, v4, v5, v6 = la - rl, ra - rl, ll - \
    #    rl, v1 - rl, v2 - rl, v3 - rl, v4 - rl, v5 - rl, v6 - rl

    v_w = (ra + la + ll)/3

    l_1 = la - ra
    l_2 = ll - ra
    l_3 = ll - la

    a_vr = (3/2)*(ra-v_w)
    a_vl = (3/2)*(la-v_w)
    a_vf = (3/2)*(ll-v_w)

    v1, v2, v3, v4, v5, v6 = v1-v_w, v2-v_w, v3-v_w, v4-v_w, v5-v_w, v6-v_w

    return np.stack((l_1, l_2, l_3, a_vr, a_vl, a_vf,
                     v1, v2, v3, v4, v5, v6))


if __name__ == "__main__":
    pass
