"""
This function is to create resolution for vegetation or fences or a porous
area for SWASH model. The resolution will be used for vegetation
implementation
"""


def create_xp_zp_veg_sws(dx, xp_sws, zp_sws,
                         total_length, loc_veg, width_veg):
    # ------------------------------------------------------------------------
    # Package
    import numpy as np
    import matplotlib.pyplot as plt
    # ------------------------------------------------------------------------
    # Calculations
    length_0veg_begin = loc_veg - xp_sws[0]
    length_0veg_rest = total_length - length_0veg_begin - width_veg
    numb_rep_length_0veg_begin = np.round(length_0veg_begin / dx, 0)
    numb_rep_width_veg = np.round(width_veg / dx, 0)
    numb_rep_length_0veg_rest_temp = np.round(length_0veg_rest / dx, 0)
    total_numb_rep_veg_pro_temp = numb_rep_length_0veg_begin +\
        numb_rep_width_veg +\
        numb_rep_length_0veg_rest_temp
    deviation_veg_pro = len(zp_sws) - total_numb_rep_veg_pro_temp
    numb_rep_length_0veg_rest = numb_rep_length_0veg_rest_temp +\
        deviation_veg_pro

    zp_0veg_begin = np.repeat(0, numb_rep_length_0veg_begin)
    zp_width_veg = np.repeat(1, numb_rep_width_veg)
    zp_0veg_rest = np.repeat(0, numb_rep_length_0veg_rest)
    zp_veg_temp = np.append(zp_0veg_begin, zp_width_veg)
    zp_veg = np.append(zp_veg_temp, zp_0veg_rest)
    # ------------------------------------------------------------------------
    # Plots
    plt.figure(figsize=(15, 5))
    plt.plot(xp_sws, zp_sws * -1, '-', color='black', linewidth=2)
    plt.plot(xp_sws, zp_veg, '--', color='green', linewidth=2)
    plt.axis([xp_sws[0], xp_sws[-1], -zp_sws[0], -zp_sws[-1] + 2])
    plt.xlabel('x-axis (m)')
    plt.ylabel('z-axis (m)')
    # ------------------------------------------------------------------------
    numb_mesh_veg = width_veg / dx
    print(f'\nNumber of mesh of the vegetation/fence area is \
{numb_mesh_veg:.0f}.')
    return zp_veg
