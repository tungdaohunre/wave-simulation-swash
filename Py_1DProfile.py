# Check wave period and calculate resolution


def check_T_for_resolution(H, T, d0):
    import numpy as np
    # -----------------------------------------------------------------------#
    # Check whereas wave period is good to calculate proper resolution       #
    # And calculate resolution from chosen wave period                       #
    # Input: Wave period T, water depth                                      #
    # Output: Resolution dx                                                  #
    # -----------------------------------------------------------------------#
    T = np.array(T)
    count_T = len(np.atleast_1d(T))
    if count_T >= 2:
        Tmin_id = np.argmin(T)
        Tmin = T[Tmin_id]
        print(f'\nAt boundary:\
            \n Initial water depth is {d0} m.\
            \n Initial wave height is {H} m.\
            \n Minimum wave period is {Tmin} s.')
    else:
        Tmin = T
        print(f'\nAt boundary:\
            \n Initial water depth is {d0} m.\
            \n Initial wave height is {H} m.\
            \n Initial wave period is {Tmin} s.\n')

    # Calculate dx, resolution for the smallest wave period
    waterzone_deep = 0.5                     # Set min d/L for deep water
    waterzone_shallow = 0.05                 # Set max d/L for shallow water
    L0, L = wavelength(Tmin, d0)
    dL = d0/L
    if dL > waterzone_deep:
        print('Warning 1:\
            \n Deep water zone.\
            \n Resolution (dx) is too large. \
            \n Reducing water depth by half; \
            \n Calculate again wavelength')
        d1 = d0*(1/2)
        L0, L1 = wavelength(Tmin, d1)
        dx = np.round(L1/100, 2)
        print(f'Result: L = {L1:.2f} m, and dx is {dx} m.')
    elif waterzone_shallow < dL <= waterzone_deep:
        print('Warning 2:\
            \n Transitional water zone.\
            \n Resolution (dx) is still large. \
            \n Reducing water depth by half;\
            \n Calculate again wavelength')
        d2 = d0*(1/2)
        L0, L2 = wavelength(Tmin, d2)
        dx = np.round(L2/100, 2)
        print(f'Result: L = {L2:.2f} m, and dx is {dx} m.')
    else:
        dx = np.around(L/100, 2)
        print(f'Warning 3:\
            \n Shallow water zone.\
            \nResult: L = {L:.2f} m, and dx is {dx} m.')
    return dx, L0, L

# generate model resolution


def create_xp_zp_sws(dx, d0, i_slope, len_slope):
    import numpy as np
    import matplotlib.pyplot as plt
    # -----------------------------------------------------------------------#
    # Input: Resolution (dx), profile slope(s) (i_slope),                    #
    #        length of slope(s) (len_slope)                                  #
    # Output: Mesh(es) in x-cordinate (xp)                                   #
    #         Mesh(es) in z-cordinate (zp), it is also bottom level.         #
    # -----------------------------------------------------------------------#
    # Create array of xp and zp with step dx
    # Mesh in x-cordinate, xp
    i_slope = np.array(i_slope)
    i_slope = np.append(0, i_slope)
    len_slope = np.array(len_slope).astype(np.float64)
    len_slope = np.append(0, len_slope)
    x = np.cumsum(len_slope).astype(np.float64)
    xp = []
    number_slope_p = []
    i_slope_p = []
    for i1 in range(len(x)):
        xp = np.arange(0, x[i1], dx)
        dx_p = np.repeat(dx, len(xp))
        temp_number_slope_p = len(np.arange(0, len_slope[i1], dx))
        number_slope_p = np.append(number_slope_p, temp_number_slope_p)
        temp_i_slope_p = np.repeat(i_slope[i1], number_slope_p[i1])
        i_slope_p = np.append(i_slope_p, temp_i_slope_p)
    differ_number = len(i_slope_p) - len(dx_p)
    repeat_differ_dx = np.repeat(dx, differ_number)
    repeat_differ_xp = np.repeat(xp[-1], differ_number)
    dx_p = np.append(dx_p, repeat_differ_dx)         # update dx_p
    xp = np.append(xp, repeat_differ_xp)             # update xp
    # Mesh in x-cordinate, zp
    temp_z = d0
    i_dx = i_slope_p * dx_p
    zp = np.full((i_dx.shape[0], 2), np.inf, dtype=(np.float64))
    for i3 in range(len(i_dx)):
        # Calculation
        z_i = temp_z - i_dx[i3]
        # Assign
        zp[i3] = [temp_z, z_i]
        # Update for next interation
        temp_z = z_i
    zp = zp[:, 0]
    # Extend one value to make sure that the number of meshes is even
    zp = np.append(zp, zp[-1])
    xp = np.append(xp, xp[-1])
    # Print results
    print(f'\nNumber of mesh is {len(zp)}\
          \nTotal length of profile is {(x[-1])} m.')
    # PLot profile
    plt.figure(figsize=(15, 5))
    plt.plot(xp, zp * -1, '-', color='black', linewidth=2)
    plt.plot([x[0], x[-1]], [0, 0], '--', color='blue', linewidth=1.2)
    plt.axis([x[0], x[-1], -zp[0], -zp[-1] + 2])
    plt.xlabel('x-axis (m)')
    plt.ylabel('z-axis (m)')
    # ------------------------------------------------------------------------
    return xp, zp


def wavelength(T, d):
    import numpy as np
# ------------------------------------------------------------------------
# Wave length calculation based on period (T) and water depth (d)        #
# Input: wave period (T); Water depth (d)                                #
# Output: Wavelength at deep water: d0                                   #
#         Wavelength at chosen water: d                                  #
# ------------------------------------------------------------------------
    # wavelength at deep water
    g = 9.81
    L0 = (g*T**2)/(2*np.pi)
    guess = L0
    L = (g*T**2)/(2*np.pi)*np.tanh((2*np.pi)*(d/guess))
    diff = abs(L-guess)
    # wavelength at intermediate water depth and shallow water
    while diff > 0.01:
        diff = abs(L-guess)
        guess = L + (0.5*diff)
        L = (g*T**2)/(2*np.pi)*np.tanh((2*np.pi)*(d/guess))
    return L0, L
