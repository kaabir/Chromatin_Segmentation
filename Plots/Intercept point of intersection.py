from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import matplotlib.colors as mcolors
from decimal import Decimal
import numpy as np
from scipy.optimize import curve_fit

def linear_func(x, a, b): # fg = ax or b =ax
    return a * x + b


def plot_intersections(filename,f_lift_values, f_grav_values, tol=None):

    if len(F_grav)==5:
        velocity_values = [0, 2.5E-11, 5.0E-11, 7.5E-11, 1.0E-11]
        velocity = np.array(velocity_values)
    elif len(F_grav)==4:  
        velocity_values = [0, 2.5E-11, 5.0E-11, 7.5E-11]
        velocity = np.array(velocity_values)
    elif len(F_grav)==3:  
        velocity_values = [0, 2.5E-11, 5.0E-11]
        velocity = np.array(velocity_values) 
    elif len(F_grav)==2:  
        velocity_values = [0, 2.5E-11]
        velocity = np.array(velocity_values)          
        
    idx_inter = np.argwhere(np.abs(f_lift_values - f_grav_values) < tol).flatten()

    if len(idx_inter) == 0:
    # No intersection points found within the given data
        if len(f_grav_values) > 0:
        # Fit a line to the data and extrapolate to predict the intersection
            popt, _ = curve_fit(linear_func, f_lift_values, velocity)
            x_inter = linear_func(f_grav_values[0], *popt)
            print('For ' + filename + ':', x_inter)
        else:
            print('No intersection points found for ' + filename)
        return
    else:
    # Intersection point found within the given data
        x_inter = np.interp(f_grav_values[idx_inter], f_lift_values, velocity)
        print('For ' + filename + ':', x_inter[0])

    interp_func = interp1d(f_lift_values, velocity)
    y_inter = interp_func(f_grav_values[idx_inter[0]]) if len(idx_inter) > 0 else None        
    # idx_inter = np.argwhere(np.abs(F_lift - F_grav) < tol).flatten()
    # x_inter = np.interp(F_grav[idx_inter], F_lift, velocity)
    # print('For '+ filename + ':',x_inter[0])
    # interp_func = interp1d(F_lift, velocity)
    # y_inter = interp_func(F_grav[idx_inter[0]])
    
    # annot_text_x = Decimal(x_inter[0]) #Decimal("%.2f" % x_inter[0])
    # annot_text_y = Decimal(F_grav[idx_inter][0])
    
    fig, ax = plt.subplots(figsize=(6, 4))

    ax.plot(velocity, F_grav, marker='o', alpha=0.8, color="#11aa00", mec='none', ms=3, lw=1, label='F_Gravity')
    ax.plot(velocity, F_lift, marker='o', alpha=0.4, color="#11aa00",mec='none', ms=3, lw=1, label='F_Lift')
    #plt.plot(x_inter, F_grav[idx_inter], 'ms', ms=5, label='Intersection on x-axis')
    ax.plot(y_inter, F_grav[idx_inter[0]], marker='o', color="#fc033d",mec='none', label='Intersection')
    ax.set_title('Vesicle Simulation Result ' + filename)

    ax.set_ylabel(r'$Force\ (Gravity\ - Lift)\ [N] $')
    ax.set_xlabel(r'$Flow\ rate\ [m^{3}/s]$')
    #ax.text(annot_text_x, annot_text_y, '({}, {})'.format(annot_text_x, annot_text_y))

    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.legend()
    ax.spines[['right', 'top']].set_visible(False)  
    ax.grid(which='minor', color='#dbdbdb', linestyle=':')
    ax.grid(which='major', color='#c2c0c0', linestyle='--')
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    ##plt.savefig(filename+'.png', dpi=400, bbox_inches="tight", pad_inches=0.1)
    plt.show()
    ##print(f'Intersection point on (velocity,Force): ({x_inter[0]}, {F_grav[idx_inter][0]})')
    #print(f'Intersection point on y-axis: ({y_inter}, {F_grav[idx_inter][0]})')



# Values of 55um trap - 4um vesicle
F_lift = np.array([0,1.4293471958570603E-15,2.858655813332324E-15,
                   4.287925663750874E-15,5.7171565563729E-15])
F_grav = np.repeat(1.2159087840683463E-14, 5)

plot_intersections('55um_4um_vesicle',F_lift, F_grav, tol=1e-15)

# Values of 55um trap - 6um vesicle
F_lift = np.array([0,6.2659727558490745E-15,1.2531787203835203E-14,
                   1.8797442982508422E-14,2.5062939742331027E-14])
F_grav = np.repeat(4.103692146230669E-14, 5)

plot_intersections('55um_6um_vesicle',F_lift, F_grav, tol=1e-15)

# Values of 55um trap - 8um vesicle
F_lift = np.array([0,1.7348874155385515E-14,3.469737867469837E-14,
                   5.2045513377565575E-14,6.939327809547115E-14])
F_grav = np.repeat(9.727270272546771E-14, 5)

plot_intersections('55um_8um_vesicle',F_lift, F_grav, tol=1e-14)

# Values of 55um trap - 10um vesicle
F_lift = np.array([0,3.7380838330963224E-14,7.476101565739288E-14,
                   1.1214053287567826E-13,1.495193909252081E-13])
F_grav = np.repeat(1.8998574751067906E-13, 5)

plot_intersections('55um_10um_vesicle',F_lift, F_grav, tol=1e-14)

# Values of 55um trap - 12um vesicle
F_lift = np.array([0,7.001941226691884E-14,1.4003747873072283E-13,
                   2.1005420319468573E-13,2.80069589635807E-13])
F_grav = np.repeat(3.282953716984535E-13, 5)

plot_intersections('55um_12um_vesicle',F_lift, F_grav, tol=1e-14)

# Values of 55um trap - 14um vesicle
F_lift = np.array([0,1.1755007649250795E-13,2.3509861519821266E-13,
                   3.5264562600101114E-13,4.701911187612114E-13])
F_grav = np.repeat(5.213208911693032E-13, 5)

plot_intersections('55um_14um_vesicle',F_lift, F_grav, tol=1e-14)

# Values of 55um trap - 16um vesicle
F_lift = np.array([0,1.8686870006424347E-13,3.737352502211608E-13,
                   5.605996712108771E-13,7.474619834003568E-13])
F_grav = np.repeat(7.781816218037417E-13, 5)

plot_intersections('55um_16um_vesicle',F_lift, F_grav, tol=1e-14)

# Values of 55um trap - 18um vesicle
F_lift = np.array([0,2.807856813851369E-13,5.615684801148335E-13,
                   8.423484360479435E-13,1.1231255888212665E-12])
F_grav = np.repeat(1.1079968794822806E-12, 5)

plot_intersections('55um_18um_vesicle',F_lift, F_grav, tol=1e-13)

# Values of 55um trap - 20um vesicle
F_lift = np.array([0,4.1407264599218304E-13,8.28143358952309E-13,
                   1.2422122103773347E-12,1.6562792714079788E-12])
F_grav = np.repeat(1.5198859800854325E-12, 5)

plot_intersections('55um_20um_vesicle',F_lift, F_grav, tol=1e-13)

# Values of 55um trap - 22um vesicle
F_lift = np.array([0,5.849005057273633E-13,1.1697991038997446E-12,
                   1.7546959177172467E-12,2.3395910703130258E-12])
F_grav = np.repeat(2.0229682394937113E-12, 5)

plot_intersections('55um_22um_vesicle',F_lift, F_grav, tol=1e-12)

# Values of 55um trap - 24um vesicle
F_lift = np.array([0,8.169597815813077E-13,1.6339191836844742E-12,
                   2.4508784079156326E-12,3.2678376557058614E-12])
F_grav = np.repeat(2.6263629735876273E-12, 5)

plot_intersections('55um_24um_vesicle',F_lift, F_grav, tol=1e-12)

# Values of 55um trap - 26um vesicle
F_lift = np.array([0,1.1100676293390336E-12,2.220138955106236E-12,
                   3.330214301686069E-12,4.44029399271545E-12])
F_grav = np.repeat(3.3391894982476956E-12, 5)

plot_intersections('55um_26um_vesicle',F_lift, F_grav, tol=1e-12)

# Values of 55um trap - 28um vesicle
F_lift = np.array([0,1.4753143522917832E-12,2.950625151060579E-12,4.425932894421467E-12])
F_grav = np.repeat(4.170567129354429E-12,4)

plot_intersections('55um_28um_vesicle',F_lift, F_grav, tol=1e-12)

# Values of 55um trap - 30um vesicle
F_lift = np.array([0,1.9136410568381156E-12,3.827270298014552E-12,5.7408884666722955E-12])
F_grav = np.repeat(5.129615182788335E-12, 4)

plot_intersections('55um_30um_vesicle',F_lift, F_grav, tol=1e-12)

# Values of 55um trap - 32um vesicle
F_lift = np.array([0,2.4195797293313793E-12,4.8391262771260756E-12,7.258640725972604E-12])
F_grav = np.repeat(6.225452974429933E-12, 4)

plot_intersections('55um_32um_vesicle',F_lift, F_grav, tol=1e-12)

# Values of 55um trap - 34um vesicle
F_lift = np.array([0,2.9299336484883366E-12,5.859777287085218E-12,8.789532430332764E-12])
F_grav = np.repeat(7.467199820159733E-12, 4)

plot_intersections('55um_34um_vesicle',F_lift, F_grav, tol=1e-12)

# Values of 55um trap - 36um vesicle
F_lift = np.array([0,3.2276791304752964E-12,6.4551746068589114E-12,9.68248833974238E-12])
F_grav = np.repeat(8.863975035858245E-12, 4)

plot_intersections('55um_36um_vesicle',F_lift, F_grav, tol=1e-12)


# # Values of 55um trap - 38um vesicle
# F_lift = np.array([0, 3.689799909895379E-14, 7.30265454987984E-14, 1.0838563178899489E-13, 1.4297541514406575E-13])
# F_grav = np.repeat(1.0424897937405988E-11, 5)

# plot_intersections('55um_38um_vesicle',F_lift, F_grav, tol=1e-15)


import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

def linear_func(x, a, b): # fg = ax or b =ax
    return a * x + b

def plot_intersections(filename, f_lift_values, f_grav_values, tol=None):

    velocity_values = [0, 2.5E-12, 5.0E-12, 7.5E-12, 1.0E-11]
    velocity = np.array(velocity_values)

    idx_inter = np.argwhere(np.abs(f_lift_values - f_grav_values) < tol).flatten()

    if len(idx_inter) == 0:
        # No intersection points found within the given data
        if len(f_grav_values) > 0:
            # Fit a line to the data and extrapolate to predict the intersection
            popt, _ = curve_fit(linear_func, f_lift_values, velocity)
            x_inter = linear_func(f_grav_values[0], *popt)
            print('For ' + filename + ':', x_inter)
        else:
            print('No intersection points found for ' + filename)
            return
    else:
        # Intersection point found within the given data
        x_inter = np.interp(f_grav_values[idx_inter], f_lift_values, velocity)
        print('For ' + filename + ':', x_inter[0])

    interp_func = interp1d(f_lift_values, velocity)
    y_inter = interp_func(f_grav_values[idx_inter[0]]) if len(idx_inter) > 0 else None

    fig, ax = plt.subplots(figsize=(6, 4))

    ax.plot(velocity, f_grav_values, marker='o', alpha=0.8, color="#11aa00", mec='none', ms=3, lw=1, label='F_Gravity')
    ax.plot(velocity, f_lift_values, marker='o', alpha=0.4, color="#11aa00",mec='none', ms=3, lw=1, label='F_Lift')
    if y_inter is not None:
        ax.plot(y_inter, f_grav_values[idx_inter[0]], marker='o', color="#fc033d",mec='none', label='Intersection')
    ax.set_title('Vesicle Simulation Result ' + filename)

    ax.set_ylabel(r'$Force\ (Gravity\ - Lift)\ [N] $')
    ax.set_xlabel(r'$Flow\ rate\ [m^{3}/s]$')

    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.legend()
    ax.spines[['right', 'top']].set_visible(False)  
    ax.grid(which='minor', color='#dbdbdb', linestyle=':')
    ax.grid(which='major', color='#c2c0c0', linestyle='--')
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))

    plt.show()
    