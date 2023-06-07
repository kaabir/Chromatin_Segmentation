import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import matplotlib.colors as mcolors
from decimal import Decimal

def plot_intersections(filename,f_lift_values, f_grav_values, tol=None):

    if len(F_grav)==5:
        velocity_values = [0, 2.5E-12, 5.0E-12, 7.5E-12, 1.0E-11]
        velocity = np.array(velocity_values)
    elif len(F_grav)==4:  
        velocity_values = [0, 2.5E-12, 5.0E-12, 7.5E-12]
        velocity = np.array(velocity_values)
    elif len(F_grav)==3:  
        velocity_values = [0, 2.5E-12, 5.0E-12]
        velocity = np.array(velocity_values)    
        
    idx_inter = np.argwhere(np.abs(F_lift - F_grav) < tol).flatten()
    x_inter = np.interp(F_grav[idx_inter], F_lift, velocity)
    print('For '+ filename + ':',x_inter[0])
    interp_func = interp1d(F_lift, velocity)
    y_inter = interp_func(F_grav[idx_inter[0]])
    
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

# Values of 43um trap - 4um vesicle
F_lift = np.array([0,3.818327355512234E-15,7.636565334253353E-15,1.1454714163530053E-14,1.527277407339555E-14])
F_grav = np.repeat(1.2159087840683463E-14, 5)

plot_intersections('43um_4um_vesicle',F_lift, F_grav, tol=1e-15)

# Values of 43um trap - 6um vesicle
F_lift = np.array([0,1.5592233140140826E-14, 3.118416934397372E-14, 4.677581053096509E-14])
F_grav = np.repeat(4.103692146230669E-14, 4)

plot_intersections('43um_6um_vesicle',F_lift, F_grav, tol=1e-14)

# Values of 43um trap - 8um vesicle
F_lift = np.array([0,4.155550038943612E-14, 8.311034744879125E-14, 1.2466454822807502E-13])
F_grav = np.repeat(9.727270272546771E-14, 4)

plot_intersections('43um_8um_vesicle',F_lift, F_grav, tol=1e-13)

# Values of 43um trap - 10um vesicle
F_lift = np.array([0,8.659729650830929E-14, 1.731937421682183E-13, 2.5978935517324326E-13])
F_grav = np.repeat(1.8998574751067906E-13, 4)

plot_intersections('43um_10um_vesicle',F_lift, F_grav, tol=1e-13)

# Values of 43um trap - 12um vesicle
F_lift = np.array([0,1.5650358011481465E-13, 3.1300637994429583E-13, 4.695084371861038E-13])
F_grav = np.repeat(3.282953716984535E-13, 4)

plot_intersections('43um_12um_vesicle',F_lift, F_grav, tol=1e-13)

# Values of 43um trap - 14um vesicle
F_lift = np.array([0,2.546849731839487E-13, 5.093695671715059E-13, 7.64053853093081E-13])
F_grav = np.repeat(5.213208911693032E-13, 4)

plot_intersections('43um_14um_vesicle',F_lift, F_grav, tol=1e-13)

# Values of 43um trap - 16um vesicle
F_lift = np.array([0, 3.8838939170646027E-13, 7.767792251761394E-13, 1.1651696252070833E-12])
F_grav = np.repeat(7.781816218037417E-13, 4)

plot_intersections('43um_16um_vesicle',F_lift, F_grav, tol=1e-13)

# Values of 43um trap - 18um vesicle
F_lift = np.array([0, 5.686363607836688E-13, 1.1372756979425596E-12])
F_grav = np.repeat(1.1079968794822806E-12, 3)

plot_intersections('43um_18um_vesicle',F_lift, F_grav, tol=1e-13)

# Values of 43um trap - 20um vesicle
F_lift = np.array([0, 7.982561059852858E-13, 1.5965193166022287E-12])
F_grav = np.repeat(1.5198859800854325E-12, 3)

plot_intersections('43um_20um_vesicle',F_lift, F_grav, tol=1e-13)

# Values of 43um trap - 22um vesicle
F_lift = np.array([0, 1.0907047292278014E-12, 2.1814226933722858E-12])
F_grav = np.repeat(2.0229682394937113E-12, 3)

plot_intersections('43um_22um_vesicle',F_lift, F_grav, tol=1e-12)

# Values of 43um trap - 24um vesicle
F_lift = np.array([0, 1.441835665408444E-12, 2.88368284173721E-12])
F_grav = np.repeat(2.6263629735876273E-12, 3)

plot_intersections('43um_24um_vesicle',F_lift, F_grav, tol=1e-12)

# Values of 43um trap - 26um vesicle
F_lift = np.array([0, 1.859142729101619E-12, 3.718287498560558E-12])
F_grav = np.repeat(3.3391894982476956E-12, 3)

plot_intersections('43um_26um_vesicle',F_lift, F_grav, tol=1e-12)

# Values of 43um trap - 28um vesicle
F_lift = np.array([0, 2.330813571228374E-12, 4.661594522894895E-12])
F_grav = np.repeat(4.170567129354429E-12, 3)

plot_intersections('43um_28um_vesicle',F_lift, F_grav, tol=1e-12)

# Values of 43um trap - 30um vesicle
F_lift = np.array([0, 2.810451958885426E-12, 5.620813779266104E-12])
F_grav = np.repeat(5.129615182788335E-12, 3)

plot_intersections('43um_30um_vesicle',F_lift, F_grav, tol=1e-12)

# Values of 43um trap - 32um vesicle
F_lift = np.array([0, 3.27187135741756E-12, 6.543531800872797E-12])
F_grav = np.repeat(6.225452974429933E-12, 3)

plot_intersections('43um_32um_vesicle',F_lift, F_grav, tol=1e-12)

# Values of 43um trap - 34um vesicle
F_lift = np.array([0, 3.548080731552862E-12, 7.0957744871031625E-12, 1.0643085337462863E-11])
F_grav = np.repeat(7.467199820159733E-12, 4)

plot_intersections('43um_34um_vesicle',F_lift, F_grav, tol=1e-12)

# Values of 43um trap - 36um vesicle
F_lift = np.array([0, 3.378764181857352E-12, 6.7569086751144605E-12, 1.013443794706224E-11])
F_grav = np.repeat(8.863975035858245E-12, 4)

plot_intersections('43um_36um_vesicle',F_lift, F_grav, tol=1e-11)

# Values of 43um trap - 38um vesicle
F_lift = np.array([0, 3.689799909895379E-14, 7.30265454987984E-14, 1.0838563178899489E-13, 1.4297541514406575E-13])
F_grav = np.repeat(1.0424897937405988E-11, 5)

plot_intersections('43um_38um_vesicle',F_lift, F_grav, tol=1e-15)