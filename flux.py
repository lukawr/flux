
"""
 Program for calculation of neutron flux based on experimental reaction rates 
 and drawing comparison between experimental values and results of simulation

 Version 0.1, Sun 24/07/2016

 Lukas Zavorka, lukas.zavorka@gmail.com
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline

# DATA FILES

home_dir = "/home/luke/python/data/"

# Cross sections, Talys 1.8:

input_cs_n_talys_27Al  = home_dir + "Al27_n_Talys18.txt"
input_cs_n_talys_59Co  = home_dir + "Co59_n_Talys18.txt"
input_cs_n_talys_197Au = home_dir + "Au197_n_Talys18.txt"
input_cs_n_talys_209Bi = home_dir + "Bi209_n_Talys18.txt"

# Spectra MCNPX (neutrons, protons)

input_flux_n_N14 = home_dir + "Spectra_n_N14.txt"
input_flux_n_A15 = home_dir + "Spectra_n_A15.txt"
input_flux_p_N14 = home_dir + "Spectra_p_N14.txt"
input_flux_p_A15 = home_dir + "Spectra_p_A15.txt"

# Experimental reaction rates

input_exp_rr = home_dir + "Experimental_Results.txt"

# FUNCTIONS

def open_file (file_name, column_name, number_of_columns):
    """Function for opening files and reading data into numpy arrays"""

    col_name = column_name
    num_cols = number_of_columns
    
    # Create list of column names
    
    scope = globals()
    variable_list = []
    
    for i in range (num_cols):
        input_data = []
        scope[col_name + str(i)] = input_data
        variable_list.append(input_data)

    # Read data
        
    with open (file_name, 'r') as f:
        lines = f.readlines()[2:]  #First two commented lines are skiped

    for line in lines:
        p = line.split()
        for i in range (num_cols):
            variable_list[i].append(float(p[i]))
 
    return 0 # Loaded values are stored in a list of global variables

def x_for_interp (inp_energy):

    in_ene = np.array([])
    in_ene = np.append(in_ene, inp_energy)
    in_ene_sel_index = np.where(in_ene > 1.) # only energies above 1 MeV

    in_ene_sel = np.array([])
    for indices in in_ene_sel_index:
        in_ene_sel = np.append(in_ene_sel, in_ene[indices])
    
    out_ene = np.array([]) # array of energies in which CS will be calculated
    out_num_of_int = np.array([]) # array with number of intervals
    
    for i in range (1,len(in_ene_sel)):
        if in_ene_sel[i] < 187:
            division_10 = (in_ene_sel[i] - in_ene_sel[i-1])/10. # 1/10 of each interval
            mid_point_10 = division_10 / 2.
            mid_points  = np.linspace( in_ene_sel[i-1] + mid_point_10 , in_ene_sel[i] - mid_point_10,10) # array of midpoints
            out_num_of_int = np.append(out_num_of_int, 10)

        else:
            number_of_intervals = int(in_ene_sel[i]-in_ene_sel[i-1]) / 2 # intervals a approx. 2 MeV
            step_2_MeV = float((in_ene_sel[i]-in_ene_sel[i-1]) / number_of_intervals)
            mid_points  = np.linspace( in_ene_sel[i-1] + (step_2_MeV / 2.) , in_ene_sel[i] - (step_2_MeV / 2.), number_of_intervals)
            out_num_of_int = np.append(out_num_of_int, number_of_intervals)
        out_ene = np.append(out_ene, mid_points)

    delta_e_ar = [] 
    delta_e_np = np.array([]) # delta E for each bin between boundaries 
    for i in range (1,len(in_ene_sel)):
        delta_e_ar.append((in_ene_sel[i] - in_ene_sel[i-1])/out_num_of_int[i-1])
    delta_e_np = np.append(delta_e_np, delta_e_ar)
        
    return (out_ene, in_ene_sel, out_num_of_int, delta_e_np)

def xs_interp (inp_ene, inp_xs, inp_ene_interp, plot_cs):

    inp_ene = inp_ene # energies from Talys
    inp_xs  = inp_xs  # xs from talys
    inp_ene_interps = inp_ene_interp # energies for interpolation
    out_xs_A = []
    out_xs = np.array([]) # iterpolated xs
    plot_fig = plot_cs
    
    x_ene = np.linspace (0,660,3301)
    
    spl = UnivariateSpline(inp_ene, inp_xs, s = 0.25)
    y_xs = spl(x_ene)

    for inp_ene_interp in inp_ene_interps:
        out_xs_A.append(spl.__call__(inp_ene_interp))
    
    out_xs = np.append(out_xs, out_xs_A)

    # optional_plot

    if plot_fig:
        plt.plot (inp_ene, inp_xs, 'ro', ms = 5)
        plt.plot (x_ene, y_xs, lw = 3, c = 'g', alpha = 0.6)
        plt.plot (inp_ene_interps, out_xs, 'o', ms = 3)
        plt.show()
    
    return out_xs

# MAIN BODY

""" xs_Al_[i] column structure: E(MeV),XS(mb),XS(mb)...
    0:E, 1:(n,a)24Na, 2:(n,2na)22Na, 3:(n,x)7Be """

open_file (input_cs_n_talys_27Al, "xs_Al_", 4)

""" xs_Co_[i] column structure: E(MeV),XS(mb),XS(mb)...
    0:E, 1:(n,g)60Co, 2:(n,2n)58Co, 3:(n,3n)57Co, 4:(n,4n)56Co, 5:(n,5n)55Co, 6:(n,p)59Fe"""

open_file (input_cs_n_talys_59Co, "xs_Co_", 7)

""" xs_Au_[i] column structure: E(MeV),XS(mb),XS(mb)...
    0:E, 1:(n,g)198Au, 2:(n,2n)196Au, 3:(n,3n)195Au, 4:(n,4n)194Au, 5:(n,5n)193Au, 6:(n,6n)192Au, 7:(n,7n)191Au"""

open_file (input_cs_n_talys_197Au, "xs_Au_", 8)

""" xs_Bi_[i] column structure: E(MeV),XS(mb),XS(mb)...
    0:E, 1:(n,3n)207Bi, 2:(n,4n)206Bi, 3:(n,5n)205Bi, 4:(n,6n)204Bi, 5:(n,7n)203Bi, 6:(n,8n)202Bi, 7:(n,9)201Bi, 8:(n,10n)200Bi"""

open_file (input_cs_n_talys_209Bi, "xs_Bi_", 9)

""" rr_[i] column structure: RR, dRR, RR, dRR...
    O:RR(1), 1:dRR(1), 2:RR(2), 3:dRR(2)..."""

open_file (input_exp_rr, "rr_", 40)

# Names of the reaction products for the rows in RR_[i] are as follows:

rr_products = ['24Na', '22Na', '7Be', '60Co', '59Fe', '58Co', '57Co', '56Co', '55Co', '198Au', '196Au', '195Au', '194Au', '193Au', '192Au', '191Au', '207Bi', '206Bi', '205Bi', '204Bi', '203Bi', '202Bi', '201Bi']

""" Flux_N14_[i] column structure: E(MeV), flux(cm-2particle-1), dflux(), ...
   0:E(1), 1:flux(1), 2:flux(1), 0:E(2), 1:flux(2), 2:flux(2), ... 
   Last few lines in some columns replaced by -1.0 to maintain the size"""

open_file (input_flux_n_N14, "flux_n14_", 60)

""" Flux_N14_[i] column structure: E(MeV), flux(cm-2particle-1), dflux(), ...
   0:E(1), 1:flux(1), 2:flux(1), 0:E(2), 1:flux(2), 2:flux(2), ...
   Last few lines in some columns replaced by -1.0 to maintain the size"""

open_file (input_flux_n_A15, "flux_n15_", 48)

""" Flux_N14_[i] column structure: E(MeV), flux(cm-2particle-1), dflux(), ...
   0:E(1), 1:flux(1), 2:flux(1), 0:E(2), 1:flux(2), 2:flux(2), ...  """

open_file (input_flux_p_N14, "flux_p14_", 60)

""" Flux_N14_[i] column structure: E(MeV), flux(cm-2particle-1), dflux(), ...
   0:E(1), 1:flux(1), 2:flux(1), 0:E(2), 1:flux(2), 2:flux(2), ...  """

open_file (input_flux_p_A15, "flux_p15_", 48)

# Collection of necessary data for further calculation
# Sample = Al3 = cell 373 N14

input_flux_ene = flux_n14_36 # A propper file with energies <flux_p/n14/15_i> i = 0,3,6,...
input_reactions = [xs_Co_6, xs_Al_1, xs_Au_2, xs_Co_2, xs_Bi_1, xs_Co_3, xs_Au_4, xs_Au_5, xs_Co_4, xs_Au_6, xs_Au_7]
input_rr_exp  = [rr_4[4],rr_4[0],rr_4[10],rr_4[5],rr_4[16],rr_4[6],rr_4[12],rr_4[13],rr_4[7],rr_4[14],rr_4[15]]
input_drr_exp = [rr_5[4],rr_5[0],rr_5[10],rr_5[5],rr_5[16],rr_5[6],rr_5[12],rr_5[13],rr_5[7],rr_5[14],rr_5[15]]

# Abstraction of energy points for interpolation of cross sections, arrays of energy boundaries between reactions and number of points per each interval 

x_for_interpol = np.array([])    # energies (mid-bins) for which interpolation is calculated
energy_boundaries = np.array([]) # energy boundaries of xs (thresholds) and neutron flux
energy_intervals = np.array([])  # number of energy bins between the given boundaries

x_for_interpol, energy_boundaries, energy_intervals, energy_bin_widths = x_for_interp(input_flux_ene) 
energy_intervals_cumulative = np.cumsum(energy_intervals)
energy_intervals_cumulative = np.insert(energy_intervals_cumulative,0,0)


# Interpolation of cross section for given eneries

plot_yn = False # on/off plot of interpolated cross sections
input_energy = xs_Al_0 # energies from talys

scope = globals()
list_xs_interp = [] # list members are np.arrays of interpolated xs, order as in "input_reactions"

for i in range(len(input_reactions)):
    xs_intrplt = np.array([])
    scope['xs_interpolated_' + str(i)] = xs_intrplt
    list_xs_interp.append(xs_intrplt)

for i in range(len(list_xs_interp)):
    list_xs_interp[i] = xs_interp (input_energy, input_reactions[i], x_for_interpol, plot_yn)

    
# EXPERIMENTAL FLUX CALCULATION

scope = globals()
flux_list = [] # declaration of list of experimental fluxes 

for i in range(len(input_reactions)):
    flux_i = float() 
    scope['flux_' + str(i)] = flux_i
    flux_list.append(flux_i)

for i in range(len(flux_list)-1,9,-1): # descending order, does not include (n,g) reaction
    #first attempt - the last reacion
    print(i)
    sum_xsE = ( np.sum(list_xs_interp[i][energy_intervals_cumulative[i]:energy_intervals_cumulative[i+1]]) ) * energy_bin_widths[i]
    #print sum_xsE
    flux_list[i] = (input_rr_exp[i]/sum_xsE)*1e27

    
# TEST:

print (energy_intervals)
print (energy_intervals_cumulative)
print (energy_bin_widths)
print (len(flux_list))
print (len(list_xs_interp[10]))
print (flux_list)


    
