from numpy import *
from curve_is_fine_gas_amplification import *
from histogram_fit_model import *
from general_fit import *
from zero_signal import *
from matplotlib.pyplot import *
from sys import argv
from scipy.stats import *
import sys
import os

def exponential_ln_model(parameters, data):
    alpha = parameters[0]
    beta = parameters[1]
    return alpha*exp(beta * log(data))

def exponential_ln_confidence_band_point(x_point, parameters, covariance_matrix, reduced_chi_squared):
    alpha = parameters[0]
    beta = parameters[1]
    f1 = (exp(beta * log(x_point)))**2
    f2 = exp(beta * log(x_point)) * alpha * log(x_point) * exp(beta * log(x_point))
    f3 = f2
    f4 = (alpha * log(x_point) * exp(beta * log(x_point)))
    summing = reduced_chi_squared **2 * f1 * covariance_matrix[0][0] + 2*f2*covariance_matrix[0][1] + f4 * covariance_matrix[1][1]
    return sqrt(abs(summing))

def exponential_ln_confidence_band(x_data, y_data, parameters, covariance_matrix, reduced_chi_squared, degrees_of_freedom):
    alpha = 0.05
    tval = t.ppf(1.0 - alpha/2.0, degrees_of_freedom)
    sigma = zeros_like(x_data)
    for i in xrange(x_data.shape[0]):
        sigma[i] = exponential_ln_confidence_band_point(x_data[i], parameters, covariance_matrix, reduced_chi_squared)
    return y_data + sigma, y_data - sigma

wire_radius = 25e-6

channels_to_be_considered = [1]

time_uncertainty = 15e-12 # According to the documentation
# of the DRS4, after calibration the error in the time-axis
# is reduced to at most 15 ps. This is a systematic error,
# since the measurement of each time-event can be done once.
# This error will therefore have to be added at the very end.

# All times here are written in nanoseconds, all distance in
# milimeters, as otherwise numerical errors start to appear.
# The conversion for the different results will be done at
# the very end.
# In particular, for practical purposes, we need to rewrite
# the time_uncertainty in ns:
time_uncertainty = time_uncertainty * 1.0e9

voltage_zero_level, voltage_uncertainty = find_noise_level() # Voltage uncertainty for the DRS
voltage_anode = 3.0 # Anode voltage for the drift velocity experiment
voltage_cathode = 4.0 # Cathode voltage for the strontium experiments
voltage_uncertainty_anode = sqrt(2) * 0.1 # Voltage uncertainty for the anode or cathode, in kV.
# Since we are adding the anode voltage and the cathode voltage, the error has to be mutliplied with
# sqrt(2), since we are not going to use the error for anything else than the sum of it

to_treat_batch1 = ["Strontium_2000_500", "Strontium_2200_250", "Strontium_2400_250", "Strontium_2500_400",
            "Strontium_2700_500"]
to_treat_batch2 = ["Strontium_2900_500_1100mbar","Strontium_3000_500_1100mbar", "Strontium_3100_1000_1100mbar"]
to_treat_batch3 = ["Strontium_2900_500_1140mbar",  "Strontium_3000_500_1140mbar", "Strontium_3100_500_1140mbar"]
to_treat_batch_drift = ["Drift_1000_1000", "Drift_2000_1000", "Drift_2500_1000", "Drift_3000_1000", "Drift_3500_1000", "Drift_3800_1000", "Drift_4000_1000"]


if len(argv) == 1:
    print "No arguments given. You have to use it as follows:"
    print "python gas_amplification.py base_name_directory1 base_name_directory2 ..."
    print " or"
    print "python gas_amplification.py -all"
    exit()

if argv[1] != "-all":
    to_treat = argv[1:]
    flag_all = False
else:
    flag_all = True
    mean_peak_height_all = []
    mean_peak_height_error_all = []
    max_peak_height_all = []
    max_peak_height_error_all = []
    voltages_all = []
    amount_of_peaks = []
    amount_of_peaks_proportion = []
    to_treat = to_treat_batch1 + to_treat_batch2 + to_treat_batch3 + to_treat_batch_drift
    
    selection_to_treat = [to_treat_batch1, to_treat_batch2, to_treat_batch3, to_treat_batch_drift]
    selection_fit1 = [[0,1,2,3,4], [0], [0], []]
    selection_fit3 = [[],[],[],[0,1,2,3,4,5,6]]


for z, main_filename in enumerate(to_treat):
    _ = os.system('clear')
    print "Starting now with " + main_filename
    sys.stdout.flush()
    data_path = "../Gas_Amplification/" + main_filename + "_Data/"
    histogram_path = data_path + main_filename + "_Gas_Amplification/"
    if "Drift" in main_filename:
        data_path = data_path = "../Drift_Velocities/" + main_filename + "_Data/"
    histogram_path = data_path + main_filename + "_Gas_Amplification/"
    if not os.path.exists(histogram_path):
        os.makedirs(histogram_path)
    results = open(histogram_path + main_filename + "_peak_heights", "w")
    mean_peak_height_results = open(histogram_path + main_filename + "_mean_peak_height", "w")
    peak_height = []

    channels_voltages = []
    print "o Starting to read in all channels...",
    sys.stdout.flush()
    for chn in channels_to_be_considered:
        channels_voltages.append(loadtxt(data_path+main_filename + "_chn{}_v".format(chn)))
    print "Done reading."
    print "o-> Starting to process the data...",
    sys.stdout.flush()
    number_of_events = channels_voltages[0].shape[0]
    for i in xrange(number_of_events):
        total_peak_height = 0.0
        all_curves_fine = True
        for j, chn in enumerate(channels_to_be_considered):
            is_fine, height = curve_is_fine(channels_voltages[j][i])
            if is_fine:
                total_peak_height += height
            else:
                all_curves_fine = False
                break
        if all_curves_fine and total_peak_height != 0.0:
            peak_height.append(total_peak_height)
    for entry in peak_height:
        results.write(str(entry)+" ")
    results.write("\n")

    parameters_ideal, parameters_errors, parameters_p_values, SSRs, bin_heights, bin_heights_mean_center, bin_heights_std_center = histogram_fit_model(array(peak_height), 15, "Relative occurences for peak_heights", [], [], [], [], [], [], full_output = True)

    legend(loc = "best")
    if "Drift" in main_filename:
        title("Peak-height distribution\n with a Cathode voltage of " + main_filename[6:10] + r"$\pm$100V")
    else:
        title("Peak-height distribution\nwith a Cathode voltage of " + main_filename[10:14] + r"$\pm$100V")
    xlabel("Height of peak in V")
    savefig(histogram_path + main_filename + '_peak_height_distribution.png')
    clf()
        
    print "Done."
    print ""
    sys.stdout.flush()
        

    mean_peak_height = array(peak_height).mean()
    mean_peak_height_error = sem(array(peak_height))

    mean_peak_height_results.write("The mean peak height is " + str(mean_peak_height) + "+-" + str(mean_peak_height_error) + "V")
    
    mean_peak_height_all.append(mean_peak_height)
    mean_peak_height_error_all.append(mean_peak_height_error+voltage_uncertainty)
    if "Drift" in main_filename:
        voltages_all.append(float(main_filename[6:10]))
    else:
        voltages_all.append(float(main_filename[10:14]))

    # For the maximum peak: We are actually interested in the maximum peak
    # of the pdf for the distribution of the peak heights (so not the
    # peak values directly). First, we thus have to determine where this
    # maximum is located:
    if len(bin_heights) == 0:
        maximum_voltage = max(peak_height)
        maximum_voltage_error = voltage_uncertainty
        # Basically, this is only needed if the data is extremely bad, 
        # as can happen for the higher voltages
    else:
        maximum_point = where(bin_heights == max(bin_heights))[0][0]
        # Typically, there should only be one value. If not, the first maximum
        # is taken because, hey, why not.

        # Now, we get the actual voltage corresponding to this peak
        maximum_voltage = bin_heights_mean_center[maximum_point]
        maximum_voltage_error = bin_heights_std_center[maximum_point]

    max_peak_height_all.append(maximum_voltage)
    max_peak_height_error_all.append(maximum_voltage_error + voltage_uncertainty)

    amount_of_peaks.append(len(peak_height))
    amount_of_peaks_proportion.append(int(len(peak_height)/(number_of_events+0.0)*100.0))

    results.close()
    mean_peak_height_results.close()

    print "ooo Done with " + main_filename + ". ooo"
    print ""
    print ""
    print ""
    sys.stdout.flush()
_ = os.system('clear')

if flag_all:
    cut1 = len(to_treat_batch1)
    cut2 = len(to_treat_batch1) + len(to_treat_batch2)
    cut3 = len(to_treat_batch1) + len(to_treat_batch2) + len(to_treat_batch3)
    cut4 = len(to_treat_batch1) + len(to_treat_batch2) + len(to_treat_batch3) + len(to_treat_batch_drift)

    figure(figsize = (12, 8))
    suptitle('Mean peak height of Strontium signals on channel {}\nfor a cathode voltage of 4000'.format(channels_to_be_considered[0]) + r'$\pm$100 V')
    subplot(1, 2, 1)

    voltages_all = (array(voltages_all) / 1000.0).tolist()


    for i, main_filename in enumerate(to_treat):

        if main_filename in to_treat_batch1:
            if main_filename == to_treat_batch1[-1]:
                errorbar(array(voltages_all[:cut1]) + voltage_cathode, array(mean_peak_height_all[:cut1]), xerr = zeros_like(array(voltages_all[:cut1])) + voltage_uncertainty_anode, yerr = array(mean_peak_height_error_all[:cut1]), color = 'b', linestyle = '', marker = 'x', label = 'Total mean peak heights for uncontrolled pressure')
            else:
                errorbar(array(voltages_all[:cut1]) + voltage_cathode, array(mean_peak_height_all[:cut1]), xerr = zeros_like(array(voltages_all[:cut1])) + voltage_uncertainty_anode, yerr = array(mean_peak_height_error_all[:cut1]), color = 'b', linestyle = '', marker = 'x')

        elif main_filename in to_treat_batch2:
            if main_filename == to_treat_batch2[-1]:
                errorbar(array(voltages_all[len(to_treat_batch1):cut2]) + voltage_cathode, array(mean_peak_height_all[len(to_treat_batch1):cut2]), xerr = zeros_like(array(voltages_all[len(to_treat_batch1):cut2])) + voltage_uncertainty_anode, yerr = array(mean_peak_height_error_all[len(to_treat_batch1):cut2]), color = 'r', linestyle = '', marker = '.', label = 'Total mean peak heights for 1.100 bar')
            else:
                errorbar(array(voltages_all[len(to_treat_batch1):cut2]) + voltage_cathode, array(mean_peak_height_all[len(to_treat_batch1):cut2]), xerr = zeros_like(array(voltages_all[len(to_treat_batch1):cut2])) + voltage_uncertainty_anode, yerr = array(mean_peak_height_error_all[len(to_treat_batch1):cut2]), color = 'r', linestyle = '', marker = '.')

        elif main_filename in to_treat_batch3:
            if main_filename == to_treat_batch3[-1]:
                errorbar(array(voltages_all[cut2:cut3]) + voltage_cathode, array(mean_peak_height_all[cut2:cut3]), xerr = zeros_like(array(voltages_all[cut2:cut3])) + voltage_uncertainty_anode, yerr = array(mean_peak_height_error_all[cut2:cut3]), color = 'g', linestyle = '', marker = 's', label = 'Total mean peak heights for 1.140 bar')
            else:
                errorbar(array(voltages_all[cut2:cut3]) + voltage_cathode, array(mean_peak_height_all[cut2:cut3]), xerr = zeros_like(array(voltages_all[cut2:cut3])) + voltage_uncertainty_anode, yerr = array(mean_peak_height_error_all[cut2:cut3]), color = 'g', linestyle = '', marker = 's')

    # Setting up first fit
    fit1_xdata = []
    fit1_ydata = []
    fit1_xerror = []
    fit1_yerror = []
    for i, batch in enumerate(selection_fit1):
        if len(batch) == 0:
            continue
        for event in batch:
            name = selection_to_treat[i][event]
            index = to_treat.index(name)
            
            fit1_xdata.append(voltages_all[index] + voltage_cathode)
            fit1_ydata.append(mean_peak_height_all[index])
            fit1_xerror.append(voltage_uncertainty_anode)
            fit1_yerror.append(mean_peak_height_error_all[index])
    fit1_xdata = array(fit1_xdata)
    fit1_ydata = array(fit1_ydata)
    fit1_xerror = array(fit1_xerror)
    fit1_yerror = array(fit1_yerror)

    #fit1_parameters, fit1_parameters_error, fit1_reject, fit1_ssr, fit1_parameters_covariance = general_fit(fit1_xdata, fit1_ydata, pure_exponential_model, [-2.0/3.0, exp(2.0)], x_err = fit1_xerror, y_err = fit1_yerror, return_covariance = True)
    fit1_parameters, fit1_parameters_error, fit1_p_value, fit1_ssr, fit1_parameters_covariance = general_fit(fit1_xdata, fit1_ydata, exponential_ln_model, [0.1, 1.0], x_err = fit1_xerror, y_err = fit1_yerror, return_covariance = True)


    fit1_residuals = fit1_ydata - exponential_ln_model(fit1_parameters, fit1_xdata)

    plot(fit1_xdata, exponential_ln_model(fit1_parameters, fit1_xdata), 'k-', label = 'Exponential fit')
    #fit1_y_upper, fit1_y_lower = exponential_ln_confidence_band(fit1_xdata, fit1_ydata, fit1_parameters, fit1_parameters_covariance, fit1_ssr/(fit1_xdata.shape[0]-2.0), fit1_xdata.shape[0]-2.0)
    #plot(fit1_xdata, fit1_y_upper, 'r--', label = '95% confidence band')
    #plot(fit1_xdata, fit1_y_lower, 'r--')

    title('Peak heights')
    xlabel('Total applied voltage in kV')
    ylabel('Average total peak height in V')
    xscale('log')
    yscale('log')
    art = []
    lgd = legend(loc = 9, bbox_to_anchor=(0.5, -0.15))
    art.append(lgd)

    subplot(1, 2, 2)
    plot(array(voltages_all[:cut3]) + voltage_cathode, zeros_like(array(voltages_all[:cut3])), 'b-')
    plot(fit1_xdata, fit1_residuals, 'kx', linestyle = '', label = 'Exponential fit. SSR = {}'.format(fit1_ssr) + '\np-value: {}'.format(fit1_p_value))
    title('Residuals of fits')
    xlabel('Total applied voltage in kV')
    ylabel('Error of fits in V') 
    legend(loc = 9, bbox_to_anchor=(0.5, -0.15))

    savefig('mean_peak_heights.png', additional_artists=art, bbox_inches="tight")
    clf()

    
    figure(figsize = (12, 8))
    suptitle('Mean peak height of muon signals on channel {}\nfor an anode voltage of 3000'.format(channels_to_be_considered[0]) + r'$\pm$100 V')
    subplot(1, 2, 1)

    for i, main_filename in enumerate(to_treat):
        if main_filename in to_treat_batch_drift:
            if main_filename == to_treat_batch_drift[-1]:
                errorbar(array(voltages_all[cut3:cut4]) + voltage_anode, array(mean_peak_height_all[cut3:cut4]), xerr = zeros_like(array(voltages_all[cut3:cut4])) + voltage_uncertainty_anode, yerr = array(mean_peak_height_error_all[cut3:cut4]), color = 'k', linestyle = '', marker = 'x', label = 'Total mean peak heights for drift experiments')
            else:
                errorbar(array(voltages_all[cut3:cut4]) + voltage_anode, array(mean_peak_height_all[cut3:cut4]), xerr = zeros_like(array(voltages_all[cut3:cut4])) + voltage_uncertainty_anode, yerr = array(mean_peak_height_error_all[cut3:cut4]), color = 'k', linestyle = '', marker = 'x')

    # Setting up third fit
    fit3_xdata = []
    fit3_ydata = []
    fit3_xerror = []
    fit3_yerror = []
    for i, batch in enumerate(selection_fit3):
        if len(batch) == 0:
            continue
        for event in batch:
            name = selection_to_treat[i][event]
            index = to_treat.index(name)
            
            fit3_xdata.append(voltages_all[index] + voltage_anode)
            fit3_ydata.append(mean_peak_height_all[index])
            fit3_xerror.append(voltage_uncertainty_anode)
            fit3_yerror.append(mean_peak_height_error_all[index])
    fit3_xdata = array(fit3_xdata)
    fit3_ydata = array(fit3_ydata)
    fit3_xerror = array(fit3_xerror)
    fit3_yerror = array(fit3_yerror)

    fit3_parameters, fit3_parameters_error, fit3_p_value, fit3_ssr, fit3_parameters_covariance = general_fit(fit3_xdata, fit3_ydata, exponential_ln_model, [0.1, 1.0], x_err = fit3_xerror, y_err = fit3_yerror, return_covariance = True)

    fit3_residuals = fit3_ydata - exponential_ln_model(fit3_parameters, fit3_xdata)

    plot(fit3_xdata, exponential_ln_model(fit3_parameters, fit3_xdata), 'b-', label = 'Exponential fit')

    #fit3_y_upper, fit3_y_lower = exponential_ln_confidence_band(fit3_xdata, fit3_ydata, fit3_parameters, fit3_parameters_covariance, fit3_ssr/(fit3_xdata.shape[0]-2.0), fit3_xdata.shape[0]-2.0)
    #plot(fit3_xdata, fit3_y_upper, 'r--', label = '95% confidence band')
    #plot(fit3_xdata, fit3_y_lower, 'r--')

    title('Mean peak height')
    legend(loc = 9, bbox_to_anchor=(0.5, -0.15))
    xlabel('Total applied voltage in kV')
    ylabel('Average total peak height in V')
    yscale('log')
    xscale('log')

    subplot(1, 2, 2)
    plot(array(voltages_all[cut3:]) + voltage_anode, zeros_like(array(voltages_all[cut3:])), 'b-')
    plot(fit3_xdata, fit3_residuals, 'kx', linestyle = '', label = 'Exponential fit. SSR = {}'.format(fit3_ssr) + '\np-value: {}'.format(fit3_p_value))
    title('Residuals of fit')
    xlabel('Total applied voltage in kV')
    ylabel('Error of fits in V') 
    art = []
    lgd = legend(loc = 9, bbox_to_anchor=(0.5, -0.15))
    art.append(lgd)
    
    savefig('mean_peak_heights_drift.png', additional_artists=art, bbox_inches="tight")
    clf()
	
    fit_results_mean = open("_fit_results_mean", "w")
    fit_results_mean.write("The computed gas amplification constants\nfor the case where the mean of the peak heights is taken are:")
    fit_results_mean.write("For the first fit: " + str(fit1_parameters[0]) + "+-" + str(fit1_parameters_error[0]) + " SI\n")
    fit_results_mean.write("and: " + str(fit1_parameters[1]) + "+-" + str(fit1_parameters_error[1]) + " SI\n")
    fit_results_mean.write("The Townsend coefficient is therefore approximately: " + str(fit1_parameters[1]/wire_radius) + "+-" + str(fit1_parameters_error[1]/wire_radius) + " m^-1\n\n")
    fit_results_mean.write("For the third fit: " + str(fit3_parameters[0]) + "+-" + str(fit3_parameters_error[0]) + " SI\n")
    fit_results_mean.write("and: " + str(fit3_parameters[1]) + "+-" + str(fit3_parameters_error[1]) + " SI\n")
    fit_results_mean.write("The Townsend coefficient is therefore approximately: " + str(fit3_parameters[1]/wire_radius) + "+-" + str(fit3_parameters_error[1]/wire_radius) + " m^-1\n\n")
    fit_results_mean.close()
    
   
print "All data analyzed. For this measurements, the zero-level of the DRS-voltage indication as well as its uncertainty has be taken to be"
print "v_0 = ", voltage_zero_level, "+-", voltage_uncertainty, " V"
                           

