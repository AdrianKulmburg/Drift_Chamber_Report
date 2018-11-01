from numpy import *
from curve_is_fine import *
from peak_topography import *
from pdf_drift import *
from general_fit import *
from histogram_fit_model import *
from matplotlib.pyplot import *
from sys import argv
import sys
import os

#########
# Height still has to be checked!!!!!!


trigger_time = 1645 # In ns.
noise_treshhold = 0.02
channels_to_be_considered = [1, 2, 3, 4, 7, 8]
strange_occurences = 0

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

voltage_uncertainty = 100.0 # This is the uncertainty for the
# cathode voltage, in V

# Dimensions of the drift chamber
D = 98.5 #105 #94.3 # Width (i.e. Cathode to Anode)
L = 104.0 #245 # Height
anode_padding = 4.2
# We assume there is no error on the distance measurements,
# except possibly for the height measurement.

# Vertical 'height' on which electrons would go to one particular anode.
# Corresponds to each channel.
chn1_zone = [30., 52.]
chn2_zone = [20., 30.]
chn3_zone = [10., 20.]
chn4_zone = [0.0, 10.]
chn5_zone = [-10., 0.0]
chn6_zone = [-20., -10.]
chn7_zone = [-30., -20.]
chn8_zone = [-52., -30.]

chn_zones = [chn1_zone, chn2_zone, chn3_zone, chn4_zone, chn5_zone, chn6_zone, chn7_zone, chn8_zone]

to_treat = ["Drift_2000_1000", "Drift_2500_1000", "Drift_3000_1000",
            "Drift_3500_1000", "Drift_3800_1000", "Drift_4000_1000"]
#to_treat = ["Drift_1000_1000", "Drift_2000_1000", "Drift_2500_1000", "Drift_3000_1000",
#            "Drift_3500_1000", "Drift_3800_1000", "Drift_4000_1000"]

if len(argv) == 1:
    print "No arguments given. You have to use it as follows:"
    print "python drift_velocities.py base_name_directory1 base_name_directory2 ..."
    print " or"
    print "python drift_velocities.py -all"
    exit()

if argv[1] != "-all":
    to_treat = argv[1:]
    flag_all = False
else:
    flag_all = True
    mean_velocities_all = []
    mean_of_mean_velocities_all = []
    mean_velocities_error_all = []
    mean_of_mean_velocities_error_all = []
    voltages_all = []

for main_filename in to_treat:
    _ = os.system('clear')
    print "Starting now with " + main_filename
    sys.stdout.flush()
    data_path = "../Drift_Velocities/" + main_filename + "_Data/"
    histogram_path = data_path + main_filename + "_Histograms_and_Crosschecks/"
    if not os.path.exists(histogram_path):
        os.makedirs(histogram_path)
    results = open(histogram_path + main_filename + "_start_times", "w")
    parameters_file = open(histogram_path + main_filename + "_parameters", "w")
    mean_velocities = []
    upper_velocities = []
    lower_velocities = []
    mean_velocities_errors = []
    upper_velocities_errors = []
    lower_velocities_errors = []
    voltages = []
    times = []

    for chn in channels_to_be_considered:
        print "o Starting to read in channel {}...".format(chn),
        sys.stdout.flush()
        voltages.append(loadtxt(data_path + main_filename + "_chn{}_v".format(chn)))
        times.append(loadtxt(data_path + main_filename + "_chn{}_t".format(chn)))
        print "Done reading."
        sys.stdout.flush()
    number_of_events = voltages[0].shape[0]


    print "o-> Starting to check the peaks...",
    sys.stdout.flush()

    good_start_times_all_channels = []
    for chn in channels_to_be_considered:
        good_start_times_all_channels.append([])

    for i in xrange(number_of_events):
        all_curves_are_fine = True
        start_times = []
        for j, chn in enumerate(channels_to_be_considered):
            curve = voltages[j][i]
            time = times[j][i]
            if curve_is_fine(curve) == False:
                all_curves_are_fine = False
                break
            main, left, right = peak_topography(curve)
            if left == -1:
                all_curves_are_fine = False
                break
            if left < 50:
                print ""
                print "o<--> Something strange happened here: chn{}".format(chn)
                print "o<--> Entry {}".format(i)
                print "o<--> The starting time was evaluated to be <50 units."
                print "o<--> If this happens too often, you should be alerted."
                print "o<--> Otherwise, you can safely ignore this message"
                print "o<--> and discard the measurement in question."
                print "o<--> Continuing to read channel {}...".format(chn),
                sys.stdout.flush()
                strange_occurences += 1
                all_curves_are_fine = False
                break
            start_times.append(time[left])
            ## For Debugging
            #if left < 50:
            #    plot(time, curve)
            #    plot(time[left], curve[left], "x")
            #    plot(time[main], curve[main], "gx")
            #    show()
        if all_curves_are_fine == False:
            continue

        results.write("For event " + str(i+1) + " the starting times were: ")
        for entry in start_times:
            results.write(str(entry)+" ")
        results.write("\n")

        for j, chn in enumerate(channels_to_be_considered):
            good_start_times_all_channels[j].append(start_times[j])



    print "Done."
    print "o-> Processing and saving data...",
    sys.stdout.flush()

    for i, chn in enumerate(channels_to_be_considered):

        start_times = good_start_times_all_channels[i] # Yup, as you may have noticed, I needed to rewrite the code preceding this line

        min_zone = chn_zones[chn-1][0]
        max_zone = chn_zones[chn-1][1]
        
        models_to_fit = [pdf_drift_model(min_zone, max_zone, D, L, anode_padding)]
        parameter_guesses = [[3000.0/2.0, D/3000.0]]
        model_linestyle = ['--']
        model_marker = ['']
        model_label = ["Fit for channel zone"]
        model_color = ['k']

        figure(figsize = (11, 4.5))
        suptitle("Peak-starting times for channel {} with a Cathode voltage of ".format(chn) + main_filename[6:10] + r"$\pm$100V")
        subplot(1, 2, 1)
        parameters_ideal, parameters_errors, p_values, SSRs, bin_heights, bin_heights_mean_center, bin_heights_std_center = histogram_fit_model(array(start_times), 10, "Relative occurences for signal starting time", models_to_fit, parameter_guesses, model_linestyle, model_color, model_marker, model_label, full_output = True)
        
        legend(loc = 9, bbox_to_anchor=(0.5, -0.15))
        xlabel("Starting time of peaks in ns")
        title("Distribution of starting times")

        subplot(1, 2, 2)
        plot(array(start_times), zeros_like(array(start_times)), 'b', linestyle = '-', marker = '')
        plot(bin_heights_mean_center, bin_heights - pdf_drift_model(min_zone, max_zone, D, L, anode_padding)(parameters_ideal[0], bin_heights_mean_center), 'rx', linestyle = '', label = 'Fit Error. SSR = {}\n'.format(SSRs[0]) + 'p-value : {}'.format(p_values[0]))
        xlabel('Starting times in ns (in this case, center of bins)')
        title('Residuals of fit')

        art = []
        lgd = legend(loc = 9, bbox_to_anchor=(0.5, -0.15))
        art.append(lgd)
        savefig(histogram_path + main_filename + '_chn{}_starting_times.png'.format(chn), additional_artists=art,
    bbox_inches="tight")
        clf()
        

        parameters_file.write("_________________________________________________\n")
        parameters_file.write("Parameters for channel {}:\n".format(chn))
        parameters_file.write("o-> Velocity: " + str(parameters_ideal[0][1]) + "+-" + str(parameters_errors[0][1]) + "mm/ns\n")
        parameters_file.write("o-> Displacement: " + str(parameters_ideal[0][0]) + "+-" + str(parameters_errors[0][0]) + "ns\n")
        parameters_file.write("o-> p-value: " + str(p_values[0]) + "\n")
        parameters_file.write("o-> Sum of squared residuals: " + str(SSRs[0]) + "\n")
        parameters_file.write("\n")
        parameters_file.write("\n")
      
        if p_values[0] < 0.90:
            strange_occurences +=1
 
        mean_velocities.append(parameters_ideal[0][1])
        mean_velocities_errors.append(parameters_errors[0][1])

    print "Done."
    print ""
    sys.stdout.flush()


    mean_of_mean_velocities = sum(array(mean_velocities)/array(mean_velocities_errors)**2)/sum(1.0/array(mean_velocities_errors)**2)
    mean_of_mean_velocities_error = sqrt(1.0/sum(1.0/array(mean_velocities_errors)**2))

    systematic_error_velocity = sqrt(2)*mean_of_mean_velocities**2*time_uncertainty/D
    parameters_file.write("Taking the weighted mean of all velocities for the mean channel yields:\n")
    parameters_file.write("Velocity = " + str(mean_of_mean_velocities) + "+-" + str(mean_of_mean_velocities_error + systematic_error_velocity) + "mm/ns\n")
    parameters_file.write("or\n")
    parameters_file.write("Velocity = " + str(mean_of_mean_velocities*1e6) + "+-" + str((mean_of_mean_velocities_error + systematic_error_velocity)*1e6) + "m/s\n\n")
    parameters_file.write("The error due to systematic errors is equal to " + str(systematic_error_velocity) + "mm/ns, or " + str(systematic_error_velocity*1e6) + "m/s.\n")
    parameters_file.write("EOF")

    parameters_file.close()
    results.close()

    errorbar(array(channels_to_be_considered), array(mean_velocities), yerr = array(mean_velocities_errors), color = 'b', linestyle = '', marker = 'x', label = "Mean velocity", capsize = 10)
    errorbar(array(channels_to_be_considered), zeros_like(channels_to_be_considered)+mean_of_mean_velocities, yerr = zeros_like(channels_to_be_considered) + mean_of_mean_velocities_error + systematic_error_velocity, color='k', linestyle = '-', marker = '.', label = "Estimated velocity")
    
    title("Overview of velocity estimates for a Cathode voltage of " + main_filename[6:10] + r"$\pm$100 V")
    legend(loc = 'best')
    xlabel("Channels")
    ylabel(r"Velocity in mm$\cdot$ns$^{-1}$")
    savefig(histogram_path + main_filename + '_Velocity_Overview.png')
    clf()

    print "ooo Done with " + main_filename + ". ooo"
    print ""
    print ""
    print ""
    sys.stdout.flush()

    if flag_all:
        mean_velocities_all.append(array(mean_velocities))
        mean_velocities_error_all.append(array(mean_velocities_errors) + systematic_error_velocity)
        mean_of_mean_velocities_all.append(mean_of_mean_velocities)
        mean_of_mean_velocities_error_all.append(mean_of_mean_velocities_error + systematic_error_velocity)
        voltages_all.append(float(main_filename[6:10]))

if flag_all:
    voltages_all = array(voltages_all)
    errorbar(voltages_all, array(mean_of_mean_velocities_all)*1e6, xerr = zeros_like(voltages_all)+voltage_uncertainty, yerr = array(mean_of_mean_velocities_error_all)*1e6, color = 'r', linestyle = '', marker = 'o', label = "Mean velocity over all channels")
    mean_velocities_all = array(mean_velocities_all)*1e6
    mean_velocities_error_all = array(mean_velocities_error_all)*1e6
    for i in xrange(len(channels_to_be_considered)-1):
        errorbar(voltages_all, mean_velocities_all[:, i], xerr = zeros_like(voltages_all)+voltage_uncertainty, yerr = mean_velocities_error_all[:, i], color = 'k', linestyle = '', marker = '.', alpha = 0.2)
    errorbar(voltages_all, mean_velocities_all[:,-1], xerr = zeros_like(voltages_all)+voltage_uncertainty, yerr = mean_velocities_error_all[:,-1], color = 'k', linestyle = '', marker = '.', label = "Velocity for each channel", alpha = 0.2)
    xlabel("Cathode voltage in V")
    ylabel(r"Velocity in m$\cdot$s$^{-1}$")
    title(r"Drift velocity $v_D$ for different cathode voltages" + "\n" + "and an anode voltage of $3000\pm100$ V")
    legend(loc='upper center')
    
    savefig("drift_velocities.png")
    clf()

_ = os.system('clear')
print "All data analyzed."
print "Number of strange occurences: ", strange_occurences
                           

