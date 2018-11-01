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

channels_to_be_considered = [1, 2, 3, 4, 7, 8]

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
D = 98.5
L = 104.0
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

def drift_velocities(main_filename, event_number):
    print "Starting to compute the starting times for event ", event_number, " of " + main_filename + "..."
    sys.stdout.flush()
    data_path = "../Drift_Velocities/" + main_filename + "_Data/"
    
    event_start_times = []

    voltages = []
    times = []
    for chn in channels_to_be_considered:
        print "o Starting to read in channel {}...".format(chn),
        sys.stdout.flush()
        voltages.append(loadtxt(data_path + main_filename + "_chn{}_v".format(chn)))
        times.append(loadtxt(data_path + main_filename + "_chn{}_t".format(chn)))
        print "Done reading."
        sys.stdout.flush()

    print "o-> Starting to check the peaks...",
    sys.stdout.flush()

    for j, chn in enumerate(channels_to_be_considered):
        curve = voltages[j][event_number-1]
        time = times[j][event_number-1]
        if curve_is_fine(curve) == False:
            all_curves_are_fine = False
            return False, [], 0, 0, 0, 0
        main, left, right = peak_topography(curve)
        if left == -1:
            all_curves_are_fine = False
            return False, [], 0, 0, 0, 0
        event_start_times.append(time[left])


    print "Done."
    print "Now computing the drift velocities and t_0 for that cathode voltage..."
    sys.stdout.flush()

    drift_velocities = []
    drift_velocities_errors = []
    t_0s = []
    t_0s_errors = []
  
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
                sys.stdout.flush()
                strange_occurences += 1
                all_curves_are_fine = False
                break
            start_times.append(time[left])
        if all_curves_are_fine == False:
            continue

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
        
        parameters_ideal, parameters_errors, p_values, SSRs, bin_heights, bin_heights_mean_center, bin_heights_std_center = histogram_fit_model(array(start_times), 10, "Relative occurences for signal starting time", models_to_fit, parameter_guesses, model_linestyle, model_color, model_marker, model_label, full_output = True)
        
        clf()

        systematic_error_velocity = sqrt(2)*parameters_ideal[0][1]**2*time_uncertainty/D
        drift_velocities.append(parameters_ideal[0][1])
        drift_velocities_errors.append(parameters_errors[0][1] + systematic_error_velocity)
        t_0s.append(parameters_ideal[0][0])
        t_0s_errors.append(parameters_errors[0][0])

    mean_drift_velocity = sum(array(drift_velocities)/array(drift_velocities_errors)**2)/sum(1.0/array(drift_velocities_errors)**2)
    mean_drift_velocity_error = sqrt(1.0/sum(1.0/array(drift_velocities_errors)**2)) + sqrt(2)*mean_drift_velocity**2*time_uncertainty/D

    t_Rights = array(t_0s) + D/(2.0 * array(drift_velocities))
    t_Rights_errors = sqrt(array(t_0s_errors)**2 + D**2*array(drift_velocities_errors)**2/array(drift_velocities)**4)

    t_Right = sum(t_Rights/t_Rights_errors**2)/sum(1.0/t_Rights_errors**2)
    t_Right_error = sqrt(1.0/sum(1.0/t_Rights_errors**2) + time_uncertainty)

    print "Done."
    print ""
    sys.stdout.flush()


    print "ooo Done with computing drift velocities."
    print ""
    print ""
    print ""
    sys.stdout.flush()
    return True, event_start_times, mean_drift_velocity, mean_drift_velocity_error, t_Right, t_Right_error

    
