from numpy import *
from matplotlib.pyplot import *
from scipy.stats import *

from general_fit import *

def histogram_fit_model(data, bin_number, hist_label, models_to_fit, parameter_guesses, model_linestyle, model_color, model_marker, model_label, full_output = False):
    """ Plots the histogram of some data, computes and plots the mean
    and error for each bin, and fits it with some models, if given.
    Note that we assume there is no error on the height of the bins,
    i.e. if there is a systematic error in the data, it will have
    to be added at the very end of this whole computation.
    By default, the x-lim and y-lim are set to be data[0]-data[1]
    and 0-height of bins. The full_output makes this function return
    the height the of bins, their mean and error (excluding bins
    which contain only one data-point."""

    precision = 1000 # This is just for the plotting of the fits,
    # to know how precise they should be, i.e. on how many points
    # they should be evaluated.
    
    # Plotting histogram, and computing the centers and
    # standard deviations for each bin.
    bin_heights, bin_borders, _ = hist(data, bin_number, density = True, label = hist_label, alpha = 0.5)
    bin_heights_mean_center = []
    bin_heights_std_center = []

    max_length = bin_heights.shape[0]

    to_delete = []

    for i in xrange(max_length):
        if i != max_length - 1:
            data_values = data[where( (data >= bin_borders[i]) & (data < bin_borders[i+1]))[0]]
        else:
            data_values = data[where( (data >= bin_borders[i]) & (data <= bin_borders[i+1]))[0]]
        
        if data_values.shape[0] <= 1:
            print ""
            print "/!\\ Care! NaN produced while computing standard error."
            print "Throwing the measurement away. This is problematic if this"
            print "happens too often. In this case, you might want to reduce"
            print "the amount of bins."
            print ""
            to_delete.append(i)
        else:
            bin_heights_std_center.append(sem(data_values))
            bin_heights_mean_center.append(data_values.mean())
    bin_heights = delete(bin_heights, to_delete)
    bin_heights_mean_center = array(bin_heights_mean_center)
    bin_heights_std_center = array(bin_heights_std_center)

    # Plotting the means, as well as their error
    errorbar(bin_heights_mean_center, bin_heights, xerr = bin_heights_std_center, color = 'red', marker = 'x', linestyle = '', label = "Mean of each bin")

    # Doing a fit for each model
    parameters_ideal = []
    parameters_errors = []
    rejectings = []
    SSRs = []
    for i, model_to_fit in enumerate(models_to_fit):
        # Doing the fitting
        parameter_ideal, parameter_error, reject, SSR = general_fit(bin_heights_mean_center, bin_heights, model_to_fit, parameter_guesses[i], x_err = bin_heights_std_center)
        
        # Preparing to plot the fit
        fit_function = lambda x : model_to_fit(parameter_ideal, x)
        start = min(data)
        end = max(data)
        x_range = linspace(start, end, precision)
        plot(x_range, fit_function(x_range), linestyle = model_linestyle[i], color = model_color[i], marker = model_marker[i], label = model_label[i])
        xlim(start, end)
        ylim(0.0, 1.1*max(bin_heights))
        parameters_ideal.append(parameter_ideal)
        parameters_errors.append(parameter_error)
        rejectings.append(reject)
        SSRs.append(SSR)
    if full_output:
        return parameters_ideal, parameters_errors, rejectings, SSRs, bin_heights, bin_heights_mean_center, bin_heights_std_center
    else:
        return parameters_ideal, parameters_errors, rejectings, SSRs
        






