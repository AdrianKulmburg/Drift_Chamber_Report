from numpy import *

# This program computes an approximation of the standard error
# due to the noise of the DRS. To do so, it takes one
# predetermined file that has been deemed reasonable
# (i.e. without any signal), reads all the different
# channels out of it, takes the mean of all the values, and then
# computes the maximum deviation from the mean (i.e. not the
# standard error, but indeed the maximal deviation).
# That way we are sure to get a realistic estimate for the error
# that one might get for the voltage measurements.

reference_file = "../Drift_Velocities/Drift_4000_1000_Data/"
reference_file_main_filename = "Drift_4000_1000"
reference_event_number = 12

channels_to_consider = [1, 2, 3, 4, 7, 8]

def find_noise_level():
    voltages_all = array([])
    for chn in channels_to_consider:
        voltage = loadtxt(reference_file + reference_file_main_filename + "_chn{}_v".format(chn))
        voltage = voltage[reference_event_number - 1]
        voltages_all = concatenate((voltages_all, voltage))
    return voltages_all.mean(), max(abs(voltages_all - voltages_all.mean()))

