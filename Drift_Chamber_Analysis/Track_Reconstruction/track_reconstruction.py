from numpy import *
from matplotlib.pyplot import *
from drift_velocities_light import *

from general_fit import *

def linear_model(parameters, data):
    a = parameters[0]
    b = parameters[1]
    return a*data + b

event_data_folder = '../Drift_Velocities/Drift_3800_1000_Data/'
event_main_filename = 'Drift_3800_1000'
event_number = 54
voltage = 3800

channels_to_be_considered = [1, 2, 3, 4, 7, 8]

time_uncertainty = 15e-12 # In m
time_uncertainty = time_uncertainty * 1e9 # Needs to be in ns

# All distances here are in mm
D = 98.5
L = 104.0
anode_padding = 4.2

chn1_zone = [30., 52.]
chn2_zone = [20., 30.]
chn3_zone = [10., 20.]
chn4_zone = [0.0, 10.]
chn5_zone = [-10., 0.0]
chn6_zone = [-20., -10.]
chn7_zone = [-30., -20.]
chn8_zone = [-52., -30.]

chn1_center = 35.0
chn2_center = 25.0
chn3_center = 15.0
chn4_center = 5.0
chn5_center = -5.0
chn6_center = -15.0
chn7_center = -25.0
chn8_center = -35.0

chn_zones = [chn1_zone, chn2_zone, chn3_zone, chn4_zone, chn5_zone, chn6_zone, chn7_zone, chn8_zone]
chn_centers = [chn1_center, chn2_center, chn3_center, chn4_center, chn5_center, chn6_center, chn7_center, chn8_center]

peaks_correct, event_start_times, drift_velocity, drift_velocity_error, t_Right, t_Right_error = drift_velocities(event_main_filename, event_number)

if peaks_correct == False:
    print "Something went wrong. The peaks were found to be incorrect. Closing program."
    exit()

x = []
y = []
x_error = []
y_error_upper = []
y_error_lower = []
y_error_min = []
for i, chn in enumerate(channels_to_be_considered):
    y.append(chn_centers[chn-1])
    y_error_upper.append(chn_zones[chn-1][1] - chn_centers[chn-1])
    y_error_lower.append(chn_centers[chn-1] - chn_zones[chn-1][0])
    y_error_min.append(min(y_error_upper[-1], y_error_lower[-1]))

    x.append(anode_padding + D - (t_Right - event_start_times[i]) * drift_velocity)
    x_error.append(sqrt(t_Right_error**2*drift_velocity**2 + time_uncertainty**2*drift_velocity**2 + drift_velocity_error**2 * (t_Right - event_start_times[i])**2))

x = array(x) - D/2.0
y = array(y)
x_error = array(x_error)
y_error = vstack((y_error_lower, y_error_upper))

anodes_y = array(chn_centers)
anodes_x = zeros_like(anodes_y) + anode_padding - D/2.0

ground_plate_y = linspace(-L/2.0, L/2.0, 10)
ground_plate_x = zeros_like(ground_plate_y) - D/2.0

cathode_y = linspace(-L/2.0, L/2.0, 10)
cathode_x = zeros_like(cathode_y) + D/2.0

odr_parameter_ideal, odr_parameter_error, odr_p_value, odr_SSR = general_fit(x, y, linear_model, [0.0, 0.0], x_err = x_error, y_err = array(y_error_min))

X = linspace(min(x)-10.0, max(x)+10.0, 100)


plot(ground_plate_x, ground_plate_y, color = 'k', marker = '', linestyle = '-', linewidth = 5.0, label = 'Ground plate')
plot(cathode_x, cathode_y, color = 'g', marker = '', linestyle = '-', linewidth = 5.0, label = 'Cathode plate')
errorbar(x, y, xerr = x_error, yerr = y_error, marker = 'x', color = 'r', linestyle = '', label = 'Impact position of the muon')
plot(X, linear_model(odr_parameter_ideal, X), marker = '', color = 'k', linestyle = '--', label = 'Linear fit, p-value: {}'.format(odr_p_value))
plot(anodes_x, anodes_y, marker = 'o', color = 'b', linestyle = '', label = 'Anode positions')
xlabel(r'$x$-position in the measuring zone, in mm')
ylabel(r'$y$-position in the measuring zone, in mm')
title('Track reconstruction for a cathode voltage of {}'.format(voltage) + r'$\pm$100V.' + '\nEvent number: ' + str(event_number))
ylim(-L/2.0 - 5.0, L/2.0 + 5.0)

art = []
lgd = legend(loc = 9, bbox_to_anchor=(0.5, -0.15))
art.append(lgd)
savefig('track_reconstruction_voltage_' + str(voltage) + '_event_' + str(event_number) + '.png',  additional_artists=art, bbox_inches="tight")







