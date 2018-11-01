from numpy import *
from matplotlib.pyplot import *

event_data_folder = '../Drift_Velocities/Drift_4000_1000_Data/'
event_main_filename = 'Drift_4000_1000'
event_number = 12

channels_to_be_considered = [7]

for chn in channels_to_be_considered:
    voltages = loadtxt(event_data_folder + event_main_filename + "_chn{}_v".format(chn))
    times = loadtxt(event_data_folder + event_main_filename + "_chn{}_t".format(chn))

    plot(times[event_number-1], voltages[event_number-1], marker = '', linestyle = '-', color = 'k')

title('Signal observed by channel 7\nfor event 12 in Drift_4000_1000_Data')
xlabel('Time in ns')
ylim((-0.1, 0.1))
ylabel('Signal height in V')
savefig('specific_noise.png')


