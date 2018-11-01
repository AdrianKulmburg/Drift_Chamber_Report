from numpy import *
from matplotlib.pyplot import *

from pdf_drift import pdf_drift_zone

D = 98.5
L = 104.0
anode_padding = 4.2


chn1_zone = [30., 52.]
#chn2_zone = [20., 30.]
#chn3_zone = [10., 20.]
chn4_zone = [0.0, 10.]
#chn5_zone = [-10., 0.0]
chn6_zone = [-20., -10.]
#chn7_zone = [-30., -20.]
chn8_zone = [-52., -30.]

#zones = [chn1_zone, chn2_zone, chn3_zone, chn4_zone, chn5_zone, chn6_zone, chn7_zone, chn8_zone]
zones = [chn1_zone, chn4_zone, chn6_zone, chn8_zone]

channels_to_be_considered = [1, 4, 6, 8]

x = linspace(-D/2.0, D/2.0, 1000)

for i, chn in enumerate(channels_to_be_considered):
    pdf_model = pdf_drift_zone(zones[i][0], zones[i][1], D, L, anode_padding)
    plot(x, pdf_model(x), label = 'Channel zone of channel {}'.format(chn) + ', from ' + str(zones[i][0]) + 'mm to ' + str(zones[i][1]) + 'mm')
    
title("Theoretical pdfs for the measured signal starting times.\nD = 98.5 mm, L = 104.0 mm and anode padding is 4.2 mm.")
xlabel("Horizontal distance to center of drift chamber, in mm")
ylabel("Vertical distance to center of drift chamber, in mm")
art = []
lgd = legend(loc = 9, bbox_to_anchor=(0.5, -0.15))
art.append(lgd)
savefig('models.png', additional_artists=art, bbox_inches="tight")
    
