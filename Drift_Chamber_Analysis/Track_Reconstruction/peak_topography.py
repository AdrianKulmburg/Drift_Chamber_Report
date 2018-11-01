from numpy import *
from scipy.signal import find_peaks

def peak_topography(curve):
    peak_treshhold = 0.99

    min_size = min(curve) # This gives us the lowest point of the curve
    noise_estimate = min(curve[:50]) # This gives an extremely rough estimate
    # of the background noise. More is not needed though

    # In case the voltage due to background on the right is
    # higher than on the left, we need to take this into account for
    # the prominence.
    noise_estimate_right = curve[-1]

    prominence_approx = min(abs(min_size - noise_estimate), abs(min_size - noise_estimate_right))

    
    # This function assumes that the curve only has one peak below 0
    main_peak, _ = find_peaks(-curve, prominence = (prominence_approx, None), distance = 50)
    # The function has to be flipped, as find_peaks only searches for
    # maxima
    main_peak = main_peak[0]

    left_peak_height = abs(min_size - noise_estimate)
    right_peak_height = min(abs(min_size - noise_estimate_right), left_peak_height)

    #cut = where(curve<=-noise_treshhold)[0][0]
    cut = where(-curve[:main_peak] <= left_peak_height*(1-peak_treshhold) + abs(noise_estimate))[0]
    #left, _ = find_peaks(-curve[:cut], prominence = (0, 0.01), distance = 50)
    #if where(left < cut)[0].shape[0] != 0:
        #left = left[where(left < cut)[0][-1]]
    if cut.shape[0]!=0:
        left = cut[-1]
    else:
        left = -1

    #right, _ = find_peaks(-curve, prominence = (0, 1e-4), distance = 50)
    # Again, the prominence here has been determined empirically
    #if where(right > main_peak)[0].shape[0] != 0:
    #    right = right[where(right > main_peak)[0][0]]
    cut = where(-curve[main_peak:] <= right_peak_height*(1-peak_treshhold) + abs(noise_estimate_right))[0]
    if cut.shape[0]!=0:
        right = cut[0]
    else:
        right = -1
    # The determination of the right point is far more incacurate than for the left one

    return main_peak, left, right

    
