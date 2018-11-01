from numpy import *
from scipy.signal import find_peaks

# This is a modified version of the curve_is_fine routine used for the
# drift velocities. Unlike its predecessor, it also returns the height of
# the peak in case there is one, and returns the approximated level of the
# noise in case no peak was found (not unique peak, just if no peak was found
# at all, i.e. the signal is 'zero').

def curve_is_fine(curve):
    # Testing sizes
    #if min(curve) > noise_treshhold:
    #    return False
    if min(curve) < -0.499:
        return False, -1

    # Now we know that the curve has the right size.

    # Now, we check if it has the right shape; for that,
    # it is enough to check that it has exactly one peak.

    min_size = min(curve) # This gives us the lowest point of the curve
    noise_estimate = min(curve[:50]) # This gives an extremely rough estimate
    # of the background noise. More is not needed though

    # The following few lines had to be added, since otherwise null-signals
    # would not make it through
    if abs(min_size) <= 0.02:
        return True, 0

    # In case the voltage due to background on the right is
    # higher than on the left, we need to take this into account for
    # the prominence.
    noise_estimate_right = curve[-1]

    prominence_approx = min(abs(min_size - noise_estimate), abs(min_size - noise_estimate_right))
    
    # Now, if we search for peaks with this prominence, we'll find
    # just the highest peak. But we want to be sure that they are no
    # other significant peaks. For this, we search for all peaks that
    # have a significant prominance compared to the main peak.
    # This is a bit arbitrary, but a secondary peak will have a prominence
    # of about 1/10 of the main peak. The rest will be noise, whose
    # prominence is far lower.
    significance = 10.0
    peaks, _ = find_peaks(-curve, prominence = (prominence_approx/significance, None), distance = 50)
    # The function has to be flipped, as find_peaks only
    # searches for local maxima
    if peaks.shape[0] > 1:
        return False, -1
    elif peaks.shape[0] == 0:
        return True, 0
    elif peaks.shape[0] == 1:
        return True, abs(min_size) - abs(noise_estimate)
    return False, -1
