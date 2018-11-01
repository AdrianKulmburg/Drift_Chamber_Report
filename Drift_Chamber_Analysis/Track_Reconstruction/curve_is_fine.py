from numpy import *
from scipy.signal import find_peaks

def curve_is_fine(curve):
    # Testing sizes
    #if min(curve) > noise_treshhold:
    #    return False
    if min(curve) < -0.499:
        return False

    # Now we know that the curve has the right size.

    # Now, we check if it has the right shape; for that,
    # it is enough to check that it has exactly one peak.

    min_size = min(curve) # This gives us the lowest point of the curve
    noise_estimate = min(curve[:50]) # This gives an extremely rough estimate
    # of the background noise. More is not needed though

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
    if peaks.shape[0] != 1:
        return False
    return True
