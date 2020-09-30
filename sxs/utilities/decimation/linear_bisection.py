"""Simple greedy algorithm from NINJA2 project"""

def minimal_grid(t, amp, phi, amp_tol=1e-5, phi_tol=1e-5):
    """Select a minimal grid from the input data allowing linear interpolation to given tolerances

    Storing the input data on the output grid will allow linear interpolation that
    agrees with the input data to within the given tolerances.

    This function is essentially a direct translation of the one originally written
    as part of the NINJA2 project <https://arxiv.org/abs/0709.0093>, which was
    provided as a utility named MinimizeGrid in the utilities/GridRefinement
    directory of the ninja2 svn repository.

    """
    import numpy as np

    # The objective here will be to create a vector of `bool`s, the same length as
    # t.  The truth value will correspond to whether or not that time step should
    # be kept in the final data.  We begin by assuming that the very first and last
    # steps should obviously be kept.  Then, there are two stages.  First is a
    # coarse stage, which steps through the data making intervals small enough to
    # reproduce the phi data at the interval's midpoint to within phi_tol, but no
    # smaller.  Second is the finer stage, which goes through each interval,
    # checking that every single point in the input data can be reproduced to
    # within phi_tol and amp_tol.  If that's not true, the interval is split evenly
    # into two, and the algorithm proceeds with the earlier interval.

    i0 = 0
    i1 = ((t.size-1) >> 1)  # = midpoint of the input data set
    n_points = 2
    include_indices = np.zeros(t.size(), dtype=bool)
    include_indices[0] = True
    include_indices[-1] = True

    # Returns True if the input (x,y) data can be interpolated between indices i0
    # and i1 to within a tolerance of tol at the midpoint between i0 and i1.
    def check(x, y, i0, i1, tol):
        if i1==i0+1:
            return True
        i_mid = ((i0+i1) >> 1)  # = floor(avg(i0+i1))
        if abs(y[i0] + (x[i_mid]-x[i0])*(y[i1]-y[i0])/(x[i1]-x[i0]) - y[i_mid]) > tol:
            return False  # Interval is not fine enough
        return True  # Interval is fine enough

    # Given data (t,phi), some initial index i0, and a guess for i1, this function
    # outputs the optimal index i1 so that (t,phi) can be interpolated between i0
    # and i1 to just within phi_tol of the full input data set.  Compare Numerical
    # Recipes's 'hunt' function; this is basically a hunt for that optimal index.
    def hunt(t, phi, phi_tol, i0, i1):
        inc = 1
        i1lo = i1
        i1hi = i1+1

        # Bracket the optimal i1 between i1lo and i1hi
        if check(t, phi, i0, i1lo, phi_tol):
            while check(t, phi, i0, i1hi, phi_tol):
                i1lo = i1hi;
                inc *= 2;
                i1hi += inc;
                if i1hi>int(t.size):
                    i1hi = t.size
                    break
        else:
            if i1lo<=i0+2:
                return i0+2
            i1hi = i1lo
            i1lo -= 1
            while not check(t, phi, i0, i1lo, phi_tol):
                i1hi = i1lo
                inc *= 2
                if inc > i1hi:
                    i1lo = i0+2
                    break
                else:
                    i1lo = i1hi-inc

        # Now use bisection between i1lo and i1hi
        while i1hi-i1lo != 1:
            i1m = ((i1hi+i1lo)>>1)
            if check(t, phi, i0, i1m, phi_tol):
                i1lo = i1m
            else:
                i1hi = i1m

        return i1lo

    # Coarse -- check only phi at midpoints of each interval This loop starts from
    #   the beginning of the data set, and forms the smallest interval such that
    #   the phi tolerance is achieved by linear interpolation.  Then, it moves to
    #   the end of that interval to find the next interval, etc.
    while ((i0+i1)>>1) < include_indices.size-1:
        # hunt for optimal i1
        i1 = hunt(t, phi, phi_tol, i0, i1)

        if not include_indices[i1]:
            include_indices[i1] = True
            n_points += 1
        i0 = i1
        i1 = 2*i1 - i0
        if i1<i0+2:
            i1 = i0+2

    # Fine -- check amp and phi at every point This loop goes through each of the
    #   intervals found above, and makes sure that every data point in both phi and
    #   amp can be reconstructed to within the given tolerances.  If not, it just
    #   adds the midpoint of the interval, and checks the new interval again.
    i0 = 0
    i1 = 1
    i = 1
    while i<t.size:
        if i>i1:  # This could happen below
            while i>i1:
                i1 += 1
            i0 = i
            while not include_indices[i0]:
                i0 -= 1
        while not include_indices[i1]:
            i1 += 1
        if i != i0 and i != i1:
            if (abs(1-(amp[i0]+(t[i]-t[i0])*(amp[i1]-amp[i0])/(t[i1]-t[i0]))/amp[i]) > amp_tol
                or abs(phi[i0]+(t[i]-t[i0])*(phi[i1]-phi[i0])/(t[i1]-t[i0])-phi[i]) > phi_tol):
                i1 = ((i0+i1)>>1)
                if not include_indices[i1]:
                    # if(i1==i0+1) { Caution(amp_tol, phi_tol); }
                    include_indices[i1] = True
                    n_points += 1
                continue
        if i==i1:
            i0 = i1
            i1 += 1
        i += 1

    return include_indices
