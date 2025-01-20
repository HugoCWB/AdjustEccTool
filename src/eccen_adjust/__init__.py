"""The eccentricity adjustment package.

This package contains code to adjust the eccentricity of a predicted pRF map to
match the eccentricity distribution implied by Horton and Hoyt (1991).
"""

def hh91_scale(surface_area, shape=0.75, min_eccen=0, max_eccen=90):
    '''Returns the `scale` parameter of Horton & Hoyt's magnification model.

    Returns the value of `scale` in the cortical magnification model of Horton
    and Hoyt (1991) that is appropriate for a V1 with the given surface-area
    (`area`), assuming that the maximum eccentricity of V1 is given by the
    optional parameter `max_eccen` (default: 90).

    Horton and Hoyt's original equation was:
      `m(r) = (scale / (shape + r))**2`
    where `r` is the eccentricity in degrees, `shape` is 0.75 degrees, and
    `scale` is 17.3 mm.
    '''
    shape_maxecc = shape + max_eccen
    shape_minecc = shape + min_eccen
    #den = np.pi * (np.log(shape_maxecc / shape) - max_eccen / shape_maxecc)
    den = np.pi * (
        min_eccen/shape_minecc 
        - max_eccen/shape_maxecc 
        + np.log(shape_maxecc/shape_minecc))
    return np.sqrt(surface_area / den)
def hh91_match(eccen, sarea, shape=0.75, min_eccen=0, max_eccen=90,
               max_steps=100, atol=1e-05, rtol=1e-08):
    """Returns eccentricity values that match the Horton and Hoyt (1991)
    distribution of eccentricity values.

    The optional parameter `shape` is the offset of the denominator in
    Horton & Hoyt's equation: `m(r) = (scale / (shape + r))**2` where `r`
    is the eccentricity in degrees, `shape` is 0.75 degrees, and scale is
    17.3 mm. The parameter `scale` is determined by the total surface area
    of V1 and thus is not required.

    The functions returns a vector `matched_eccen` that is similar to the
    argument `eccen`; in fact, `argsort(eccen)` must be equal to
    `argsort(matched_eccen)`. However, the values will be conformed to
    match the eccentricity distribution implied by Horton & Hoyt.
    """
    # Start by figuring out the scale parameter.
    total_sarea = np.sum(sarea)
    scale = hh91_scale(
        total_sarea,
        shape=shape, 
        min_eccen=min_eccen,
        max_eccen=max_eccen)
    # We can integrate H&H over the visual field to get the equation for the
    # amount of surface area of V1 from the foveal confluence up to a particular
    # eccentricity (r):
    #    f(r) = pi * scale^2 * (log((r + shape)/shape) - r/(r + shape)).
    # Note that this is derived by integrating Horton & Hoyt's function:
    #    m(r) = (scale / (shape + r))^2
    # over half of the visual field (for unilateral V1) via a double-integral:
    #    SS q m(q) dq dt; from q=0 to r, t=-pi/2 to pi/2 (t is the polar angle)
    #    (SS stands for the double-integral symbol.)
    # This yeilds the formula for f(r) above.
    # Suppose that the function g(a) is the inverse of f(r) (i.e., g(a) = r if
    # and only if f(r) = a). The function g(a) tells us what the maximum eccen
    # is for the most foveal `a` square-mm of V1. Unfortunately I wasn't able to
    # find a closed form for g(a) (but we can find it numerically).
    # If we sort all of the eccentricity values in V1 from the Benson14
    # retinotopic template then take the first (e.g.) 50 vertices, they will
    # represent the most foveal patch of V1; suppose their area is `a0`. If we
    # then add one more vertex to this patch, the area will be slightly larger;
    # we can call the new area `a` with `a > a0` and `a - a0` equal to the
    # surface area of the added vertex. The value `g(a)` represents the 
    # eccentricity radius of the new patch, so a good estimate of the new
    # vertex's eccentricity is `g((a + a0)/2)`. The argument here, `(a + a0)/2`
    # is easy to calculate; we do that here:
    ordering = np.argsort(eccen)
    sarea = sarea[ordering]
    a = np.cumsum(sarea)
    a0 = a - sarea
    arg = (a + a0)/2
    # We now just need to invert the function f. As mentioned above, this can
    # be done numerically. Given `f(r) = a`, we know `f(r) - a = 0`, so we can
    # find the zeros of the expression `f(r) - a`. This expression is
    # monotonically increasing for all non-negative values of `r`, so we can use
    # a very simple search-by-halving algorithm to find these.
    # First, we set up the function f(r); note that this version has been tweaked
    # slightly to handle the case where we limit the minimum eccentricity to
    # something other than 0, but otherwise matches the equations in the comments
    # above:
    c1 = np.pi * scale**2
    c2 = min_eccen + shape
    c3 = min_eccen / c2
    f = lambda r: c1 * (c3 + np.log((r + shape)/c2) - r/(r + shape))
    # We start by assuming that the solution is between the min and the max
    # eccentricities and each step we check the middle of the range and
    # eliminate the half of the range that cannot match. (This usually converges
    # quickly; probably in 20 steps or fewer depending on the tolerance.)
    rmin = np.full_like(arg, min_eccen)
    rmax = np.full_like(arg, max_eccen)
    for stepno in range(max_steps):
        r_estimate = (rmin + rmax) / 2
        arg_estimate = f(r_estimate)
        # If all of the values are within the error tolerance, we can end early.
        if np.isclose(arg, arg_estimate, atol=atol, rtol=rtol).all():
            print("Finished after", stepno, "steps.")
            break
        eps = arg_estimate - arg
        # Anywhere that f(r_estimate) is greater than arg, we know the max
        # possible value for r must be below the current value of the estimate.
        ii = eps > 0
        rmax[ii] = r_estimate[ii]
        # Conversely, anywhere that f(r_estimate) is less than arg, we know that
        # the min possible value of r is above the current estimate.
        ii = eps < 0
        rmin[ii] = r_estimate[ii]
    # Our best guess is the current mean of rmin and rmax!
    r = (rmin + rmax) / 2
    # We need to unsort these, however!
    r[ordering] = np.array(r)
    return r
