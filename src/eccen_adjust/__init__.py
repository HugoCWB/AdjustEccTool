"""The eccentricity adjustment package.

This package contains code to adjust the eccentricity of a predicted pRF map to
match the eccentricity distribution implied by Horton and Hoyt (1991).
"""
import numpy as np

# def load_v1(sub,hemlabel,fname):
#     # fname = str(benson_path) + '/' + hemlabel + '-adjusted-eccen-VtxLabels.txt'

#     hem = sub.hemis[hemlabel]
#     print(hem)
#     # If we want to run this on another subject, we can apply the the Benson14
#     # template using this code::
#     # ret = ny.vision.predict_retinotopy(hem, sym_angle=False)
#     # hem = hem.with_prop(
#     #     b14_polar_angle=ret['angle'],
#     #     b14_eccentricity=ret['eccen'],
#     #     b14_visual_area=ret['varea'])

#     # Pick out just V1 on the white surface, using Benson atlas
#     v1 = hem.white_surface.submesh(hem.mask(('b14_visual_area', 1)))
#     # Display the size of V1 mask
#     print('Size V1 mask=',len(hem.mask(('b14_visual_area', 1))),sep='') 
#     # Display the information about V1 mesh, such as number of faces and vertices
#     print('V1=',v1,sep='')
#     # The polar angle values are fine we assume; we just want to rescale
#     # the eccentricity.
#     # print(dir(hem.white_surface)) 
#     # help(v1)
#     np.set_printoptions(threshold=np.inf)
#     # Save the V1 vertices index into a txt file
#     ny.save(fname, v1.labels)
#     return hem,v1

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

def adjust_eccen_in_v1(v1_eccen, v1_surface_area, shape=0.75,min_eccen=0,max_eccen=90):
   scale = hh91_scale(
       np.sum(v1_surface_area),
       shape=shape,
       min_eccen=min_eccen,
       max_eccen=max_eccen)
   # Run the matching algorithm and report on it:
   r = hh91_match(
       v1_eccen, v1_surface_area,
       shape=shape, min_eccen=min_eccen, max_eccen=max_eccen)
   # Return the adjusted eccen:
   return r,scale

# def adjust_eccen(v1,hemlabel,fname,shape=0.75,min_eccen=0,max_eccen=90):
#     # Go ahead and adjust the eccentricity using the function above.
#     # We can pick a shape parameter; Horton & Hoyt found it to be 0.75 on
#     # average (and we found that to be correct on average), but there is also
#     # a fair amount of variance in the population w.r.t. this parameter,
#     # anecdotally. We also assume that the visual field goes out to 90 degrees 
#     # (technically not correct, but works well enough for the model). Finally, 
#     # we assume that V1 starts at 0° of eccentricity.

#     # Our initial eccentricity values are the Benson 14 template predictions.
#     r0 = v1.prop('b14_eccentricity')
#     sarea = v1.prop('midgray_surface_area')
#     # We can calculae the scale using the surface area, shape, and max-eccen:
#     r,scale = adjust_eccen_in_v1(r0, sarea, shape=shape, min_eccen=min_eccen, max_eccen=max_eccen)
#     # Save adjusted eccentricities
#     ny.save(fname, r)
#     print(hemlabel, ' - Benson14 Eccentricity: min=', np.min(r0), ', max=', np.max(r0), sep='')
#     print(hemlabel, ' - Adjusted Eccentricity: min=', np.min(r), ', max=', np.max(r), sep='') 
#     # Return the initial, adjusted eccentricities and scale variable
#     return r0,r,scale

def ring_area_deg2(min_eccen, max_eccen, hemifield=True):
    """Computes the area (in square degrees) of a ring in the visual field."""
    if hemifield:
        return np.pi/2 * (max_eccen**2 - min_eccen**2)
    else:
        return np.pi * (max_eccen**2 - min_eccen**2)
def cmag(eccen, sarea, hwidth=0.075):
    n = len(eccen)
    if not isinstance(hwidth, int):
        hwidth = int(np.floor(n*hwidth))
    width = 2*hwidth + 1
    ii = np.argsort(eccen)
    eccen = eccen[ii]
    sarea = sarea[ii]
    ii_min = np.arange(0, n - width)
    ii_max = np.arange(width, n)
    out_ecc = eccen[hwidth:n-hwidth-1]
    cumsarea = np.concatenate([np.cumsum(sarea), [0]])
    rings_srfarea = cumsarea[ii_max] - cumsarea[ii_min - 1]
    rings_visarea = ring_area_deg2(eccen[ii_min], eccen[ii_max])
    return (out_ecc, rings_srfarea / rings_visarea)

# Plotting Function ############################################################

def plot_originalvsadjusted(r0, r, fname=None):
    # fname = str(Save_DIR / subs) + '_' + hemlabel + '_BensonVsAdjusted.png'

    # What did the above cell do? Let's plot it:
    (fig,ax) = plt.subplots(1,1, figsize=(3,3), dpi=128)
    ax.loglog(r0, r, 'k.')
    ax.loglog([0.1,90], [0.1,90], 'r:', zorder=-1)
    ax.set_xlabel('Benson14 Eccentricity [deg]')
    ax.set_ylabel('Adjusted Eccentricity [deg]')
    ax.set_xlim([0.1,100])
    ax.set_ylim([0.1,100])
    ax.set_xticks([0.1, 1, 10, 100])
    ax.set_yticks([0.1, 1, 10, 100])
    ax.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    ax.get_yaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    if fname is not None:
        plt.savefig(fname, dpi=256, bbox_inches = "tight")
    return fig

def plot_distributionECCvalues(r0,r,scale,shape=0.75,min_eccen=0,max_eccen=90,fname=None):
    #fname = str(Save_DIR / subs) + '_' + hemlabel + '_DistributionECCValues.png'

    # What is the distribution of values like before and after adjustment?
    # The red dotted line is the H&H prediction (specifically for histogram
    # bins of the size below).
    (fig,axs) = plt.subplots(1,2, figsize=(5,2), dpi=256)
    fig.subplots_adjust(0,0,1,1,0.35,0)
    for (ax,prop,tag) in zip(axs, [r0, r], ['Benson14', 'Adjusted']):
        nbins = 50
        ax.hist(prop, bins=nbins, range=(min_eccen, max_eccen), weights=sarea/100)
        ax.set_xlabel(f'{tag} Eccentricity [deg]')
        ax.set_ylabel(f'Surface Area [cm$^2$]')
        bin_hwidth = ((max_eccen - min_eccen) / (nbins + 1)) / 2
        x = np.linspace(min_eccen+bin_hwidth, max_eccen - bin_hwidth, 250)
        b0 = x - bin_hwidth
        b1 = x + bin_hwidth
        b0[b0 < min_eccen] = min_eccen
        b1[b1 > max_eccen] = max_eccen
        y = (scale / (shape + x))**2 * np.pi*(b1**2 - b0**2)/2 / 100
        ax.plot(x, y, 'r:')
        ax.set_xscale('log')
        ax.set_xlim(left=0.5)
        ax.set_xticks([1,10,100])
        ax.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())

        if fname is not None:
            plt.savefig(fname, dpi=256, bbox_inches = "tight")
        return fig

def plot_comparisonCorticalECC(hem,v1,r0,r,existpRF=1,fname=None):
    # fname = str(Save_DIR / subs) + '_' + hemlabel + '_ComparisonCorticalECCmaps.png'

    # Let's plot comparison maps:
    if existpRF == 1:
        (fig,axs) = plt.subplots(1, 3, figsize=(6, 2), dpi=256)
    else:
        (fig,axs) = plt.subplots(1, 2, figsize=(4, 2), dpi=256)

    fig.subplots_adjust(0,0,1,1,0,0)
    # Make a flatmap:
    flatmap = hem.mask_flatmap(
        ('parcellation', 43),
        map_right='right',
        radius=np.pi/2)
    ecc = np.zeros(flatmap.vertex_count)
    if existpRF == 1:
        cod = flatmap.prop('prf_variance_explained')
        for (ax, prop) in zip(axs, [r0, r, v1.prop('prf_eccentricity')]):
            ecc[flatmap.tess.index(v1.labels)] = prop
            ny.cortex_plot(
                flatmap,
                color=ecc,
                mask=((ecc > 0) & (ecc < 25) & (cod > 0.04)),
                cmap='eccentricity', vmin=0, vmax=90,
                axes=ax)
            ax.axis('off')
        axs[0].set_title('Benson14')
        axs[1].set_title('Adjusted')
        axs[2].set_title('pRFs')
    else:
        for (ax, prop) in zip(axs, [r0, r]):
            ecc[flatmap.tess.index(v1.labels)] = prop
            ny.cortex_plot(
                flatmap,
                color=ecc,
                mask=((ecc > 0) & (ecc < 25)),
                cmap='eccentricity', vmin=0, vmax=90,
                axes=ax)
            ax.axis('off')
        axs[0].set_title('Benson14')
        axs[1].set_title('Adjusted')
        
        if fname is not None:
            plt.savefig(fname, dpi=256, bbox_inches = "tight")
        return fig

def plot_comparisonECCvsNative(v1,r0,r,existpRF=1,fname=None):
    if existpRF == 1:
        # fname = str(Save_DIR / subs) + '_' + hemlabel + '_ComparisonECCvsNative.png'
        
        density_plot = True
        (fig,axs) = plt.subplots(1,2, figsize=(5,2), dpi=256)
        fig.subplots_adjust(0,0,1,1,0.35,0)
        ecc = v1.prop('prf_eccentricity')
        cod = v1.prop('prf_variance_explained')
        #ii = (ecc < 20) & (cod > 0.04)
        ii = cod > 0.04
        (pmin,pmax) = (0.1, 100)
        (logpmin, logpmax) = (np.log(pmin), np.log(pmax))
        logpdif = logpmax - logpmin
        yy, xx = np.mgrid[logpmin:logpmax:250j, logpmin:logpmax:250j]
        positions = np.vstack([xx.ravel(), yy.ravel()])
        for (ax,prop,tag) in zip(axs, [r0, r], ['Benson14', 'Adjusted']):
            if density_plot:
                values = np.vstack([prop[ii], ecc[ii]])
                kernel = sp.stats.gaussian_kde(np.log(values), 0.1)
                f = np.reshape(kernel(positions).T, xx.shape)
                ax.imshow(f, cmap='Blues', vmin=0, vmax=0.8)
                ax.plot([0,249], [0,249], 'r:', lw=0.5)
                # Draw a line at 20°
                ln = 249*(np.log([20,20]) - logpmin)/logpdif
                ax.plot([0,250], ln, '--', c='0.5', lw=0.25)
                ax.plot(ln, [0,250], '--', c='0.5', lw=0.25)
                ax.invert_yaxis()
                ax.set_xticks(250*(np.log([0.1, 1, 10, 100]) - logpmin)/logpdif)
                ax.set_yticks(250*(np.log([0.1, 1, 10, 100]) - logpmin)/logpdif)
                ax.set_xticklabels([0.1,1,10,100])
                ax.set_yticklabels([0.1,1,10,100])
                ax.set_xlim([0,250])
                ax.set_ylim([0,250])
            else:
                ax.loglog(prop[ii], ecc[ii], 'k.', alpha=0.15)
                ax.loglog([0.1,100], [0.1,100], 'r:', zorder=-1)
                ax.set_xlim([0.1,100])
                ax.set_ylim([0.1,100])
                ax.set_xticks([0.1, 1, 10, 100])
                ax.set_yticks([0.1, 1, 10, 100])
                ax.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
                ax.get_yaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
            ax.set_xlabel(f'{tag} Eccentricity [deg]')
            ax.set_ylabel(f'pRF Eccentricity [deg]')
        
        if fname is not None:
            plt.savefig(fname, dpi=256, bbox_inches = "tight")
        return fig
    else:
        print("Data from pRF mapping experiment required!")

def plot_CMF(v1,r0,r,scale,shape=0.75,existpRF=1,fname=None):
    # Check the Cortical Magnification; a good way to calculate this is
    # to sort the vertices by eccentricity then look at a sliding window of
    # the vertices at a time.

    # fname = str(Save_DIR / subs) + '_' + hemlabel + '_CMF.png'

    (fig,ax) = plt.subplots(1,1, figsize=(5,3), dpi=256)
    if existpRF == 1:
        r_prf = v1.prop('prf_eccentricity')
        sarea = v1.prop('midgray_surface_area')
        eccs = [r0, r, r_prf]
        tags = ['Benson14','Adjusted','pRF']
        clrs = ['r','c','0.5']
    else:
        sarea = v1.prop('midgray_surface_area')
        eccs = [r0, r]
        tags = ['Benson14','Adjusted']
        clrs = ['r','c']

    for (prop,tag,c) in zip(eccs, tags, clrs):
        (x, cm) = cmag(prop, sarea)
        ax.loglog(x, cm, '-', c=c, lw=0.5, label=tag)
        ax.set_xlabel('Eccentricity [deg]')
        ax.set_ylabel(r'Cortical Magnification [mm$^2$/deg$^2$]')
    # We also plot the Horton and Hoyt cmag as a reference.
    x = np.linspace(0.25, 20, 500)
    hh_cmag = ((scale / (x + shape)))**2
    ax.loglog(x, hh_cmag, 'k:', label='H&H')
    ax.set_xlim([0.5, 21])
    ax.set_ylim([0.25, 200])
    ax.set_xticks([0.5,1,2,5,10,20])
    ax.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    ax.get_yaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    plt.legend()

    if fname is not None:
        plt.savefig(fname, dpi=256, bbox_inches = "tight")
    return fig
