from functions import *
import matplotlib.pyplot as plt
from matplotlib import cm as cm
from matplotlib import colors as colors

# stop flow with reaction in 50 ul chip (and film if any)
# irradiation at each time point, then uv-vis immediately after in the same chip
# flushed with fresh starting reaction mixture and repeated for new time point

tR = "$t_\\mathrm{R}$" # tR
evold = (r"[EV$^{+\!\cdot\!}$]", r"[EV$^{2\!+\!}$]$_0$", 600, 1.22e6, 700e-6)
ev = (r"[EV$^{+\!\cdot\!}$]", r"[EV$^{2\!+\!}$]$_0$", 600, 1.22e6, 800e-6)
bv = (r"[BV$^{+\!\cdot\!}$]", r"[BV$^{2\!+\!}$]$_0$", 535, 1.4e6, 800e-6)
area = np.pi*(9.5e-3/2)**2 # m-2, area of thorlabs power meter detector

# simple selector to determine which viologen labels and values to use 
def viologen(x):
    if x == 'E':
        return ev
    elif x == 'B':
        return bv
    elif x == 'Eold':
        return evold

# time point intervals in s between flow rate switches
def flowpoints(tRs, tbuffer, vrxn, vfull):
    flowrates = [flowprop(tR=tR, vol=vrxn) for tR in tRs] # flow rates for each tR
    points = []
    for i in range(len(tRs)):
        points.append(round(flowprop(flowrate=flowrates[i], vol=vfull) + 2*tbuffer))
    return points

# time points in s to record absorbances, cumulative from the start of experiment
def recpoints(tRs, tbuffer, vrxn, vfull):
    flowrates = [flowprop(tR=tR, vol=vrxn) for tR in tRs] # flow rates for each tR
    points = []
    for i in range(len(tRs)):
        if i == 0:
            prevtime = 0
        else:
            prevtime = points[-1] + tbuffer
        points.append(round(prevtime + flowprop(flowrate=flowrates[i], vol=vfull) + tbuffer))
    return points

# concentration vs time point plot with linear regression
def concvstimeplot(data, tRs, points, start, c0, x, film, light, wavelength, power,
                   istart=1, iend=1e10, i=None, r2threshold=0, interceptthreshold=1e10):

    XV, XV0, peak, e, l = viologen(x)
    c = concvstimepoints(data, peak, e, l, points=points, start=start)
    cscale = [conc*1e3 for conc in c] # plot concentrations in mM
    n = ratephotons(power, wavelength, l*area*1e3) # M s-1 of photons

    fig, ax = plt.subplots()
    ax.plot(tRs, cscale, 'x', color='black')
    fit, i = linreg(tRs, c, istart=istart, iend=iend, i=i, r2threshold=r2threshold, interceptthreshold=interceptthreshold)
    x = np.linspace(tRs[0], tRs[i], 100) # plot linear regression line over the initial linear region
    y = (fit.slope*x + fit.intercept)*1e3
    ax.plot(x, y, color='black')
    ax.text(0.95, 0.05,
            rf'$\Phi$ = $\mathrm{{\frac{{rate_{{reaction}}}}{{rate_{{photon}}}}}}$ = {errformat(fit.slope/n, fit.stderr/n, 1e-2)} %''\n'
            rf'$k$ = {errformat(fit.slope/c0, fit.stderr/c0, prefix=1e-6)} $\times$ 10$^{{-6}}$ s$^{{-1}}$''\n'
            rf'rate = slope = {errformat(fit.slope, fit.stderr, prefix=1e-6)} $\mathrm{{\mu}}$M s$^{{-1}}$''\n'
            f'intercept = {errformat(fit.intercept, fit.intercept_stderr, prefix=1e-3)} mM\n'
            f'$r^2$ = {fit.rvalue**2:.3f}',
            transform=ax.transAxes, ha="right", va="bottom")

    fig.canvas.draw() # calculate ticks preliminarily
    xscale = np.diff(ax.get_xticks())[0]
    yscale = np.diff(ax.get_yticks())[0]
    xmin = np.floor(min(tRs)/xscale)*xscale
    xmax = np.ceil(max(tRs)/xscale)*xscale
    ymin = min(np.floor((min(cscale)-0.2*yscale)/yscale)*yscale, 0)
    ymax = np.ceil((max(cscale)+0.2*yscale)/yscale)*yscale

    ax.set_xlim(xmin, xmax)
    ax.set_xticks(np.arange(xmin, xmax+abs(0.0001*xmax), xscale))
    ax.set_ylim(ymin, ymax)
    ax.set_yticks(np.arange(ymin, ymax+abs(0.0001*ymax), yscale))
    ax.set_xlabel(rf'{tR} / s')
    ax.set_ylabel(rf'{XV} / mM')
    ax.set_title(rf'{XV} / mM against {tR} / s for {XV0} = {c0*1e3:.0f} mM''\n'
                 rf'with {film} film under {light} ({int(wavelength*1e9)} nm) light at {power*1e3:.1f} mW')
    plt.show()

# lnc vs time point plot for a first order reaction with linear regression
def lncvstimeplot(data, tRs, points, start, c0, x, film, light, wavelength, power,
                  istart=1, iend=1e10, i=None, r2threshold=0, interceptthreshold=1e10):

    XV, XV0, peak, e, l = viologen(x)
    c = concvstimepoints(data, peak, e, l, points=points, start=start)
    lnc = [np.log(c0-conc) for conc in c]
    n = ratephotons(power, wavelength, l*area*1e3) # M s-1 of photons
    
    fig, ax = plt.subplots()
    ax.plot(tRs, lnc, 'x', color='black')
    fit, i = linreg(tRs, lnc, istart=istart, iend=iend, i=i, r2threshold=r2threshold, interceptthreshold=interceptthreshold)
    x = np.linspace(tRs[0], tRs[i], 100) # plot linear regression line over the initial linear region
    y = (fit.slope*x + fit.intercept)
    ax.plot(x, y, color='black')
    ax.text(0.95, 0.95,
            rf'$\Phi$ = $\mathrm{{\frac{{rate_{{reaction}}}}{{rate_{{photon}}}}}}$ = {errformat(-fit.slope*c0/n, fit.stderr*c0/n, 1e-2)} %''\n'
            rf'$k$ = $-$slope = {errformat(-fit.slope, fit.stderr, prefix=1e-6)} $\times$ 10$^{{-6}}$ s$^{{-1}}$''\n'
            f'intercept = {errformat(fit.intercept, fit.intercept_stderr, sci=False)}\n'
            f'$r^2$ = {fit.rvalue**2:.3f}',
            transform=ax.transAxes, ha="right", va="top")

    fig.canvas.draw() # calculate ticks preliminarily
    xscale = np.diff(ax.get_xticks())[0]
    yscale = np.diff(ax.get_yticks())[0]
    xmin = np.floor(min(tRs)/xscale)*xscale
    xmax = np.ceil(max(tRs)/xscale)*xscale
    ymin = np.floor((min(lnc)-0.2*yscale)/yscale)*yscale
    ymax = np.ceil((max(lnc)+0.2*yscale)/yscale)*yscale

    ax.set_xlim(xmin, xmax)
    ax.set_xticks(np.arange(xmin, xmax+abs(0.0001*xmax), xscale))
    ax.set_ylim(ymin, ymax)
    ax.set_yticks(np.arange(ymin, ymax+abs(0.0001*ymax), yscale))
    ax.set_xlabel(rf'{tR} / s')
    ax.set_ylabel(rf'ln({XV0}$-${XV} / M)')
    ax.set_title(rf'ln({XV0}$-${XV} / M) against {tR} / s for {XV0} = {c0*1e3:.0f} mM''\n'
                 rf'with {film} film under {light} ({int(wavelength*1e9)} nm) light at {power*1e3:.1f} mW')
    plt.show()

# absorbance spectra for each t
def spectrumvstimeplot(data, tRs, points, start, c0, x, film, light, wavelength, power, lmin=450, lmax=700):

    _, XV0, _, _, _ = viologen(x)
    fig, ax = plt.subplots(figsize=(8, 5))
    cmap = plt.get_cmap("viridis_r")
    col = cmap(np.linspace(0, 1, len(tRs)))

    l, a = spectrumvstimepoints(data, points=points, start=start)
    mask = (lmin <= l) & (l <= lmax)
    l = l[mask]
    a = a[:, mask]

    for i in range(len(tRs)):
        ax.plot(l, a[i], color=col[i])

    norm = colors.Normalize(0, len(tRs)-1)
    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    cbar = fig.colorbar(sm, ax=ax)
    cbar.set_ticks(range(len(tRs)))
    cbar.set_ticklabels(np.asarray(tRs)/60)
    cbar.set_label(rf"{tR} / min")

    fig.canvas.draw() # calculate ticks preliminarily
    yscale = np.diff(ax.get_yticks())[0]
    ymax = np.ceil((np.max(a)+0.2*yscale)/yscale)*yscale

    ax.set_xlim(lmin, lmax)
    ax.set_ylim(0, ymax)
    ax.set_yticks(np.arange(0, ymax+abs(0.0001*ymax), yscale))
    ax.set_xlabel(r'$\lambda$ / nm')
    ax.set_ylabel('absorbance')
    ax.set_title(rf'UV-Vis spectra against {tR} / s for {XV0} = {c0*1e3:.0f} mM''\n'
                 rf'with {film} film under {light} ({int(wavelength*1e9)} nm) light at {power*1e3:.1f} mW')
    plt.show()

# plots what the uv-vis detector sees over the course of the whole run
def detectorconcplot(data, points, start, c0, x, film, light, wavelength, power, duration=1e10):

    XV, XV0, peak, e, l = viologen(x)
    t, c = concvstime(data, peak, e, l, start=start, duration=duration)
    cscale = [conc*1e3 for conc in c] # plot concentrations in μM

    fig, ax = plt.subplots()
    ax.plot(t, cscale, color="black", alpha=0.2) # plot the overall detector 
    for point in points: 
        ax.axvline(x=point, color='black', linestyle='--') # plot vertical lines for each tR

    fig.canvas.draw() # calculate ticks preliminarily
    xscale = np.diff(ax.get_xticks())[0]
    yscale = np.diff(ax.get_yticks())[0]
    xmin = np.floor(min(t)/xscale)*xscale
    xmax = np.ceil(max(t)/xscale)*xscale
    ymin = np.floor((min(cscale)-0.2*yscale)/yscale)*yscale
    ymax = np.ceil((max(cscale)+0.2*yscale)/yscale)*yscale
    if abs(min(cscale)) < yscale: # if the smallest c value is approximately zero
        ymin = 0

    ax.set_xlim(xmin, xmax)
    ax.set_xticks(np.arange(xmin, xmax+abs(0.001*xmax), xscale))
    ax.set_ylim(ymin, ymax)
    ax.set_yticks(np.arange(ymin, ymax+abs(0.001*ymax), yscale))
    ax.set_xlabel(r'time / s')
    ax.set_ylabel(rf'{XV} / mM')
    ax.set_title(rf'Detected {XV} / mM over time / s for {XV0} = {c0*1e3:.0f} mM''\n'
                 rf'with {film} film under {light} ({int(wavelength*1e9)} nm) light at {power*1e3:.1f} mW')
    plt.show()

# uv-vis absorption spectrum over time for each run
def detectorspectrumplot(data, points, start, c0, x, film, light, wavelength, power, duration=1e10, 
                         lmin=450, lmax=700, vmin=None, vmax=None):

    XV, XV0, _, _, _ = viologen(x)
    t, l, a = spectrumvstime(data, start=start, duration=duration)
    mask = (lmin <= l) & (l <= lmax) 
    a = a[mask, :]
    l = l[mask]

    if vmin == None:
        vmin = np.floor(a.min()/0.1)*0.1 # round down to nearest 0.1
    if vmax == None:
        vmax = np.ceil(a.max()/0.1)*0.1 # round up to nearest 0.1
    
    fig, ax = plt.subplots(figsize=(8,3))
    cmap = plt.get_cmap("viridis_r")
    ax.contourf(t, l, a, levels=100, vmin=vmin, vmax=vmax, cmap=cmap)
    for point in points: 
        ax.axvline(x=point, color='black', linestyle='--') # plot vertical lines for each tR

    norm = colors.Normalize(vmin=vmin, vmax=vmax)
    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_ticks(np.arange(vmin, vmax*1.001, 0.1))
    cbar.set_label("absorbance")

    fig.canvas.draw() # calculate ticks preliminarily
    # xscale = np.diff(ax.get_xticks())[0]
    # xmin = np.floor(min(t)/xscale)*xscale
    # xmax = np.ceil(max(t)/xscale)*xscale

    # ax.set_xlim(xmin, xmax)
    # ax.set_xticks(np.arange(xmin, xmax+abs(0.0001*xmax), xscale))
    ax.set_ylim(lmin, lmax)
    ax.set_xlabel("time / s")
    ax.set_ylabel(r'$\lambda$ / nm')
    ax.set_title(rf"UV-Vis spectra of {XV} over time / s for {XV0} = {c0*1e3:.0f} mM""\n"
                 rf'with {film} film under {light} ({int(wavelength*1e9)} nm) light at {power*1e3:.1f} mW')
    plt.show()