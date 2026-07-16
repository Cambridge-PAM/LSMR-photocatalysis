from functions import *
import matplotlib.pyplot as plt
from matplotlib import cm as cm
from matplotlib import colors as colors

# stop flow with reaction in 50 ul chip and PtOEP/DPA film
# irradiation at each time point, then uv-vis immediately after in the same chip
# flushed with fresh starting reaction mixture and repeated for new time point

ev = (r"[EV$^{+\!\cdot\!}$]", r"[EV$^{2\!+\!}$]$_0$", 600, 1.22e6, 700e-6)
bv = (r"[BV$^{+\!\cdot\!}$]", r"[BV$^{2\!+\!}$]$_0$", 535, 1.4e6, 700e-6)
area = np.pi*(9.5e-3/2)**2 # m-2, area of thorlabs power meter detector

# simple selector to determine which viologen labels and values to use 
def viologen(x):
    if x == 'E':
        return ev
    elif x == 'B':
        return bv

# concentration vs time plot with linear regression
def concvstimeplot(data, t, c0, x, film, light, wavelength, power,
                   istart=1, iend=1e10, i=None, r2threshold=0, interceptthreshold=1e10):

    XV, XV0, peak, e, l = viologen(x)
    c = concvstimepoints(data, peak, e, l)
    cscale = [conc*1e3 for conc in c] # plot concentrations in mM
    n = photons(power, wavelength, l*area*1e3) # M s-1 of photons

    fig, ax = plt.subplots()
    ax.plot(t, cscale, 'x', color='black')
    fit, i = linreg(t, c, istart=istart, iend=iend, i=i, r2threshold=r2threshold, interceptthreshold=interceptthreshold)
    x = np.linspace(t[0], t[i], 100) # plot linear regression line over the initial linear region
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

    ax.set_xlim(np.floor(min(t)/xscale)*xscale, np.ceil(max(t)/xscale)*xscale)
    ax.set_ylim(min(np.floor((min(cscale)-0.2*yscale)/yscale)*yscale, 0), np.ceil((max(cscale)+0.2*yscale)/yscale)*yscale)
    ax.set_xlabel(r'time / s')
    ax.set_ylabel(rf'{XV} / mM')
    ax.set_title(rf'{XV} / mM over time / s for {XV0} = {c0*1e3:.0f} mM''\n'
                 rf'with {film} film under {light} ({int(wavelength*1e9)} nm) light at {power*1e3:.1f} mW')
    plt.show()

# lnc vs time plot for a first order reaction with linear regression
def lncvstimeplot(data, t, c0, x, film, light, wavelength, power,
                  istart=1, iend=1e10, i=None, r2threshold=0, interceptthreshold=1e10):

    XV, XV0, peak, e, l = viologen(x)
    c = concvstimepoints(data, peak, e, l)
    lnc = [np.log(c0-conc) for conc in c]
    n = photons(power, wavelength, l*area*1e3) # M s-1 of photons
    
    fig, ax = plt.subplots()
    ax.plot(t, lnc, 'x', color='black')
    fit, i = linreg(t, lnc, istart=istart, iend=iend, i=i, r2threshold=r2threshold, interceptthreshold=interceptthreshold)
    x = np.linspace(t[0], t[i], 100) # plot linear regression line over the initial linear region
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

    ax.set_xlim(np.floor(min(t)/xscale)*xscale, np.ceil(max(t)/xscale)*xscale)
    ax.set_ylim(np.floor((min(lnc)-0.2*yscale)/yscale)*yscale, np.ceil((max(lnc)+0.2*yscale)/yscale)*yscale)
    ax.set_xlabel(r'time / s')
    ax.set_ylabel(rf'ln({XV0}$-${XV} / M)')
    ax.set_title(rf'ln({XV0}$-${XV} / M) over time / s for {XV0} = {c0*1e3:.0f} mM''\n'
                 rf'with {film} film under {light} ({int(wavelength*1e9)} nm) light at {power*1e3:.1f} mW')
    plt.show()

# absorbance spectra for each t
def spectrumvstimeplot(data, t, c0, x, film, light, wavelength, power):

    XV, XV0, _, _, _ = viologen(x)

    fig, ax = plt.subplots(figsize=(8, 5))
    cmap = plt.get_cmap("viridis_r")
    col = cmap(np.linspace(0, 1, len(t)))

    l, a = spectrumvstimepoints(data)
    mask = (l >= 450) & (l <= 700)
    l = l[mask]
    a = a[:, mask]

    for i in range(len(t)):
        ax.plot(l, a[i], color=col[i])

    norm = colors.Normalize(0, len(t)-1)
    sm = cm.ScalarMappable(norm=norm, cmap=cmap)

    cbar = fig.colorbar(sm, ax=ax)
    cbar.set_ticks(range(len(t)))
    cbar.set_ticklabels(np.asarray(t)/60)
    cbar.set_label(f"time / min")

    fig.canvas.draw() # calculate ticks preliminarily
    yscale = np.diff(ax.get_yticks())[0]

    ax.set(xlim=(450, 700), ylim=(0, np.ceil((np.max(a)+0.2*yscale)/yscale)*yscale),
           xlabel=r'$\lambda$ / nm', ylabel='absorbance')
    ax.set_title(rf'UV-Vis spectra over time for {XV0} = {c0*1e3:.0f} mM''\n'
                 rf'with {film} film under {light} ({int(wavelength*1e9)} nm) light at {power*1e3:.1f} mW')
    plt.show()