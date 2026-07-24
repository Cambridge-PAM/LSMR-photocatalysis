from functions import *
import matplotlib.pyplot as plt
from matplotlib import cm as cm
from matplotlib import colors as colors

# functions for batch reaction plotting

ev = (r"[EV$^{+\!\cdot\!}$]", r"[EV$^{2\!+\!}$]$_0$", 600, 1.22e6, 1e-2)
bv = (r"[BV$^{+\!\cdot\!}$]", r"[BV$^{2\!+\!}$]$_0$", 535, 1.4e6, 1e-2)
rubpy = (r"[Ru(bpy)$_3^{2\!+\!}$]", r"[Ru(bpy)$_3^{2\!+\!}$]$_0$", 453, 1.482e6, 1e-2)

cstock = 1e-3 # [] in stock solution
vsol = 3e-3 # total volume of reaction mixture

# simple selector to determine which viologen labels and values to use 
def species(x):
    if x == 'E':
        return ev
    elif x == 'B':
        return bv
    elif x == 'Rubpy':
        return rubpy

# concentration vs time plot with linear regression
def concvstimeplot(data, v, start, x, duration=1e10, reg=True,
                   istart=1, iend=1e10, i=None, r2threshold=0, interceptthreshold=1e10):

    X, X0, peak, e, l = species(x)
    c0 = solprop(cstock=cstock, valiquot=v*1e-6, vsol=vsol) # initial concentration
    t, c = concvstime(data, peak, e, l, start=start, duration=duration)
    cscale = [conc*1e6 for conc in c] # plot concentrations in μM

    fig, ax = plt.subplots()
    ax.plot(t, cscale, color="black", alpha=0.2)
    if reg:
        fit, i = linreg(t, c, istart=istart, iend=iend, i=i, r2threshold=r2threshold, interceptthreshold=interceptthreshold)
        x1 = np.linspace(0, t[i], 100) # plot linear regression line over the initial linear region
        y1 = (fit.slope*x1 + fit.intercept)*1e6
        ax.plot(x1, y1, color="black")
        x2 = np.linspace(t[i], t[i]*2, 100) # extrapolate linear regression line
        y2 = (fit.slope*x2 + fit.intercept)*1e6
        ax.plot(x2, y2, color="black", linestyle="--")
        ax.text(0.95, 0.05,
                rf'$k$ = {errformat(fit.slope/c0, fit.stderr/c0, prefix=1e-6)} $\times$ 10$^{{-6}}$ s$^{{-1}}$''\n'
                rf'rate = slope = {errformat(fit.slope, fit.stderr, prefix=1e-9)} nM s$^{{-1}}$''\n'
                f'intercept = {errformat(fit.intercept, fit.intercept_stderr, prefix=1e-6)} µM\n'
                f'$r^2$ = {fit.rvalue**2:.3f}',
                transform=plt.gca().transAxes, ha="right", va="bottom")

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
    ax.set_ylabel(rf'{X} / μM')
    ax.set_title(rf'{X} / μM over time / s for {X0} = {c0*1e6:.0f} μM')
    plt.show()

# lnc vs time plot for a first order reaction with linear regression
def lncvstimeplot(data, v, start, x, duration=1e10, reg=True,
                  istart=1, iend=1e10, i=None, r2threshold=0, interceptthreshold=1e10):

    X, X0, peak, e, l = species(x)
    c0 = solprop(cstock=cstock, valiquot=v*1e-6, vsol=vsol) # initial concentration
    t, c = concvstime(data, peak, e, l, start=start, duration=duration)
    lnc = [np.log(c0-conc) for conc in c]

    fig, ax = plt.subplots()
    ax.plot(t, lnc, color="black", alpha=0.2) # plot time in s, concentrations in M
    if reg:
        fit, i = linreg(t, lnc, istart=istart, iend=iend, i=i, r2threshold=r2threshold, interceptthreshold=interceptthreshold)
        x1 = np.linspace(0, t[i], 100) # plot linear regression line over the initial linear region
        y1 = (fit.slope*x1 + fit.intercept)
        ax.plot(x1, y1, color="black")
        x2 = np.linspace(t[i], t[i]*2, 100) # extrapolate linear regression line
        y2 = (fit.slope*x2 + fit.intercept)
        ax.plot(x2, y2, color="black", linestyle="--")
        ax.text(0.95, 0.95,
            rf'$k$ = $-$slope = {errformat(-fit.slope, fit.stderr, prefix=1e-6)} $\times$ 10$^{{-6}}$ s$^{{-1}}$''\n'
            f'intercept = {errformat(fit.intercept, fit.intercept_stderr, sci=False)}\n'
            f'$r^2$ = {fit.rvalue**2:.3f}',
            transform=plt.gca().transAxes, ha="right", va="top")

    fig.canvas.draw() # calculate ticks preliminarily
    xscale = np.diff(ax.get_xticks())[0]
    yscale = np.diff(ax.get_yticks())[0]
    xmin = np.floor(min(t)/xscale)*xscale
    xmax = np.ceil(max(t)/xscale)*xscale
    ymin = np.floor((min(lnc)-0.2*yscale)/yscale)*yscale
    ymax = np.ceil((max(lnc)+0.2*yscale)/yscale)*yscale

    ax.set_xlim(xmin, xmax)
    ax.set_xticks(np.arange(xmin, xmax+abs(0.001*xmax), xscale))
    ax.set_ylim(ymin, ymax)
    ax.set_yticks(np.arange(ymin, ymax+abs(0.001*ymax), yscale))
    ax.set_xlabel(r'time / s')
    ax.set_ylabel(rf'ln({X0}$-${X} / M)')
    ax.set_title(rf'ln({X0}$-${X} / M) over time / s for {X0} = {c0*1e6:.0f} μM')
    plt.show()

# uv-vis absorption spectrum over time for each run
def spectrumvstimeplot(data, v, start, x, duration=1e10, lmin=450, lmax=700, vmin=None, vmax=None):

    X, X0, _, _, _ = species(x)
    c0 = solprop(cstock=cstock, valiquot=v*1e-6, vsol=vsol) # initial concentration
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

    norm = colors.Normalize(vmin=vmin, vmax=vmax)
    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_ticks(np.arange(vmin, vmax*1.001, 0.1))
    cbar.set_label("absorbance")

    fig.canvas.draw() # calculate ticks preliminarily
    xscale = np.diff(ax.get_xticks())[0]
    xmin = np.floor(min(t)/xscale)*xscale
    xmax = np.ceil(max(t)/xscale)*xscale

    ax.set_xlim(xmin, xmax)
    ax.set_xticks(np.arange(xmin, xmax+abs(0.001*xmax), xscale))
    ax.set_ylim(lmin, lmax)
    ax.set_xlabel("time / s")
    ax.set_ylabel(r'$\lambda$ / nm')
    ax.set_title(rf"UV-Vis spectra of {X} over time / s for {X0} = {c0*1e6:.0f} μM")
    plt.show()

# uv-vis absorption spectrum at time t for each run
def spectrumvsconcplot(root, runs, x, t, lmin=450, lmax=700):

    X, X0, _, _, _ = species(x)
    fig, ax = plt.subplots(figsize=(8,5))
    cmap = plt.get_cmap("viridis_r")
    col = cmap(np.linspace(0, 1, len(runs)))

    amax = 0
    for i, (v, start) in enumerate(runs): # exclude the run on 260625 which did not include Rubpy in background measurement
        data = f'{root}{v}_EDTA_Buffered'
        l, a = spectrumvstimepoints(data, [t], start=start)
        mask = (lmin <= l) & (l <= lmax)
        l = l[mask]
        a = a[:, mask]
        ax.plot(l, a[0], color=col[i])
        
        if np.max(a) > amax:
            amax = np.max(a)

    norm = colors.Normalize(vmin=0, vmax=len(runs)-1)
    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_ticks(range(len(runs)))
    cbar.set_ticklabels([f'{solprop(cstock=cstock, valiquot=v*1e-6, vsol=vsol)*1e6:.0f}' for v, _ in runs])
    cbar.set_label(f"{X0} / μM")

    fig.canvas.draw() # calculate ticks preliminarily
    yscale = np.diff(ax.get_yticks())[0]
    ymax = np.ceil((amax+0.2*yscale)/yscale)*yscale

    ax.set_xlim(lmin, lmax)
    ax.set_ylim(0, ymax)
    ax.set_yticks(np.arange(0, ymax+abs(0.001*ymax), yscale))
    ax.set_xlabel(r'$\lambda$ / nm')
    ax.set_ylabel(r'absorbance')
    ax.set_title(f'UV-Vis spectra of {X} after {t/60} min of irradiation')
    plt.show()