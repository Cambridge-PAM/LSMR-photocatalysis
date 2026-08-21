import numpy as np
import scipy.stats as stats
import scipy.constants as constants
from pathlib import Path

# provide three of the four parameters (a, e(L mol-1 m-1), c(mol/L), l(m)) to calculate the missing one for the Beer-Lambert law.
def beerlambert(a=None, e=None, c=None, l=None):
    try:
        if a is None:
            return e * c * l
        elif e is None:
            return a / (c * l)
        elif c is None:
            return a / (e * l)
        elif l is None:
            return a / (e * c)
        else:
            raise ValueError
    except:
        print("Error: three of the four parameters (a, e(L mol-1 m-1), c(mol/L), l(m)) needed to calculate the missing one.")
        return None

# provide two of the three parameters (m(g), M(g/mol), mol) to calculate the missing one for molar calculations.
def molcalc(m=None, M=None, mol=None):
    try:
        if m is None:
            return mol * M
        elif M is None:
            return m / mol
        elif mol is None:
            return m / M
        else:
            raise ValueError
    except:
        print("Error: two of the three parameters (m(g), M(g/mol), mol) needed to calculate the missing one.")
        return None

# provide the three parameters (mol, v(L), c(mol/L)) to calculate the missing one for solution calculations.
def solcalc(mol=None, v=None, c=None):
    try:
        if mol is None:
            return v * c
        elif v is None:
            return mol / c
        elif c is None:
            return mol / v
        else:
            raise ValueError
    except:
        print("Error: two of the three parameters (mol, v(L), c(mol/L)) needed to calculate the missing one.")
        return None

# provide the four parameters (M(g/mol), mstock(g), vstock(L), cstock(mol/L), valiquot(L), vsol(L), csol(mol/L)) to calculate the missing one for solution preparation calculations.
def solprop(M=None, mstock=None, vstock=None, cstock=None, valiquot=None, vsol=None, csol=None):
    try:
        if cstock is None:    
            cstock = solcalc(mol=molcalc(m=mstock, M=M), v=vstock) # concentration of stock solution
    except:
        print("Error: (M(g/mol), mstock(g), vstock(L)) or cstock(mol/L) needed")
        return None
    try:
        if csol is None:
            return solcalc(mol=solcalc(c=cstock, v=valiquot), v=vsol) # concentration of overall solution
        elif vsol is None:
            return solcalc(mol=solcalc(c=cstock, v=valiquot), c=csol) # volume in overall solution
        elif valiquot is None:
            return solcalc(mol=solcalc(c=csol, v=vsol), c=cstock) # volume of aliquot taken from stock solution
    except:
        print("Error: two of the three parameters (valiquot(L), vsol(L), csol(mol/L)) needed to calculate the missing one.")
        return None

# provide two of the three parameters (flowrate, tR, vol) in any consistent units to calculate the missing one for flow calculations.
def flowprop(flowrate=None, tR=None, vol=None):
    try:
        if flowrate is None:
            return vol / tR
        elif tR is None:
            return vol / flowrate
        elif vol is None:
            return flowrate * tR
    except:
        print("Error: two of the three parameters (flowrate, tR, vol) needed to calculate the missing one.")
        return None

# return M s-1 of photons for a given power(W) and wavelength(m) and volume(L)
def ratephotons(power, wavelength, volume):
    h = constants.physical_constants['Planck constant'][0]
    c = constants.physical_constants['speed of light in vacuum'][0]
    n = constants.physical_constants['Avogadro constant'][0]
    return power/(volume*n*h*c/wavelength)

# return the elapsed time between two timestamps in the format hh-mm-ss-ms.
def elapsed(start, end):

    def parse(time):
        parts = time.split("-")
        if len(parts) != 4:
            raise ValueError("Timestamps must use the format hh-mm-ss-ms.")

        h, m, s, ms = (int(part) for part in parts)
        if not (0 <= m < 60 and 0 <= s < 60 and 0 <= ms < 1000):
            raise ValueError("Minutes, seconds, and milliseconds must be valid time units.")

        return h * 3600 + m * 60 + s + ms * 0.001

    return parse(end) - parse(start)

# read absorbances at peak wavelength and within peak width
def readpeak(file, peak, pwidth):
    absorbances = [] # list of absorbances to average
    for line in file:
        try:
            wavelength = float(line.split("	")[0])
            if peak-pwidth/2 <= wavelength <= peak+pwidth/2: # search for absorbance peak
                absorbances.append(float(line.split("	")[1]))
            elif peak+pwidth/2 <= wavelength: # finish search for absorbance peak
                break
        except:
            continue
    return absorbances

# concentration of compound over time from oceanoptics data, optional linear regression for initial first order rate
def concvstime(data, peak, e, l, start="", end="", duration=1e10, pwidth=10):

    a = [] # array of absorbances
    t = [] # array of time in seconds
    read = False

    for f in sorted(Path(data).glob('*.txt'), key=lambda f: f.stem[-12:]): # sort files by timestamp
        timestamp = f.stem[-12:]
        if start == "":
            start = timestamp # if no start time defined, start with the first file
        
        time = elapsed(start, timestamp) # convert timestamp to time in seconds

        if time >= duration or timestamp == end:
            read = False # stop reading files
            break
        elif time >= 0:
            read = True # start reading files
        
        if read:
            t.append(time) # list of time points
            with f.open("r") as file:
                absorbances = readpeak(file, peak, pwidth)
                a.append(sum(absorbances)/len(absorbances)) # average the absorbances

    c = [beerlambert(a=absorbance, e=e, l=l) for absorbance in a] # convert absorbance to concentration using Beer-Lambert law
    return t, c

# concentration of compound at given points from oceanoptics data
def concvstimepoints(data, peak, e, l, points=None, start="", pwidth=10, twidth=5):
    
    a = [] # array of absorbances
    files = sorted(Path(data).glob('*.txt'), key=lambda f: f.stem[-12:]) # sort files by timestamp
    
    for i in range(len(files)):
        
        if points: # read files from certain time points (in seconds)
            timestamp = files[i].stem[-12:]
            
            if start == "": 
                start = timestamp # if no start time defined, start with the first file
            
            time = round(elapsed(start, timestamp)) # convert timestamp to time in seconds
            
            if time in points:
                absorbances = []
                for j in range(-twidth//2 + 1, twidth//2 + 1): # search for absorbance peak in files within twidth seconds of timestamp
                    with files[i+j].open("r") as file:
                        absorbances += readpeak(file, peak, pwidth)
                a.append(sum(absorbances)/len(absorbances)) # average the absorbances

        else: # read all files
            with files[i].open("r") as file:
                absorbances = readpeak(file, peak, pwidth)
                a.append(sum(absorbances)/len(absorbances)) # average the absorbances

    c = [beerlambert(a=absorbance, e=e, l=l) for absorbance in a] # convert absorbance to concentration using Beer-Lambert law
    return c

# reads the uv-vis spectrum and appends data to wavelengths l and absorbances a
def readspectrum(file, l, a):
    a.append([]) # list of absorbances for this time point
    for line in file:
        try: # line contains spectrum data
            if len(a) == 1: # if first instance
                l.append(float(line.split("	")[0])) # create list of wavelengths
            a[-1].append(float(line.split("	")[1])) # append absorbances for this time point
        except: # line does not contain spectrum data
            continue

# uv-vis absorbance spectrum over time from oceanoptics data
def spectrumvstime(data, start="", end="", duration=1e10):

    t = [] # time in seconds
    l = [] # wavelengths
    a = [] # absorbances
    read = False

    for f in sorted(Path(data).glob('*.txt'), key=lambda f: f.stem[-12:]): # sort files by timestamp
        timestamp = f.stem[-12:]
        if start == "":
            start = timestamp # if no start time defined, start with the first file
        
        time = elapsed(start, timestamp) # convert timestamp to time in seconds

        if time >= duration or timestamp == end:
            read = False # stop reading files
            break
        elif time >= 0:
            read = True # start reading files
        
        if read:
            t.append(time) # list of time points 
            with f.open("r") as file:
                readspectrum(file, l, a)
    
    return np.asarray(t), np.asarray(l), np.array(a).T # convert lists to numpy array and transpose a to have wavelengths as rows and time points as columns

# uv-vis absorbance spectrum at given time points (in seconds) from oceanoptics data
def spectrumvstimepoints(data, points=None, start=""):

    l = [] # wavelengths
    a = [] # absorbances
    
    for f in sorted(Path(data).glob('*.txt'), key=lambda f: f.stem[-12:]): # sort files by timestamp
        
        if points: # read files from certain time points (in seconds)
            timestamp = f.stem[-12:]
            if start == "":
                start = timestamp # if no start time defined, start with the first file
 
            time = round(elapsed(start, timestamp)) # convert timestamp to time in seconds

            if time in points:
                with f.open("r") as file:
                    readspectrum(file, l, a)

        else: # read all files
            with f.open("r") as file:
                readspectrum(file, l, a)

    return np.asarray(l), np.asarray(a)

# extract uv-vis absorption spectrum from oceanoptics data
def spectrum(data):
    
    l = [] # wavelengths
    a = [] # absorbances

    with open(data, "r") as file:
        for line in file:
            try: # line contains spectrum data
                l.append(float(line.split("	")[0]))
                a.append(float(line.split("	")[1]))
            except: # line does not contain spectrum data
                continue
    
    return l, a

# running linear regression analysis, starts selecting points from istart, adds more points until iend or r2 and intercept exceed threshold
def linreg(t, c, istart=1, iend=1e10, i=None, r2threshold=0, interceptthreshold=1e10):
    if iend == 1e10:
        iend = len(t)-1
    if i:
        istart = i
        iend = i
    for i in range(istart, iend+1):
        fit = stats.linregress(t[:i+1], c[:i+1])
        if fit.rvalue**2 <= r2threshold or abs(fit.intercept) >= abs(interceptthreshold): # stop when r2 or intercept exceeds thresholds
            fit = stats.linregress(t[:i], c[:i]) # return the last fit
            return fit, i
    return fit, i

# format numbers and error to appropriate decimal places
def errformat(val, err, prefix=1, sci=True):

    val = val*(1/prefix) if not np.isnan(val) else 0 # scale by the unit prefix
    err = err*(1/prefix) if not np.isnan(err) else 0

    i = int(f"{err:.0e}".split("e")[-1])+1 # find the exponent such that the error has 1sf in the first decimal place

    if sci and prefix==1: # scientific notation to 1dp by default
        return rf"${val*(10**-i):.1f}$ ± ${err*(10**-i):.1f}$ $\times$ 10$^{{{i}}}$"
    else: # use more dp (implied default if a prefix is specified)
        return rf"${val:.{max(-i+1,0)}f}$ ± ${err:.{max(-i+1,0)}f}$"