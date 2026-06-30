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

# concentration of compound over time from oceanoptics data
def concvstime(data, peak, width, l, e, start="", end="", duration=None):

    a = [] # array of absorbances
    t = [] # array of time in seconds

    absorbances = [] # list of absorbances to average
    read = False

    for f in Path(data).glob('*.txt'):
        timestamp = f.name[-16:-4]
        if start == "": 
            start = timestamp # if no start time defined, start with the first file
        
        time = elapsed(start, timestamp) # convert timestamp to time in seconds

        if time >= duration or  timestamp == end:
            read = False # stop reading files
            break
        elif time >= 0:
            read = True # start reading files
        
        if read:
            t.append(time)
            with f.open("r") as file:
                absorbances = [] # list of absorbances to average
                for line in file:
                    try:
                        wavelength = float(line.split("	")[0])
                        if peak-width <= wavelength <= peak+width: # search for absorbance peak
                            absorbances.append(float(line.split("	")[1]))
                        elif peak+width <= wavelength: # finish search for absorbance peak
                            a.append(sum(absorbances)/len(absorbances)) # average the absorbances 
                            absorbances = [] # reset absorbances list
                            break
                    except:
                        pass
    
    c = [beerlambert(a=absorbance, l=l, e=e) for absorbance in a] # convert absorbance to concentration using Beer-Lambert law
    return t, c

# extract uv-vis absorption spectrum from oceanoptics data
def uvvis(data):
    
    l = []
    a = []

    with open(data, "r") as file:
        for line in file:
            try: # line contains spectrum data
                l.append(float(line.split("	")[0]))
                a.append(float(line.split("	")[1]))
            except: # line does not contain spectrum data
                continue
    
    return l, a