# Provide three of the four parameters (a, e(L mol-1 cm-1), c(mol/L), l(m)) to calculate the missing one for the Beer-Lambert law.
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

# Provide two of the three parameters (m(g), M(g/mol), mol) to calculate the missing one for molar calculations.
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

# Provide the three parameters (mol, v(L), c(mol/L)) to calculate the missing one for solution calculations.
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

# Provide the four parameters (M(g/mol), mstock(g), vstock(L), cstock(mol/L), valiquot(L), vsol(L), csol(mol/L)) to calculate the missing one for solution preparation calculations.
def solprop(M=None, mstock=None, vstock=None, cstock=None, valiquot=None, vsol=None, csol=None):
    try:
        if cstock is None:    
            cstock = solcalc(mol=molcalc(m=mstock, M=M), v=vstock) # concentration of stock solution
        
        if csol is None:
            return solcalc(mol=solcalc(c=cstock, v=valiquot), v=vsol) # concentration of overall solution
        elif vsol is None:
            return solcalc(mol=solcalc(c=cstock, v=valiquot), c=csol) # volume in overall solution
        elif valiquot is None:
            return solcalc(mol=solcalc(c=csol, v=vsol), c=cstock) # volume of aliquot taken from stock solution
    except:
        print("Error: four of the seven parameters (M(g/mol), mstock(g), vstock(L), cstock(mol/L), valiquot(L), vsol(L), csol(mol/L)) needed to calculate the missing one.")
        return None