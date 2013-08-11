from scipy import pi, vectorize

def inch2meter(inch):
    return (inch * 0.0254)
#inch2meter = vectorize(_inch2meter)

def meter2inch(meter):
    return (meter / 0.0254)
#meter2inch = vectorize(_meter2inch)

def meter2foot(meter):
    return (meter / 0.0254 / 12.0)
#meter2foot = vectorize(_meter2foot)

def foot2meter(foot):
    return (foot * 12.0 * 0.0254)
#foot2meter = vectorize(_foot2meter)

def lb2kgram(lb):
    return (lb * 0.45359)
#lb2kgram = vectorize(_lb2kgram)

def kgram2lb(kgram):
    return (kgram / 0.45359)
#kgram2lb = vectorize(_kgram2lb)

def newton2lb(newton):
    return (newton * 0.22480894)
#newton2lb = vectorize(_newton2lb)

def lb2newton(lbf):
    return (lbf * 4.4482216)
#newton2lb = vectorize(_newton2lb)

steel_density    = 7850.   # kg/m^3
aluminum_density = 2700.   # kg / m^3
pine_density     =  521.   # kg/m^3
air_density      =  1.225  # kg/m^3

pine_modulus     =   8.5e9 * 1.1  # eastern white pine modulus [Pa]
aluminum_modulus =  69e9          # [Pa]
steel_modulus    = 200e9          # [Pa]

