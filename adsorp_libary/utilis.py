import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline

def get_gas_constant(gas_type='N2'):
    """Get gas constant of two different gas

    Args:
        gas_type: str(), the type of the gas
    Returns:
        const: dict(), contain the A and Vmol 
    """

    const = dict()
    # liquid molar volume [cm^3/mol]
    if gas_type == 'N2':
        const['A'] = 9.53
        const['Vmol'] = 34.67
    elif gas_type == 'Ar':
        const['A'] = 10.44  # 10.44(best)
        const['Vmol'] = 28 #28(best), 22.56(calculated) 
    else:
        print('wrong gas type')
        const=-1
    return const

def use_my_pressure_points(p_exp, Q_exp, gas_type = 'N2'):
    """Define my own radius points, use raw pressure and adosorption quantity,
        to interpolate the correponding pressure and adsorption quantity

    Args:
        p_exp: 
    Returns:

    """
    const = get_gas_constant(gas_type)
    # predefine the radius
    radius_s = np.array([8,9,10,11,12,13,14,16,19,27,42,75,140,279,360,460,828,1600])
    # convert radius to pressure
    p_s = radius_to_pressure(radius_s,const)
    func_spline = get_spline(p_exp,Q_exp,2)
    # get Q_s from spline function
    Q_s = func_spline(p_s)
	
    return p_s,Q_s
	
def get_spline(p_rels,Q,order=2):
    """Calculate spline base function using p_rels and Q

    Args:
        p_rels: numpy 1D array, constain a list of relative pressure, each element has no unit  
        Q: numpy 1D array, constain a list of adsorption quantity responding to p_rels
    """
    s = InterpolatedUnivariateSpline(p_rels,Q,k=order)
    return s
	
def insert_zero(a):
    """ Insert zero in 1D array, so to make the starting points start from 1 

    Args:
        a: np.array()
    Returns:
        np.array()
    """
    return np.insert(a,0,0)
	
def restrict_isotherm( P, Q, Pmin, Pmax ):
    """Restrict the isotherm Q, P to pressures between min and max.

    Q: Quanity adsorbed
    P: Pressure (relative or absolute)
    Pmin: minimum pressure
    Pmax: maximum pressure

    Returns:  Qads, P restricted to the specified range
    """

    index = (P >= Pmin) & (P <= Pmax)
    #b = np.logical_and( P >= Pmin, P <= Pmax)
    return P[index], Q[index]
	
def kelvin_radius(p_rel,const):
    """ Calculate kelvin radius from pressure 

    Args:
        p_rel: float, relative pressure, no unit, value from 0 to 1
        const: dict, gas constant, 

    Returns:
        kelvin radius, unit [A] 
    """
    # Rc is in [A]
    # np.log is e base
    return -const['A'] / np.log(p_rel) 

def radius_to_pressure(Rc,const):	
    """Calculate radius to pressure

    Args:
        Rc: float, the radius of the pore, in unit of [A]
        const: dict(), gas constant, keys contains ['A']
    Return: 
        pressure, in unit 

    """

    #Rc is in [A]
    return np.exp(-const['A']/Rc)
	
def thickness_Harkins_Jura(p_rel):
    """
    Args:
        p_rel: float, relative pressure, no unit,

    Returns:
        thickness, in unit of [A]
    """

    return   (13.99 / (0.034 - np.log10(p_rel)))**0.5