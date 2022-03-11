import numpy as np 

π = np.pi 

def au_temperature(x): ### given Temperature in Kelvin
    return x/315774.64

def MaxwellBoltzmann(Mass, Temperature):
    """ Get Maxwell-Boltzmann 3D velocities given: Mass (1D np.array, in a.u.) and Temperature (np.float, in Kelvin) """
    Temperature = au_temperature(Temperature)
    
    unique_masses, indices, counts = np.unique(Mass, return_inverse=True, return_counts=True)
    
    θ = np.random.random_sample(len(Mass))*(np.pi)
    φ = np.random.random_sample(len(Mass))*(2.*np.pi)

    if Temperature == 0:
        return np.zeros((len(Mass),3))
    
    else:
        speed_i = np.zeros(len(Mass))
        for i, masstype in enumerate(unique_masses):
            v_peak = np.sqrt( 2*Temperature / masstype )

            domain  = np.linspace(0., 5.*v_peak, 100)                                                           ### Speed values                                            
            PDF     = (2/np.pi)**(0.5) * (2/v_peak**2)**(1.5) * domain**2 * np.exp( - domain**2 / v_peak**2  )  ### Maxwell-Boltzmann PDF for given: Mass & Temperature
            PDF    *= 1/(np.sum(PDF))                                                                           ### Normalize to 1 for the specified range

            speeds     = np.random.choice(domain, counts[i], p=PDF)                                             ### counts[i] = Number of Atoms of a certain Mass type i
            j          = np.where(indices == i)[0]
            speed_i[j] = speeds

        vx = speed_i*np.sin(θ)*np.cos(φ)
        vy = speed_i*np.sin(θ)*np.sin(φ)
        vz = speed_i*np.cos(θ)
            
        return np.asarray([vx, vy, vz]).T


