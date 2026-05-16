import numpy as np

# ===== #
# Tasks #
# ===== #

CP       = True  #             Plotting the [C]hemical [P]otential. (Is being computed by default.)
N_el     = True  # Computing & Plotting the [N]umber of [El]ectrons.
N_states = True  # Computing & Plotting the [N]umber of [States].
DOS      = True  # Computing & Plotting the [D]ensity [O]f [S]tates
Magn     = True  # Computing & Plotting the [Magn]etization.


# ==================================== #
# Variable Declarations for the System #
# ==================================== #

B           = 10.0                    # Magnetic Field in Tesla
Omega       = 10.0 * 10 ** (-3)       # Phonon Frequency in eV
g_int       = 5.0 * 10 ** (-3)        # Electron-Phonon Coupling Strength in eV
N_states    =  2.0                    # Number of Electronic Single-Particle States in the System

N_el        =  1.0                             # Number of electrons that we seek to have in the system
n_T         = 10                               # Number of Temperature Points

# For bette resolution at lower temperatures: 
# We create a temperature array that is more dense at low T and more coarse at high T

T_min       =  10.0 ** (-4)                    # Minimum Temperature in Kelvin
T_mid       =   2.0                            # Midpoint Temperature in Kelvin
T_max       =   5.0                            # Maximum Temperature in Kelvin

T1_arr      = np.linspace(T_min, T_mid, n_T // 2)          # First  Half of the Temperature Array
T2_arr      = np.linspace(T_mid, T_max, n_T // 2 + 1)[1:]  # Second Half of the Temperature Array
T_arr       = np.concatenate((T1_arr, T2_arr))             # Full               Temperature Array



# ====================== #
# Convergence Parameters #
# ====================== #

n_Energy = 6000              # Number of Energy Points for the Integration
E_max    = 30 * 10 ** (-3)   # Maximum Energy in eV for the Integration
eta      = 1.0 * 10 ** (-5)  # Regularization Parameter in meV for the retarded Green-function


# ================== #
# Physical Constants #
# ================== #

k_B     =  8.617 * 10 ** (-5)          # in eV / K
hbar    =  6.582119 * 10 ** (-16)      # in eV s
q_e     = -1.602176634 * 10 ** (-19)   # charge of an electron in Coulomb.
m       =  9.1093837015 * 10 ** (-31)  # mass of an electron in kg
g       =  2.00231930436256            # g-factor of an electron
gamma_e = - g * np.abs(q_e) / (2 * m)  # gyromagnetic ratio in Coulomb / kg


# ==================================== #
# Other variables for the Computations #
# ==================================== #

mu_0      = 0.0                            # Chemical Potential for the non-interacting system (temperature-independent)
epsilon_z = np.abs(gamma_e * hbar * B / 2) # Zeeman-Energies in eV
mu_z_0_gs = - gamma_e * hbar / 2           # Ground-State Magnetization along z-axis 