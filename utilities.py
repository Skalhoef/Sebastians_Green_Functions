import os
import shutil
import params
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate


def print_with_border(message):
    print("+" + "=" * (len(message) + 4) + "+")
    print("|| " + message + " ||")
    print("+" + "=" * (len(message) + 4) + "+")
    print()

def print_intro():
    print("\n" + "=" * 75)
    print("=" * 15 + " Electron-Phonon-System with Green-Functions " + "=" * 15)
    print("=" * 75 + "\n\n")
    
    # print("We consider the model Hamiltonians H = H_el + H_ph + H_int with\n")
    print("We consider the model Hamiltonians\n")
    print("    H = H_el + H_ph + H_int\n")
    print("with\n")

    print(" - H_el  = epsilon_z * (n_up - n_down)")
    print(" - H_ph  =     Omega * b^dagger * b\n")
    print(" - H_int = 0                                        (non-interacting case)")
    print(" - H_int = g_int * (n_up + n_down) + (b^dagger + b) (    interacting case)\n")

    print_with_border("Tasks")

    if params.CP:       print(" -             Plotting the Chemical Potential mu(T)")
    if params.N_el:     print(" - Computing & Plotting the Number of Electrons <N_el>(T) (should be constant)")
    if params.N_states: print(" - Computing & Plotting the Number of States N_states(T) (should be constant)")
    if params.DOS:      print(" - Computing & Plotting the Density Of States DOS(T)")
    if params.Magn:     print(" - Computing & Plotting the Magnetization M_z(T)")
    print("\n")

    print_with_border("Model Parameters")

    print(f" - Magnetic Field Strength           B = {params.B:.3f} Tesla")
    print(f" - Zeeman-Energies         +-epsilon_z = +- {params.epsilon_z * 10 ** 3:.3f} meV")
    print(f" - Phonon Energy                 Omega = {params.Omega * 10 ** 3} meV")
    print(f" - Electron-Phonon Coupling:     g_int = {params.g_int * 10 ** 3} meV\n\n")


    print_with_border("Convergence Parameters")

    print(f" - Regularization-Parameter        eta = " + str(params.eta * 10**3) + " meV")
    print(f" - Interval for the Integrations:    I = [{params.E_max * 10 ** 3}, {params.E_max * 10 ** 3}] meV\n\n")

def clean_output():

    directories = ["Plots"]

    for dir_name in directories:
        dir_path = os.path.join(os.getcwd(), dir_name)
        if os.path.exists(dir_path) and os.path.isdir(dir_path):
            for filename in os.listdir(dir_path):
                file_path = os.path.join(dir_path, filename)
                try:
                    if os.path.isfile(file_path) or os.path.islink(file_path):
                        os.unlink(file_path)
                    elif os.path.isdir(file_path):
                        shutil.rmtree(file_path)
                except Exception as e:
                    print(f"Failed to delete {file_path}. Reason: {e}")

def Bose(E, T):
    beta = 1 / (params.k_B * T)
    np.seterr(over='ignore')
    return 1 / (np.exp(beta * np.real(E))- 1)

def Fermi(E, T):
    beta = 1 / (params.k_B * T)
    np.seterr(over='ignore')
    return 1 / (np.exp(beta * np.real(E)) + 1)

def integration(integrand):
    I, Error = integrate.quad(integrand, a=- params.E_max, b=+ params.E_max, limit=params.n_Energy)
    return I

def multi_purpose_function(T, mu, Type, E, interacting = True):
    
    # E_xi = hbar * omega_xi for xi in {p, m}
    E_p = + params.epsilon_z - mu
    E_m = - params.epsilon_z - mu


    if not interacting:
        Sigma = 0.0 * np.eye(2)
    else:

        # Hartree-Contribution to the Self-Energy
        D_00    = - 2 / params.Omega
        Sigma_H = params.g_int ** 2 * D_00 * (Fermi(E_p, T) + Fermi(E_m, T)) * np.eye(2)

        # Fock-Contribution to the Self-Energy
        # [O]ccupation-[P]ole-[T]erm
        OPT_1_up = (Bose(params.Omega, T) + Fermi(E_p, T)) / (E + 1j * params.eta - E_p + params.Omega)
        OPT_2_up = (Bose(params.Omega, T) + 1 - Fermi(E_p, T)) / (E + 1j * params.eta - E_p - params.Omega)
        M_up     = OPT_1_up + OPT_2_up


        OPT_1_down = (Bose(params.Omega, T) + Fermi(E_m, T)) / (  E + 1j * params.eta - E_m + params.Omega)
        OPT_2_down = (Bose(params.Omega, T) + 1 - Fermi(E_m, T)) / (E + 1j * params.eta - E_m - params.Omega)
        M_down     = OPT_1_down + OPT_2_down

        M       = np.diag([M_up, M_down])
        Sigma_F = params.g_int ** 2 * M

        # Total Self-Energy
        Sigma   = Sigma_H + Sigma_F


    H_0 = np.diag([E_p, E_m])
    Xi  = H_0 + Sigma - mu * np.eye(2) 

    # Calculation of the Green-function
    G = np.linalg.inv( (E + 1j * params.eta) * np.eye(2) - Xi)

    if Type == "dN_over_dE":
        dN_over_dE = - 1 / np.pi * Fermi(E, T) * np.imag(np.trace(G))
        return dN_over_dE

    elif Type == "dNstates_over_dE":
        dNstates_over_dE = - 1 / np.pi * np.imag(np.trace(G))
        return dNstates_over_dE
    
    elif Type == "dMag_over_dE":
        mag_up = 1 / (2 * np.pi) * np.imag(G[0, 0]) * Fermi(E, T) * params.g * np.abs(params.q_e) * params.hbar / (2 * params.m)
        mag_down = 1 / (2 * np.pi) * np.imag(G[1, 1]) * (-1) * Fermi(E, T) * params.g * np.abs(params.q_e) * params.hbar / (2 * params.m)
        dMag_over_dE = mag_up + mag_down
        return dMag_over_dE

def find_root_monotonic(func, low, high, tolerance=1e-12, max_iterations=100):
    """
        Find the root of a monotonically increasing function within the interval [low, high].

        Parameters:
        - func: The monotonically increasing function.
        - low: The lower bound of the search interval.
        - high: The upper bound of the search interval.
        - tolerance: The desired accuracy of the root.
        - max_iterations: The maximum number of iterations.

        Returns:
        - root: The approximate root of the function within the specified interval.
        """

    print(" - Searching for the root of a monotonically increasing function.")
    print(" - The interval in which we expect the root is J = [" + str(low * 10 ** 3) + ", " + str(high * 10 ** 3) + "] meV.")

    # Check if the function values at the interval bounds have opposite signs
    if func(low) * func(high) > 0:
        print_with_border("ERROR")
        print("\n ====> The binary-search-algorithm is not applicable.")
        print(" ====> The function values at the interval bounds must have opposite signs.")
        print(f" ====> We have f(low = {low * 10 ** 3:.1f} meV) = {func(low):.3f} and f(high = {high * 10 ** 3:.1f} meV) = {func(high):.3f}\n")
        return 0
        # raise ValueError("Function values at interval bounds must have opposite signs.")

    # Perform binary search
    for _ in range(max_iterations):
        mid = (low + high) / 2
        f_mid = func(mid)
        if _ > max_iterations / 2:
            print("mu = " + str(mid) + " would give R(mu) = " + str(f_mid))

        # Check if the root is found within the tolerance
        if abs(f_mid) < tolerance:
            return mid

        # Update the search interval
        if f_mid > 0 and low != mid:
            low = mid
        elif f_mid < 0 and high != mid:
            high = mid
        else:
            print("The binary-search-algorithm does not converge.")
            return mid

    # If the maximum number of iterations is reached, raise an exception
    print("Binary search for the root didn't converge. (Should be further analysed.)")
    return mid
    # raise RuntimeError("Binary search did not converge within the specified number of iterations.")

def get_CP(T, interacting):

    def zero(mu):

        def dNoverdE(E):
            return multi_purpose_function(T, mu, "dN_over_dE", E, interacting)
        
        N_calc = integration(dNoverdE)

        return params.N_el - N_calc

    mu = find_root_monotonic(zero, - params.E_max / 5, params.E_max / 5)

    print(f" - The determined value mu = {(mu * 10 ** 3):.3f} meV yields <N>(mu) = {(zero(mu) + params.N_el):.3f}")

    return mu

def get_N_el(T, mu, interacting):

    def integrand(E):
        return multi_purpose_function(T, mu, "dN_over_dE", E, interacting)

    I = integration(integrand)

    return I

def get_norm(T, mu, interacting):

    def integrand(E):
        return multi_purpose_function(T, mu, "dNstates_over_dE", E, interacting)

    I = integration(integrand)
    print(f" - The normalization condition is fulfilled if N_states = 2 = {I:.3f}")

    return I

def get_magnetization(T, mu, interacting):

    def integrand(E):
        return multi_purpose_function(T, mu, "dMag_over_dE", E, interacting)

    I = integration(integrand)
    if interacting: print(f" - The Magnetization along the z-axis for the interacting system reads \n    <\mu_z> / < g.s. | \mu_z | g.s. > = {I / params.mu_z_0_gs:.3f}")
    if not interacting: print(f" - The Magnetization along the z-axis for the non-interacting system reads \n    <\mu_z>_0 / < g.s. | \mu_z | g.s. > = {I / params.mu_z_0_gs:.3f}")

    return I

def get_magnetization0_analytical(T):

    I = np.tanh(params.epsilon_z / (2 * params.k_B * T)) / 2.0 * params.g * np.abs(params.q_e) * params.hbar / (2 * params.m)
    print(f" - The Magnetization (analytical expression) along the z-axis for the non-interacting system reads \n    <\mu_z>_0 / < g.s. | \mu_z | g.s. > = {I / params.mu_z_0_gs:.3f} ")

    return I

def plot_DOS_for_different_T(Coarse_T_Arr, coarse_mu_arr):
        
    print(" - Plotting the Density of States")
    
    plt.figure(10)
    
    # Non-Interacting DOS is just a sum of two delta-functions
    plt.axvline(x=- params.epsilon_z * 10 ** 3, color='red', linestyle='--', label=r'$\rho_0 (E)$', linewidth=2.5)
    plt.axvline(x=+ params.epsilon_z * 10 ** 3, color='red', linestyle='--', linewidth=2.5)
    
    plt.xticks(ticks=[-params.epsilon_z * 10 ** 3, params.epsilon_z * 10 ** 3],
               labels=[r'$- \epsilon_z$', r'$+ \epsilon_z$'])
    

    for i, T in enumerate(Coarse_T_Arr):

        def dNstatesoverdE(E):
            return multi_purpose_function(T, coarse_mu_arr[i], "dNstates_over_dE", E, interacting = True)
        
        E_Arr = np.linspace( - 3 * params.epsilon_z , 3 * params.epsilon_z, params.n_Energy)
        dNoverdE_Arr = [dNstatesoverdE(E) for E in E_Arr]
        
        # Convert to numpy array
        dNoverdE_Arr = np.array(dNoverdE_Arr)

        plt.plot(E_Arr * 1e3, dNoverdE_Arr * 1e-3, label=rf"$\rho(E, T = {T:.1f}\,\mathrm{{K}})$")
    
    
    plt.xlim(-2.0, 2.0)
    plt.ylim(bottom=0)
    plt.xlabel(r"$E$ $(meV)$")
    plt.ylabel(r"$ \rho (E ,T) $ $(meV^{-1})$")
    plt.title(r"Density of States")
    plt.legend()
    plt.savefig(f"Output/DOS.png")
    plt.show()
    plt.clf()  # Clear the figure for the next iteration
    print()


