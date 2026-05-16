import numpy as np
import params
import utilities as util
import matplotlib.pyplot as plt


util.print_intro()
util.clean_output()

# Arrays to store the computed values in the interacting case
mu_arr   = np.zeros((params.n_T))
N_arr    = np.zeros((params.n_T))
norm_Arr = np.zeros((params.n_T))
magn_arr = np.zeros((params.n_T))

# Arrays to store the computed values in the non-interacting case
magn_0_arr            = np.zeros((params.n_T))
magn_0_analytical_arr = np.zeros((params.n_T))

util.print_with_border("Doing the computations for different Temperatures")

for i, T in enumerate(params.T_arr):
    print(f"\n{'=' * 75}\n\nT = {T:.1f} K\n")
    mu_arr[i]   = util.get_CP           (T, interacting = True)
    N_arr[i]    = util.get_N_el         (T, mu_arr[i], interacting = True)
    norm_Arr[i] = util.get_norm         (T, mu_arr[i], interacting = True)
    magn_arr[i] = util.get_magnetization(T, mu_arr[i], interacting = True)    
    magn_0_arr[i]            = util.get_magnetization(T, params.mu_0, interacting = False) 
    magn_0_analytical_arr[i] = util.get_magnetization0_analytical(T)


print("\n\n")

# Plotting the DOS for a few temperatures
coarse_indices = np.linspace(0, len(params.T_arr) - 1, 3, dtype=int)
coarse_T_Arr   = params.T_arr[coarse_indices]
coarse_mu_arr  = mu_arr[coarse_indices]

util.plot_DOS_for_different_T(coarse_T_Arr, coarse_mu_arr)

plt.figure(1)
plt.plot(params.T_arr, mu_arr * 10 ** 3, label = r"$\mu(T)$", linewidth=2.5)
plt.axhline(y=0.0, color='green', linestyle='--', label=r'$ \mu_0 = 0 $', linewidth=2.5)
plt.xlabel(r"$T$ in Kelvin")
plt.ylabel(r'Chemical Potential $\mu(T)$ (meV)')
plt.legend()
plt.savefig(f"Output/mu_vs_T.png")
plt.show()

plt.figure(2)
plt.plot(params.T_arr, N_arr, label=r"$ \langle N_{el} \rangle (T) $", linewidth=2.5)
plt.plot(params.T_arr, norm_Arr, label=r"$ \langle N_{states} \rangle (T) $", linewidth=2.5)
plt.xlabel(r"$T$ in Kelvin")
plt.legend()
plt.savefig(f"Output/N_and_Norm_vs_T.png")
plt.show()

plt.figure(3)
plt.plot(params.T_arr, magn_arr / params.mu_z_0_gs, label=r"$ \langle \mu_z \rangle (T) /< g.s. | \mu_z | g.s. > $", linewidth=2.5)
plt.plot(params.T_arr, magn_0_arr / params.mu_z_0_gs, label=r"$ \langle \mu_z \rangle_0 (T) /< g.s. | \mu_z | g.s. > $", linewidth=2.5)
plt.plot(params.T_arr, magn_0_analytical_arr / params.mu_z_0_gs, '--', label=r"$ \langle \mu_z \rangle_{0} (T) /< g.s. | \mu_z | g.s. > = \tanh ( \epsilon_z / (2 k_B T) ) $", linewidth=2.5)
plt.ylabel(r'Normalized Magnetization along z-axis')
plt.xlabel(r"$T$ in Kelvin")
plt.ylim(0, 1.1)
plt.legend()
plt.savefig(f"Output/Magn_vs_T.png")
plt.show()

util.print_with_border("Program Finished Successfully")