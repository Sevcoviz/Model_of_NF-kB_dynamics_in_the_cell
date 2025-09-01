
import matplotlib.pyplot as plt

def plot_nfkb_dynamics(t, K, N, I, R, G, G_mm):
    """
    Vykreslí dynamiku NF-κB systému.
    
    Parameters:
    - t: časový vektor
    - K: koncentrácia IKKα
    - N: koncentrácia NF-κB
    - I: koncentrácia IκBα
    - R: koncentrácia IκBα mRNA
    - G: génový stav
    - G_mm: Michaelis–Menten aproximácia génového stavu
    """
    
    plt.figure(figsize=(7, 5))
    plt.plot(t, K, label='K (IKKα)', color='blue')
    plt.plot(t, N, label='N (NF-κB)', color='green')
    plt.plot(t, I, label='I (IκBα)', color='red')
    plt.plot(t, R, label='R (IκBα RNA)', color='purple')
    plt.plot(t, G, label='G (gene state)', color='brown')
    plt.plot(t, G_mm, '--', label='G (MM approx)', color='orange')
    plt.xlabel('Time [h]')
    plt.ylabel('Normalized Concentration')
    plt.title('NF-κB Signaling Dynamics')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    fig, axs = plt.subplots(2,2)
    axs[0,0].plot(t, K, label='K (IKKα)', color='blue')
    axs[0,0].set_title('IKKα Dynamics')

    axs[0,1].plot(t, N, label='N (NF-κB)', color='green')
    axs[0,1].plot(t, I, label='I (IκBα)', color='red')
    axs[0,1].set_title('NF-κB and IκBα Dynamics')

    axs[1,0].plot(t, R, label='R (IκBα RNA)', color='purple')
    axs[1,0].set_title('IκBα mRNA Dynamics')

    axs[1,1].plot(t, G, label='G (gene state)', color='brown')
    axs[1,1].plot(t, G_mm, '--', label='G (MM approx)', color='orange')


    for ax in axs.flat:
        ax.set_xlabel('Time [h]')
        ax.set_ylabel('Normalized Concentration')
        ax.legend()
        ax.grid(True)

    for ax in axs.flat:
        ax.label_outer()



# A single function to visualize the simulation results.
def plot_compartmental_dynamics(t, K, N_C, N_N, I_C, I_N, NI_C, NI_N, R, G):
    """
    Plots the dynamics of the compartmentalized NF-κB system.
    """
    plt.style.use('seaborn-v0_8-whitegrid')
    fig, axs = plt.subplots(3, 1, figsize=(10, 15), sharex=True)

    # Plot NF-κB dynamics
    axs[0].plot(t, N_C, label='NF-κB (Cytoplasm)', color='blue')
    axs[0].plot(t, N_N, label='NF-κB (Nucleus)', color='green')
    axs[0].set_title('NF-κB Dynamics')
    axs[0].set_ylabel('Normalized Concentration')
    axs[0].legend()

    # Plot IκB dynamics
    axs[1].plot(t, I_C, label='IκBα (Cytoplasm)', color='red')
    axs[1].plot(t, I_N, label='IκBα (Nucleus)', color='orange')
    axs[1].set_title('IκBα Dynamics')
    axs[1].set_ylabel('Normalized Concentration')
    axs[1].legend()

    # Plot Complex and mRNA dynamics
    axs[2].plot(t, NI_C, label='Complex (Cytoplasm)', color='purple')
    axs[2].plot(t, NI_N, label='Complex (Nucleus)', color='brown')
    axs[2].plot(t, R, label='IκBα mRNA', color='cyan')
    axs[2].set_title('Complex and mRNA Dynamics')
    axs[2].set_xlabel('Time [h]')
    axs[2].set_ylabel('Normalized Concentration')
    axs[2].legend()

    plt.tight_layout()
    plt.show()




def plot_nfkb_compartmental(t, K, N_C, N_N, I_C, I_N, NI_C, NI_N, G, R):
    """
    Vykreslí dynamiku compartmentalizovaného NF-κB systému.

    Parameters:
    - t: časový vektor (sekundy)
    - K: koncentrácia IKKα
    - N_C, N_N: koncentrácia NF-κB (cytoplazma, jadro)
    - I_C, I_N: koncentrácia IκBα (cytoplazma, jadro)
    - NI_C, NI_N: komplex NF-κB:IκBα (cytoplazma, jadro)
    - G: génový stav
    - R: IκBα mRNA
    """
    # Prevedieme čas na hodiny kvôli porovnaniu s publikovanými grafmi
    t_h = t / 3600.0

    # --- Kompletný prehľad ---
    plt.figure(figsize=(8,6))
    plt.plot(t_h, K, label='IKKα (K)', color='blue')
    plt.plot(t_h, N_C, label='NF-κB cyt (N_C)', color='green')
    plt.plot(t_h, N_N, label='NF-κB nuc (N_N)', color='lime', linestyle='--')
    plt.plot(t_h, I_C, label='IκBα cyt (I_C)', color='red')
    plt.plot(t_h, I_N, label='IκBα nuc (I_N)', color='darkred', linestyle='--')
    plt.plot(t_h, R, label='mRNA (R)', color='purple')
    plt.plot(t_h, G, label='Gene state (G)', color='brown')
    plt.xlabel("Time [h]")
    plt.ylabel("Normalized concentration")
    plt.title("Compartmental NF-κB dynamics (overview)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    # --- Subplots pre detail ---
    fig, axs = plt.subplots(2, 2, figsize=(10,7))

    # IKK dynamics
    axs[0,0].plot(t_h, K, color='blue', label='IKKα')
    axs[0,0].set_title("IKKα dynamics")

    # NF-κB dynamics
    axs[0,1].plot(t_h, N_C, color='green', label='NF-κB cyt')
    axs[0,1].plot(t_h, N_N, '--', color='lime', label='NF-κB nuc')
    axs[0,1].set_title("NF-κB dynamics")

    # IκBα dynamics
    axs[1,0].plot(t_h, I_C, color='red', label='IκBα cyt')
    axs[1,0].plot(t_h, I_N, '--', color='darkred', label='IκBα nuc')
    axs[1,0].set_title("IκBα dynamics")

    # Complex + gene/mRNA
    axs[1,1].plot(t_h, NI_C, color='orange', label='NI cyt')
    axs[1,1].plot(t_h, NI_N, '--', color='darkorange', label='NI nuc')
    axs[1,1].plot(t_h, R, color='purple', label='mRNA (R)')
    axs[1,1].plot(t_h, G, color='brown', label='Gene state (G)')
    axs[1,1].set_title("Complex + gene")

    # Formatting
    for ax in axs.flat:
        ax.set_xlabel("Time [h]")
        ax.set_ylabel("Normalized concentration")
        ax.grid(True)
        ax.legend()
        ax.label_outer()

    plt.tight_layout()
    plt.show()


def plot_compartmental_dynamics(t, K, N_C, N_N, I_C, I_N, NI_C, NI_N, R, G):
    """
    Vykreslí dynamiku NF-κB systému s oddelenými priestormi.
    Teraz zobrazuje aj celkové koncentrácie.
    """
    t_h = t / 3600 # convert to hours
    
    # Calculate total concentrations
    N_total = N_C + N_N
    I_total = I_C + I_N
    NI_total = NI_C + NI_N

    # --- Hlavný graf s celkovým prehľadom ---
    plt.figure(figsize=(10, 7))
    plt.plot(t_h, K, label='IKKα', color='blue')
    plt.plot(t_h, N_total, label='Total NF-κB', color='green', linewidth=2, linestyle='-')
    plt.plot(t_h, I_total, label='Total IκBα', color='red', linewidth=2, linestyle='-')
    plt.plot(t_h, NI_total, label='Total NF-κB:IκBα Complex', color='orange', linewidth=2, linestyle='-')
    plt.plot(t_h, R, label='IκBα mRNA (R)', color='purple')
    plt.plot(t_h, G, label='Gene state (G)', color='brown')
    plt.xlabel("Time [h]")
    plt.ylabel("Normalized concentration")
    plt.title("Compartmental NF-κB dynamics (Total concentrations)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    # --- Subplots pre detail ---
    fig, axs = plt.subplots(2, 2, figsize=(10, 7))
    fig.suptitle("Compartmental NF-κB dynamics (Detail)")

    # IKK dynamics
    axs[0, 0].plot(t_h, K, color='blue', label='IKKα')
    axs[0, 0].set_title("IKKα dynamics")
    axs[0, 0].set_xlabel("Time [h]")
    axs[0, 0].set_ylabel("Concentration")

    # NF-κB dynamics
    axs[0, 1].plot(t_h, N_C, color='green', label='NF-κB cyt')
    axs[0, 1].plot(t_h, N_N, '--', color='lime', label='NF-κB nuc')
    axs[0, 1].set_title("NF-κB dynamics")
    axs[0, 1].set_xlabel("Time [h]")
    axs[0, 1].set_ylabel("Concentration")
    axs[0, 1].legend()

    # IκBα dynamics
    axs[1, 0].plot(t_h, I_C, color='red', label='IκBα cyt')
    axs[1, 0].plot(t_h, I_N, '--', color='darkred', label='IκBα nuc')
    axs[1, 0].set_title("IκBα dynamics")
    axs[1, 0].set_xlabel("Time [h]")
    axs[1, 0].set_ylabel("Concentration")
    axs[1, 0].legend()

    # NF-κB:IκBα complex dynamics
    axs[1, 1].plot(t_h, NI_C, color='orange', label='Complex cyt')
    axs[1, 1].plot(t_h, NI_N, '--', color='brown', label='Complex nuc')
    axs[1, 1].set_title("Complex dynamics")
    axs[1, 1].set_xlabel("Time [h]")
    axs[1, 1].set_ylabel("Concentration")
    axs[1, 1].legend()

    plt.tight_layout()
    plt.show()

