
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