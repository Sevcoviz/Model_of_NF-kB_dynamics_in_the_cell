# analysis_tools.py

import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.integrate import solve_ivp

# Importujeme triedy z našich existujúcich súborov
from CytoNuc_rovnice import NFkBSystemExact
from CytoNuc_params import CytoNucParamsExact

def run_bifurcation_analysis(observed_variable, bifurcation_param, param_range, y0, t_span, t_eval):
    """
    Vykoná všeobecnú bifurkačnú analýzu pre kompartmentalizovaný model.

    Args:
        observed_variable (str): Názov premennej, ktorá sa má sledovať na osi y (napr. 'Nn', 'Im', 'I').
        bifurcation_param (str): Názov parametra, ktorý sa má meniť na osi x (napr. 't3', 'k1').
        param_range (np.ndarray): Pole hodnôt pre bifurkačný parameter.
        y0 (list): Počiatočné podmienky pre systém.
        t_span (tuple): Časový interval simulácie.
        t_eval (np.ndarray): Časové body pre vyhodnotenie riešenia.
    
    Returns:
        tuple: Vráti dáta (param_values, min_values, max_values) pre prípadné ďalšie spracovanie.
    """
    print(f"Spúšťam bifurkačnú analýzu pre parameter '{bifurcation_param}', sledujem premennú '{observed_variable}'...")

    # Mapa názvov premenných na ich indexy v poli `sol.y`
    variable_map = {'N': 0, 'Nn': 1, 'I': 2, 'In': 3, 'Im': 4, 'NI': 5, 'NIn': 6}
    
    if observed_variable not in variable_map:
        raise ValueError(f"Neznáma premenná '{observed_variable}'. Dostupné možnosti: {list(variable_map.keys())}")
    
    observed_idx = variable_map[observed_variable]

    param_values, min_values, max_values = [], [], []

    for p_val in tqdm(param_range, desc=f"Analyzujem {bifurcation_param}"):
        params = CytoNucParamsExact()
        setattr(params, bifurcation_param, p_val)
        system = NFkBSystemExact(params)
        
        sol = solve_ivp(fun=system.rhs, t_span=t_span, y0=y0, t_eval=t_eval, method="LSODA")
        
        if not sol.success:
            continue

        tail_fraction = 0.5
        num_points_in_tail = int(len(sol.t) * tail_fraction)
        variable_tail = sol.y[observed_idx, -num_points_in_tail:]
        
        min_val = np.min(variable_tail)
        max_val = np.max(variable_tail)
        
        param_values.append(p_val)
        min_values.append(min_val)
        max_values.append(max_val)

    # Vykreslenie diagramu
    plt.figure(figsize=(12, 7))
    plt.plot(param_values, min_values, 'k.', markersize=2)
    plt.plot(param_values, max_values, 'k.', markersize=2)
    
    plt.title(f'Bifurkačný Diagram: Vplyv "{bifurcation_param}" na "{observed_variable}"', fontsize=16)
    plt.xlabel(f'Hodnota parametra "{bifurcation_param}"')
    plt.ylabel(f'Ustálená koncentrácia {observed_variable} [μM]')
    plt.grid(True)
    plt.show()

    return param_values, min_values, max_values