import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from sklearn.metrics import r2_score
import matplotlib
matplotlib.use('Agg')
plt.rcParams['font.sans-serif'] = ['avenir']
plt.rcParams['font.size'] = 14

# Predict modulus function
def predict_modulus(K_val, conc, other_K, init_g, comp_concs, fit_target):
    if fit_target == "Kac":
        Kac = K_val[0]
        Kab = other_K
    else:  # fitting Kab
        Kac = other_K
        Kab = K_val[0]

    Kab_app = Kab / (1 + Kac * comp_concs)
    p = (1 + (1 / (2 * conc_cross * Kab_app))) - np.sqrt((1 + (1 / (2 * conc_cross * Kab_app)))**2 - 1)
    g0 = (conc_cross / 16) * (3 - np.sqrt((4 / p) - 3))**3 * (np.sqrt((4 / p) - 3) + 1)
    g0_init = (conc_cross / 16) * (3 - np.sqrt((4 / p[0]) - 3))**3 * (np.sqrt((4 / p[0]) - 3) + 1)
    normalized_g0 = g0 / g0_init
    return init_g * normalized_g0

# Loss function
def loss(K_val, conc_cross, other_K, init_g, comp_concs, exp_moduli, fit_target):
    pred = predict_modulus(K_val, conc_cross, other_K, init_g, comp_concs, fit_target)
    return np.sum((pred - exp_moduli)**2)

# Streamlit UI
st.title("Estimate Keq from Modulus")
st.write("Fit either the competitor Keq (Kac) or the crosslink Keq (Kab) from modulus data.")

# Choose parameter to fit
fit_target = st.radio("Which parameter do you want to estimate?", ["Kac", "Kab"])

# Inputs
init_g = st.number_input("Initial modulus (G₀)", value=21.8539333)
conc_cross_mM = st.number_input("Crosslink concentration (mM)", value=80.0)

if fit_target == "Kac":
    fixed_Kab = st.number_input("Kab (Crosslink Keq, M^-1)", value=2185.0)
    guess = 1000
    fixed_param = fixed_Kab
else:
    fixed_Kac = st.number_input("Kac (Competitor Keq, M^-1)", value=1000.0)
    guess = 2185
    fixed_param = fixed_Kac

# Experimental data input
st.markdown("### Experimental Data")
default_concs = "0, 10, 20, 30, 40"
default_moduli = "21.8539333, 17.052, 14.9795, 12.54093333, 11.6045"
comp_concs_str = st.text_input("Competitor concentrations (comma-separated, in mM)", default_concs)
exp_moduli_str = st.text_input("Observed moduli (comma-separated)", default_moduli)

try:
    # Convert inputs
    comp_concs_mM = np.array([float(val.strip()) for val in comp_concs_str.split(',')])
    exp_moduli = np.array([float(val.strip()) for val in exp_moduli_str.split(',')])
    comp_concs = comp_concs_mM / 1000
    conc_cross = conc_cross_mM / 1000

    # Fit selected parameter
    result = minimize(
        loss,
        x0=[guess],
        args=(conc_cross, fixed_param, init_g, comp_concs, exp_moduli, fit_target),
        bounds=[(1, 1e25)]
    )

    fit_value = result.x[0]
    fit_vals = predict_modulus([fit_value], conc_cross, fixed_param, init_g, comp_concs, fit_target)
    r2 = r2_score(exp_moduli, fit_vals)

    st.success(f"Estimated {fit_target} (M^-1): {fit_value:.2f}")
    st.markdown(f"**R² of fit:** {r2:.3f}")

    # Smooth fit
    smooth_concs_mM = np.linspace(0, max(comp_concs_mM)*1.1, 1000)
    smooth_concs = smooth_concs_mM / 1000
    smooth_fit = predict_modulus([fit_value], conc_cross, fixed_param, init_g, smooth_concs, fit_target)

    # Plotting
    fig, ax = plt.subplots()
    ax.plot(comp_concs_mM, exp_moduli, 'o', label='Experimental', color='black')
    ax.plot(smooth_concs_mM, smooth_fit, '-', label='Model Fit', color='blue')
    ax.set_xlabel("Competitor concentration (mM)")
    ax.set_ylabel("Modulus")
    textstr = f"{fit_target} = {fit_value:.2e} M^-1\nR² = {r2:.3f}"
    ax.text(0.53, 0.75, textstr, transform=ax.transAxes,
        fontsize=12, verticalalignment='top',
        bbox=dict(boxstyle="round,pad=0.3", facecolor='white', alpha=0.8))
    ax.legend()
    st.pyplot(fig)

except Exception as e:
    st.error(f"Something went wrong: {e}")
