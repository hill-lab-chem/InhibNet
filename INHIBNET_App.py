import numpy as np
import streamlit as st
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.optimize import minimize
import plotly.graph_objects as go
import pandas as pd
from sklearn.metrics import r2_score
from PIL import Image, ImageDraw
import random
# ========================================
# Description / Welcome Page
# ========================================
def description_page():
    st.title("InhibNet (Inhibition of Polymer Networks)")
    
    st.markdown("""
    ### Welcome!
    This app allows you to explore how **competitive inhibition** affects mechanical properties of polymer networks.

    For more information about this tool feel free to read our associated publication here: **TBD Please keep your eyes peeled for updates!**

    Or for information about how the code works feel free to visit our GitHub here: https://github.com/hill-lab-chem/InhibNet/

    The purpose of this work is to quantify how a polymer network's material properties can be altered by adding a small molecule that can compete with a crosslink. It has long been shown in the literature 
    that the material properties (i.e. modulus and relaxation time) of dynamic polymer networks are determined by the dynamics of the crosslink. We propose that adding a small molecule competitor can 
    alter the ma
    We hypothesized that principles from competitive inhibition of enzymes could be adapted to dynamic hydrogels to provide a similarly simple framework for predicting how key network
    properties change in the presence of competing species. In particular, we reasoned that the apparent equilibrium ($K_{a,app}$), widely used in enzyme kinetics to capture effective affinities
    under competitive inhibition, could be directly translated to dynamic networks as a predictor of mechanical response. The crosslinks exist in a ternary equilibrium between unbound,
    crosslinked, and the crosslink-competitor complex, so we reason that the association constants for the formation of the crosslink ($K_{a,XL}$) and the association of the competitor with the
    crosslink ($K_{a,C}$) could be inputted to the $K_{a,app}$ assumption.""")

    st.image("pics/TernaryEq.jpg", 
         
         use_container_width=True)
    st.markdown("**Figure:** Ternary Equilibrium, allowing for the $K_{a,app}$ assumption")
    st.markdown("""



    We designed a quantitative model that inputs $K_{a,app}$ into crosslink conversion ($p$, eqn. 1) which can be used to further calculate the shear modulus from either
    the **affine** or **phantom** network model (eqns. 2,3). The phantom network captures defects well, especially for dilute systems, whereas the affine better depicts
    the material's properties at high concentrations.""")
    st.markdown(r"""
    ### Equations

    **1. Crosslink conversion:**

    $$
    p = \left(1 + \frac{1}{2 N_a K_{a,app}}\right) - \sqrt{\left(1 + \frac{1}{2 N_a K_{a,app}}\right)^2 - 1}
    $$

    **2. Shear modulus (affine network):**
   
    First, define the effective probability term:
    $$
    P_{out} = \sqrt{\frac{1}{p} - \frac{3}{4}} - \frac{1}{2}
    $$
    
    The let us calculate the probability 3 crosslinks will form:

    $$
    P_{3} = 4*P_{out} (1-P_{out})^3
    $$
    Then the probability 4 crosslinks will form:

    $$
    P_{4} = (1-P_{out})^4
    $$
    So the elastically active strands are:
    $$
    v_{e} = \frac{N_a}{4} (\frac{3}{2}P_3+2P_4 )     
    $$
    So modulus is:
    $$
    \frac{g_0}{K_b T} = v_{e}
    $$
    
    **2. Shear modulus (phantom network):**
    The phantom network corrects for concentration of crosslinks as:
    $$
    \frac{g_0}{K_b T} = v_{e} - \mu_{e}
    $$
    where:
    $$
    \mu_{e} = \frac{N_a}{4} (P_3+P_4 )   
    $$
    
    Where the modulus equation can be simpliefied as:
    $$
    g_0 = \frac{N_a}{16} \left(3 - \sqrt{\frac{4}{p} - 3}\right)^3 \left(\sqrt{\frac{4}{p} - 3} + 1\right)
    $$
    """)

    st.markdown("""
    **Four main tools are included:**
    1. **Modulus vs Competitor Concentration** – 2D plot to see how modulus changes with a single competitor. Input features of the crosslink, and the competitor, and the initial stiffness of your gel.
    The output will be a prediction of how modulus will change under a range of competitor concentrations. Feel free to download a CSV of your data.
    2. **Experimental Modulus and Concentration to predict $K_{a,C}$ or $K_{a,XL}$** – A tool to predict either $K_{a,C}$ or $K_{a,XL}$ from experimental modulus and concentration data. 
    3. **Relaxation time as a function of competitor** – This predicts how relaxation time changes as a function of competitor. Simply input the uninhibited relaxation time, 
    and thermodynamic information about the crosslink and competitor and outputs for relaxation time will be plotted.
    4. **Modulus Surface Visualization** – This represents the design space available in this work. You can input features like assocation of the crosslink($K_{a,XL}$),
    concentration of crosslinks(Na), association of the competitor($K_{a,C}$), concentration of competitor([C]). This plots modulus vs. $K_{a,C}$ vs. [C]. This interactive graph allows
    for you to see what modulus values a range of comeptitor strengths and concentrations will yield.
    
    ### Network Models
    You can choose between two network models:
    - **Phantom Network** – accounts for network fluctuations, this works best for dilute polymer networks (i.e. close to the overlap concentration).
    - **Affine Network** – assumes a fully connected network with affine deformations, this works best for more concentrated networks (i.e. greater than ~4x overlap concentration).
    
    ### Instructions
    - Use the sidebar to select a page.
    - Adjust the model parameters in each app.
    - For App 4, hover over the 3D plot to see modulus values at each point.
    - Use the model selector to switch between Phantom and Affine calculations.
       
    """)
# ========================================
# Phantom Network Functions
# ========================================
def phantom_network_modulus_surface(conc_cross, Kab, Kac_range, comp_conc_range, g0_init_user): #this is for the 3d plot of the phantom network model. inputs are: concentration of crosslinks, assocation of the crosslink, range of comeptitor assocaition concnstants, range of competitor concentrations, and the intial modulus.
    Kac_vals = np.linspace(Kac_range[0], Kac_range[1], 100) #this creates an evenly spaced list of 100 competitor association constant values
    comp_concs = np.linspace(comp_conc_range[0], comp_conc_range[1], 100) #this creates an evenly spaced list of 100 competitor concentration values
    CCOMP, KAC = np.meshgrid(comp_concs, Kac_vals) #this is a little matrix of the possible concentrations with association constants.

    Kab_app = Kab / (1 + KAC * (CCOMP / 1000)) #this calculates the Ka,app assumption from enzyme kinetics.
    inner = 1 / (2 * conc_cross / 1000 * Kab_app) #this is the inner term of a sqr root for calculating conversion of crosslinks
    p = (1 + inner) - np.sqrt((1 + inner)**2 - 1) #this calculates conversion of crosslinks

    sqrt_term = (4 / p) - 3 #this is a term that will be inside the square root for the shear modulus
    sqrt_term[sqrt_term < 0] = np.nan #safety to make sure we don't have negative values
    root = np.sqrt(sqrt_term) #square root

    g0 = (conc_cross / 1000 / 16) * (3 - root)**3 * (root + 1) #shear modulus
    g0[g0 < 0] = np.nan #drop any negative values because they represent no gel formed.

    g0_init = g0[0, 0] #initial modulus
    normalized_g0 = g0 / g0_init #normalizing modulus to initial value
    real_g0 = normalized_g0 * g0_init_user #by multiplying by the initial modulus value you get the output based on your initial system.

    return comp_concs, Kac_vals, real_g0, CCOMP, KAC

def phantom_2d_line(conc_cross, Kab, Kac, comp_conc_range, g0_init_user): #2d version of the 3d one
    comp_concs = np.linspace(comp_conc_range[0], comp_conc_range[1], 100)
    

    Kab_app = Kab / (1 + Kac * (comp_concs / 1000))
    inner = 1 / (2 * conc_cross / 1000 * Kab_app)
    p = (1 + inner) - np.sqrt((1 + inner)**2 - 1)

    sqrt_term = (4 / p) - 3
    sqrt_term[sqrt_term < 0] = np.nan
    root = np.sqrt(sqrt_term)

    g0 = (conc_cross / 1000 / 16) * (3 - root)**3 * (root + 1)
    g0[g0 < 0] = np.nan

    g0_init = g0[0]
    normalized_g0 = g0 / g0_init
    real_g0 = normalized_g0 * g0_init_user

    return comp_concs, normalized_g0, real_g0

def predict_phantom_modulus(K_val, conc_cross, other_K, init_g, comp_concs, fit_target): #for fitting Kac or Kab
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

# ========================================
# Affine Network Functions
# ========================================


def affine_network_modulus_surface(conc_cross, Kab, Kac_range, comp_conc_range, g0_init_user):  #this is for the 3d plot of the affine network model. inputs are: concentration of crosslinks, assocation of the crosslink, range of comeptitor assocaition concnstants, range of competitor concentrations, and the intial modulus.
    Kac_vals = np.linspace(Kac_range[0], Kac_range[1], 100) #this creates an evenly spaced list of 100 competitor association constant values
    comp_concs = np.linspace(comp_conc_range[0], comp_conc_range[1], 100)#this creates an evenly spaced list of 100 competitor concentration values
    CCOMP, KAC = np.meshgrid(comp_concs, Kac_vals)#this is a little matrix of the possible concentrations with association constants.

    Kab_app = Kab / (1 + KAC * (CCOMP / 1000))#this calculates the Ka,app assumption from enzyme kinetics.
    p = (1 + 1 / (2 * conc_cross / 1000 * Kab_app)) - np.sqrt((1 + 1 / (2 * conc_cross / 1000 * Kab_app))**2 - 1) #this calculates conversion of crosslinks

    P_out = np.sqrt(1 / p - 0.75) - 0.5 #for calculating probability of 3 and 4 crosslinks forming
    P3 = 4 * P_out * (1 - P_out)**3 #probability that 3 crosslinks form
    P4 = (1 - P_out)**4 #probability that 4 corsslinks form
    ve = conc_cross * ((3/2) * P3 + 2 * P4) #elastically active strands
    g0 = ve #for affine network g0 = rT(ve) we are doing normalized modulus so rt get calculated out

    g0[g0 < 0] = np.nan #only give safe values

    g0_init = g0[0, 0] #initial modulus
    normalized_g0 = g0 / g0_init #normalized modulus 
    real_g0 = normalized_g0 * g0_init_user #corrected for initial modulus

    return comp_concs, Kac_vals, real_g0, CCOMP, KAC #return values


def affine_2d_line(conc_cross, Kab, Kac, comp_conc_range, g0_init_user): #2d version of above
    comp_concs = np.linspace(comp_conc_range[0], comp_conc_range[1], 100)
    Kab_app = Kab / (1 + Kac * (comp_concs / 1000))
    inner = 1 / (2 * conc_cross / 1000 * Kab_app)
    p = (1 + inner) - np.sqrt((1 + inner)**2 - 1)
    P_out = np.sqrt(1 / p - 0.75) - 0.5
    P3 = 4 * P_out * (1 - P_out)**3
    P4 = (1 - P_out)**4
    ve = conc_cross * ((3/2) * P3 + 2 * P4)
    g0 = ve

    g0[g0 < 0] = np.nan

    g0_init = g0[0]
    normalized_g0 = g0 / g0_init
    real_g0 = normalized_g0 * g0_init_user

    return comp_concs, normalized_g0, real_g0

def predict_affine_modulus(K_val, conc_cross, other_K, init_g, comp_concs, fit_target): #for fitting Kac or Kab
    if fit_target == "Kac":
        Kac = K_val[0]
        Kab = other_K
    else:  # fitting Kab
        Kac = other_K
        Kab = K_val[0]

    Kab_app = Kab / (1 + Kac * comp_concs)
    p = (1 + (1 / (2 * conc_cross * Kab_app))) - np.sqrt((1 + (1 / (2 * conc_cross * Kab_app)))**2 - 1)
    P_out = np.sqrt(1 / p - 0.75) - 0.5
    P3 = 4 * P_out * (1 - P_out)**3
    P4 = (1 - P_out)**4
    ve = conc_cross * ((3/2) * P3 + 2 * P4)
    g0 = ve
    g0_int = ve[0]
    normalized_g0 = g0 / g0_int
    return init_g * normalized_g0


# ========================================
# Tau functions
# ========================================

def comp_inhib_p(conc_cross, KabMax, KacMin, c):
    Kab_app = KabMax / (1 + KacMin * c)
    p = (1 + (1 / (2 * conc_cross * Kab_app))) - np.sqrt((1 + (1 / (2 * conc_cross * Kab_app)))**2 - 1)
    return p

def compute_ve(c, conc_cross, KabMax, KacMin):
    p = comp_inhib_p(conc_cross, KabMax, KacMin, c)
    P_out = np.sqrt(1 / p - 3/4) - 0.5
    P3 = 4 * P_out * (1 - P_out)**3
    P4 = (1 - P_out)**4
    ve = conc_cross * ((3/2) * P3 + 2 * P4)
    return ve

def tau_model_physical(c, tau_0, tau_min, conc_cross, KabMax, KacMin):
    ve = compute_ve(c, conc_cross, KabMax, KacMin)
    ve_0 = compute_ve(np.array([0]), conc_cross, KabMax, KacMin)[0]
    ve_scaled = ve / ve_0
    return tau_min + (tau_0 - tau_min) * ve_scaled

# --- Fit & Plot helper ---
def fit_and_plot_physical(conc, tau, conc_cross, KabMax, KacMin):
    tau_0_fixed = tau[0]

    def model(c, tau_min):
        return tau_model_physical(c, tau_0_fixed, tau_min, conc_cross, KabMax, KacMin)

    popt, _ = curve_fit(model, conc, tau, p0=[0.1], bounds=([0.0001], [100.0]))

    # Predictions
    c_pred = np.linspace(0, conc.max(), 300)
    tau_pred = model(c_pred, *popt)

    # R²
    residuals = tau - model(conc, *popt)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((tau - np.mean(tau))**2)
    r_squared = 1 - (ss_res / ss_tot)

    return c_pred, tau_pred, popt[0], r_squared


# ========================================
# Helper function for choosing model
# ========================================
def compute_modulus_surface(conc_cross, Kab, Kac_range, comp_conc_range, g0_init_user, model_choice): #for picking 3d graph
    if model_choice == "Phantom":
        return phantom_network_modulus_surface(conc_cross, Kab, Kac_range, comp_conc_range, g0_init_user)
    else:
        return affine_network_modulus_surface(conc_cross, Kab, Kac_range, comp_conc_range, g0_init_user)
    
def compute_modulus_line(conc_cross, Kab, Kac, comp_conc_range, g0_init_user, model_choice): #for picking 2d graph
    if model_choice == "Phantom":
        return phantom_2d_line(conc_cross, Kab, Kac, comp_conc_range, g0_init_user)
    else:
        return affine_2d_line(conc_cross, Kab, Kac, comp_conc_range, g0_init_user)

def loss(K_val, conc_cross, other_K, init_g, comp_concs, exp_moduli, fit_target,model_choice): #for picking fitting model
    if model_choice == "Phantom":
        pred = predict_phantom_modulus(K_val, conc_cross, other_K, init_g, comp_concs, fit_target)
    else:
        pred = predict_affine_modulus(K_val, conc_cross, other_K, init_g, comp_concs, fit_target)
    return np.sum((pred - exp_moduli)**2)

# ==========================================================
# App 1: 3D Competitive Inhibition Modulus Prediction Tool
# ==========================================================
def app1():
    st.title("3D Competitive Inhibition Modulus Prediction Tool")

    # Sidebar model parameters
    st.sidebar.header("Model Parameters")
    Kab = st.sidebar.number_input("$K_{a,XL}$ (Keq of crosslink)", min_value=0.0, value=2185.0)
    conc_cross = st.sidebar.number_input("Crosslink Concentration (mM)", min_value=0.0, value=80.0)
    g0_init_user = st.sidebar.number_input("Initial Modulus (kPa)", min_value=0.0, value=20.0)
    model_choice = st.radio("Choose Network Model:", ["Phantom", "Affine"])

    # Competitor parameter ranges
    st.sidebar.markdown("### Competitor Parameters")
    kac_min = st.sidebar.number_input("$K_{a,C}$ min", min_value=0.0, value=0.0)
    kac_max = st.sidebar.number_input("$K_{a,C}$ max", min_value=0.0, value=2000.0)
    comp_min = st.sidebar.number_input("Competitor conc min (mM)", min_value=0.0, value=0.0)
    comp_max = st.sidebar.number_input("Competitor conc max (mM)", min_value=0.0, value=100.0)

    # Highlight options
    st.sidebar.markdown("### Highlight Options")
    highlight_mode = st.sidebar.radio("Highlight By:", ["Modulus", "Kac"])

    # Compute surface
    comp_concs, kac_vals, modulus, CCOMP, KAC = compute_modulus_surface(
        conc_cross, Kab, (kac_min, kac_max), (comp_min, comp_max), g0_init_user, model_choice
    )

    # Apply highlighting logic
    if highlight_mode == "Modulus":
        target_modulus = st.sidebar.number_input("Target Modulus (kPa)", min_value=0.0, value=7.0)
        tolerance = st.sidebar.number_input("Tolerance (±)", min_value=0.01, max_value=5.0, value=0.1)

        mask = np.abs(modulus - target_modulus) <= tolerance
        highlight_x = CCOMP[mask]
        highlight_y = KAC[mask]
        highlight_z = modulus[mask]
        highlight_name = f'Modulus ≈ {target_modulus}±{tolerance}'
        highlight_color = "red"

    else:  # Highlight by KaC
        target_kac = st.sidebar.number_input("Target $K_{a,C}$", min_value=0.0, value=1000.0)
        tolerance = st.sidebar.number_input("Tolerance (±)", min_value=0.01, max_value=500.0, value=50.0)

        mask = np.abs(KAC - target_kac) <= tolerance
        highlight_x = CCOMP[mask]
        highlight_y = KAC[mask]
        highlight_z = modulus[mask]
        highlight_name = f'Ka,C ≈ {target_kac}±{tolerance}'
        highlight_color = "red"

    # Build figure
    fig = go.Figure(data=[
        go.Surface(
            z=modulus,
            x=CCOMP,
            y=KAC,
            colorscale='Blues',
            colorbar=dict(title='Modulus (kPa)'),
            hovertemplate="Conc: %{x:.2f} mM<br>Kac: %{y:.2f} M^-1<br>Modulus: %{z:.2f} kPa<extra></extra>"
        ),
        go.Scatter3d(
            x=highlight_x,
            y=highlight_y,
            z=highlight_z,
            mode='markers',
            marker=dict(size=4, color=highlight_color),
            name=highlight_name,
        )
    ])

    fig.update_layout(
        scene=dict(
            xaxis_title='Competitor Conc (mM)',
            yaxis_title='Competitor Keq (Ka,C)',
            zaxis_title='Predicted Modulus (kPa)',
            camera=dict(eye=dict(x=2, y=2, z=1.5))
        ),
        height=700,
        title="Modulus Surface with Competitive Inhibition"
    )
    st.plotly_chart(fig, use_container_width=True)


# ==========================================================
# App 2: Modulus Competitive Inhibition Tool
# ==========================================================
def app2():
    st.title("Modulus Prediction")
    model_choice = st.radio("Choose Network Model:", ["Phantom", "Affine"])
    name = st.text_input("Input File Name", value='phantom_pred_modulus_default')
    conc_cross_mM = st.number_input("Crosslink Concentration (mM)", value=80.0)
    Kab = st.number_input("Crosslink Ka (M^-1)", value=2185.0)
    Kac = st.number_input("Competitor Ka (M^-1)", value=280.0)

    unit_options = {"Pa":1, "kPa":1e3, "MPa":1e6, "GPa":1e9}
    unit_choice = st.selectbox("Units for Initial Gel Stiffness", list(unit_options.keys()), index=1)
    init_g_input = st.number_input(f"Initial Gel Stiffness ({unit_choice})", value=21.85)
    init_g = init_g_input * unit_options[unit_choice]

    cmin_mM = st.number_input("Minimum Competitor Concentration (mM)", value=0.0)
    cmax_mM = st.number_input("Maximum Competitor Concentration (mM)", value=80.0)

    crange_mM = np.linspace(cmin_mM, cmax_mM, 100)*100

    comp_concs, normalized_g0, pred_g0 = compute_modulus_line(conc_cross_mM, Kab, Kac, crange_mM, init_g, model_choice)

##    fig_norm = go.Figure()
##    fig_norm.add_trace(go.Scatter(x=comp_concs, y=normalized_g0, mode='lines'))
##    fig_norm.update_layout(xaxis_title="Competitor Conc (mM)", yaxis_title="Normalized Modulus")
##    st.plotly_chart(fig_norm, use_container_width=True)

    fig_pred = go.Figure()
    fig_pred.update_layout(
        width=600,  # corresponds to figsize=6 inches roughly
        height=600,
        font=dict(size=25)
        )
    fig_pred.add_trace(go.Scatter(x=comp_concs, y=pred_g0 / unit_options[unit_choice]))
    fig_pred.update_layout(xaxis_title="Competitor Conc (mM)", yaxis_title=f"Predicted Modulus ({unit_choice})")

    st.plotly_chart(fig_pred, use_container_width=True)

    df = pd.DataFrame({
        'Competitor Concentration (mM)': comp_concs,
        'Normalized g₀': normalized_g0,
        f'Predicted g₀ ({unit_choice})': pred_g0 / unit_options[unit_choice]
    })
    st.download_button("Download CSV", df.to_csv(index=False), name+".csv")


# ==========================================================
# App 3: Estimate Keq from Modulus
# ==========================================================
def app3():
    st.title("Estimate Keq from Modulus")
    model_choice = st.radio("Choose Network Model:", ["Phantom", "Affine"])

    fit_target = st.radio("Which parameter to estimate?", ["Kac", "Kab"])
    init_g = st.number_input("Initial modulus (G₀)", value=21.85)
    conc_cross_mM = st.number_input("Crosslink concentration (mM)", value=80.0)

    if fit_target == "Kac":
        fixed_Kab = st.number_input("Kab (Crosslink Keq, M^-1)", value=2185.0)
        guess, fixed_param = 1000, fixed_Kab
    else:
        fixed_Kac = st.number_input("Kac (Competitor Keq, M^-1)", value=450.0)
        guess, fixed_param = 2185, fixed_Kac

    comp_concs_str = st.text_input("Competitor concs (mM, comma-separated)", "0,10,20,30,40")
    exp_moduli_str = st.text_input("Observed moduli", "21.85,17.05,14.98,12.54,11.60")

    try:
        comp_concs_mM = np.array([float(v.strip()) for v in comp_concs_str.split(',')])
        exp_moduli = np.array([float(v.strip()) for v in exp_moduli_str.split(',')])
        comp_concs = comp_concs_mM / 1000
        conc_cross = conc_cross_mM / 1000

        result = minimize(
            loss,
            x0=[guess],
            args=(conc_cross, fixed_param, init_g, comp_concs, exp_moduli, fit_target,model_choice),
            bounds=[(1, 1e25)]
        )
        fit_value = result.x[0]
        if model_choice == "Phantom":
            fit_vals = predict_phantom_modulus([fit_value], conc_cross, fixed_param, init_g, comp_concs, fit_target)
        else:
            fit_vals = predict_affine_modulus([fit_value], conc_cross, fixed_param, init_g, comp_concs, fit_target)
        

        r2 = r2_score(exp_moduli, fit_vals)

        st.success(f"Estimated {fit_target}: {fit_value:.2e} M^-1 (R²={r2:.3f})")

        smooth_concs_mM = np.linspace(0, max(comp_concs_mM)*1.1, 200)
        smooth_concs = smooth_concs_mM / 1000
        if model_choice == "Phantom":
            smooth_fit = predict_phantom_modulus([fit_value], conc_cross, fixed_param, init_g, smooth_concs, fit_target)
        else:
            smooth_fit = predict_affine_modulus([fit_value], conc_cross, fixed_param, init_g, smooth_concs, fit_target)
        

        fig, ax = plt.subplots()
        ax.plot(comp_concs_mM, exp_moduli, 'o', label='Experimental', color='black')
        ax.plot(smooth_concs_mM, smooth_fit, '-', label='Model Fit', color='blue')
        ax.set_xlabel("Competitor Conc (mM)",fontsize=18)
        ax.set_ylabel("Modulus",fontsize=18)
        ax.legend()
        st.pyplot(fig)
    except Exception as e:
        st.error(f"Error: {e}")

# ==========================================================
# App 4: Tau Prediction & Fitting
# ==========================================================
def app4():
    st.title("Relaxation Time (τ) Prediction")

    # --- User Inputs ---
    st.sidebar.header("Model Parameters")
    conc_cross = st.sidebar.number_input("Crosslink concentration (mM)", value=80.0) / 1000  # M
    Ka_XL = st.sidebar.number_input("Ka,XL (Crosslink Ka, M^-1)", value=2185.0)
    Kac = st.sidebar.number_input("Ka,C (Competitor Ka, M^-1)", value=489.0)

    st.sidebar.subheader("Experimental Data")
    comp_concs_str = st.text_input("Competitor concs (mM, comma-separated)", "0,10,20,30,40")
    exp_tau_str = st.text_input("Observed τ values (s, comma-separated)", 
                                "4.56,3.86,3.47,3.20,2.87")
    
    # --- Parse data ---
    try:
        conc = np.array([float(x) for x in comp_concs_str.split(",")]) / 1000  # mM → M
        tau = np.array([float(y) for y in exp_tau_str.replace(" ", ",").split(",")])
    except Exception as e:
        st.error(f"Error parsing input: {e}")
        return

    # --- Fit & plot ---
    c_pred, tau_pred, tau_min_fit, r2 = fit_and_plot_physical(conc, tau, conc_cross, Ka_XL, Kac)

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.errorbar(conc*1e3, tau, fmt="o", label="Experimental")
    ax.plot(c_pred*1e3, tau_pred, "-", label=f"Fit (τ_min = {tau_min_fit:.2f} s, R²={r2:.3f})")
    ax.set_xlabel("Competitor Concentration (mM)")
    ax.set_ylabel("Relaxation Time τ (s)")
    ax.legend()
    st.pyplot(fig)

# ==========================================================
# App 5: Visualization of Inhibited Network
# ==========================================================

def app5():
    st.title("Junction Grid Visualization")

    st.markdown("""
    This tool generates a grid of junction images with probabilities 
    determined by the crosslinking model. 
    """)

    # --- user inputs ---
    conc_cross = st.number_input("Crosslink concentration (mM)", value=80.0)
    Kab = st.number_input("Crosslink Ka (M⁻¹)", value=2185.0)
    Kac = st.number_input("Competitor Ka (M⁻¹)", value=500.0)
    comp_conc = st.number_input("Competitor concentration (mM)", value=10.0)
    rows = st.number_input("Grid rows", value=5, step=1)
    cols = st.number_input("Grid cols", value=5, step=1)
    

    # --- select image set depending on conditions ---
    if comp_conc == 0:
        image_paths = [
            "pics/0arms.png",
            "pics/1arm.png",
            ["pics/2arms_1.png", "pics/2arms_2.png"], 
            "pics/3arms.png",
            "pics/4arms.png"
        ]
    elif Kac == 0:
        image_paths = [
            "pics/0arms.png",
            "pics/1arm.png",
            ["pics/2arms_1.png", "pics/2arms_2.png"], 
            "pics/3arms.png",
            "pics/4arms.png"
        ]
    else:
        image_paths = [
            "pics/comp0arms.png",
            "pics/comp1arm.png",
            ["pics/comp2arms_1.png", "pics/comp2arms_2.png"], 
            "pics/comp3arms.png",
            "pics/4arms.png"
        ]

    if st.button("Generate Grid"):
        fig = junction_grid_border_overlay_streamlit(
            conc_cross, Kab, Kac, comp_conc, image_paths,
            grid_size=(rows, cols),
            border_thickness=10,
        )
        st.pyplot(fig)


def junction_grid_border_overlay_streamlit(conc_cross, Kab, Kac, comp_conc, image_paths, 
                                           grid_size, border_thickness=20, dpi=600):


    # --- probability calculations ---
    Kab_app = Kab / (1 + Kac * (comp_conc / 1000))
    p = (1 + 1 / (2 * conc_cross / 1000 * Kab_app)) - np.sqrt(
        (1 + 1 / (2 * conc_cross / 1000 * Kab_app))**2 - 1
    )
    P_out = np.sqrt(1 / p - 0.75) - 0.5

    P0 = P_out**4
    P1 = 4 * (1 - P_out) * P_out**3
    P2 = 6 * (1 - P_out)**2 * P_out**2
    P3 = 4 * P_out * (1 - P_out)**3
    P4 = (1 - P_out)**4
    probs = np.array([P0, P1, P2, P3, P4])
    probs = probs / probs.sum()
                                               
    st.markdown("### Elastically Active Strands")                                              
    st.write(f"Probability of 4 junctions = {P4:.3f}")
    st.write(f"Probability of 3 junctions = {P3:.3f}") 

    st.markdown("### Elastically Inactive Strands")
    st.write(f"Probability of 2 junctions = {P2:.3f}")
    st.write(f"Probability of 1 junctions = {P1:.3f}") 
    st.write(f"Probability of 0 junctions = {P0:.3f}")
  
    rows, cols = grid_size
    choices = np.random.choice(len(probs), size=(rows, cols), p=probs)

    fig, ax = plt.subplots(figsize=(cols, rows), dpi=dpi)
    ax.set_xlim(0, cols)
    ax.set_ylim(0, rows)
    ax.axis('off')

    for i in range(rows):
        for j in range(cols):
            idx = choices[i, j]

            if isinstance(image_paths[idx], list):
                img_path = random.choice(image_paths[idx])
            else:
                img_path = image_paths[idx]

            img = Image.open(img_path).convert("RGBA")
            draw = ImageDraw.Draw(img)
            width, height = img.size
            draw.rectangle([0, 0, width-3, height-3], outline="black", width=border_thickness)

            angle = random.choice([0, 90, 180, 270])
            img = img.rotate(angle, expand=False)

            ax.imshow(img, extent=[j, j+1, rows-i-1, rows-i])

    plt.close(fig)
                                               
    return fig




# ==========================================================
# Main Controller
# ==========================================================
page = st.sidebar.radio("Select a Page", ["Description", "Modulus Prediction", "Keq Prediction","Tau Prediction & Fitting","3D Surface Explorer","Network Visualization"])
if page == "Description":
    description_page()
elif page == "3D Surface Explorer":
    app1()
elif page == "Modulus Prediction":
    app2()
elif page == "Keq Prediction":
    app3()
elif page == "Tau Prediction & Fitting":
    app4()
else:
    app5()
