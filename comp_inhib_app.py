import streamlit as st
import numpy as np
import pandas as pd
import plotly.graph_objects as go

st.title("Gelation Competitive Inhibition Tool")
st.markdown(
    r"""
### What does this tool do?
This app predicts how **competitive inhibition** affects gel stiffness (*g₀*) in an ideal dynamic polymer network.

- **Crosslink reaction:** $A + B \rightleftharpoons AB$ with equilibrium constant $K_{ab}$
- **Competitor reaction:** $AB + C \rightleftharpoons AC + B$ with equilibrium constant $K_{ac}$
- **Goal:** See how adding competitor $C$ lowers the effective crosslink density and modulus.

Use the controls below to explore parameter space and download the calculated data for further analysis.
"""
)

with st.expander("Technical details & equations (click to expand)"):
    st.markdown(r"""
The apparent equilibrium constant in the presence of competitor based on enzyme competitive inhibition is

$$K_{ab}^{\text{app}} = \frac{K_{ab}}{1 + K_{ac}[C]}$$

The conversion of chains being cross‑linked in an ideal network is

$$p = 1 + \frac{1}{2N_aK_{ab}^{\text{app}}} - \sqrt{\left(1+\frac{1}{2N_aK_{ab}^{\text{app}}}\right)^2-1}$$

where $N_a$ is the total crosslink concentration.

The initial modulus is approximated by using ideal network theory as

$$g_0 = \frac{N_a}{16}\bigl(3-\sqrt{\tfrac{4}{p}-3}\bigr)^3\bigl(\sqrt{\tfrac{4}{p}-3}+1\bigr).$$
""")

# --- Inputs (all concentrations in mM) ---
conc_cross_mM = st.number_input("Crosslink Concentration (mM)", value=80.0)
Kab = st.number_input("Crosslink Keq (M^-1)", value=2185.0)
Kac = st.number_input("Competitor Keq (M^-1)", value=280.0)

# Units
unit_options = {"Pa":1, "kPa":1e3, "MPa":1e6, "GPa":1e9}
unit_choice = st.selectbox("Units for Initial Gel Stiffness", list(unit_options.keys()), index=1)
init_g_input = st.number_input(f"Initial Gel Stiffness ({unit_choice})", value=21.85)
init_g = init_g_input * unit_options[unit_choice]  # Convert to Pa

cmin_mM = st.number_input("Minimum Competitor Concentration (mM)", value=0.0)
cmax_mM = st.number_input("Maximum Competitor Concentration (mM)", value=80.0)

# --- Unit conversion: mM → M ---
conc_cross = conc_cross_mM / 1000  # M
crange_mM = np.linspace(cmin_mM, cmax_mM, 1000)
crange = crange_mM / 1000  # M

# --- Model Calculations ---
Kab_app = Kab / (1 + Kac * crange)
inner_term = (1 + (1 / (2 * conc_cross * Kab_app)))**2 - 1
inner_term[inner_term < 0] = np.nan
p = (1 + (1 / (2 * conc_cross * Kab_app))) - np.sqrt(inner_term)

# g₀ calculation
g0 = np.full_like(p, np.nan)
safe_indices = ((4 / p) - 3) > 0
temp = np.sqrt((4 / p[safe_indices]) - 3)
g0[safe_indices] = (conc_cross / 16) * (3 - temp)**3 * (temp + 1)

# Handle possible nan at g0[0]
if np.isnan(g0[0]):
    st.error("Initial modulus could not be calculated. Please adjust input parameters.")
    st.stop()

# Normalize and scale
normalized_g0 = g0 / g0[0]
pred_g0 = init_g * normalized_g0

# --- Interactive Plotly Graphs ---

# Normalized Modulus
st.subheader("Normalized Modulus vs. Competitor Concentration (mM)")
fig_norm = go.Figure()
fig_norm.add_trace(go.Scatter(
    x=crange_mM,
    y=normalized_g0,
    mode='lines',
    name='Normalized g₀',
    hovertemplate='[C] = %{x:.2f} mM<br>Normalized g₀ = %{y:.3f}<extra></extra>'
))
fig_norm.update_layout(
    xaxis_title="Competitor Concentration (mM)",
    yaxis_title="Normalized Modulus",
    title="Normalized Modulus vs. Competitor Concentration",
    hovermode="x",
    template="simple_white"
)
st.plotly_chart(fig_norm, use_container_width=True)

# Predicted Modulus
st.subheader(f"Predicted Modulus vs. Competitor Concentration (mM) [{unit_choice}]")
fig_pred = go.Figure()
fig_pred.add_trace(go.Scatter(
    x=crange_mM,
    y=pred_g0 / unit_options[unit_choice],
    mode='lines',
    name=f'Predicted g₀ ({unit_choice})',
    hovertemplate='[C] = %{x:.2f} mM<br>g₀ = %{y:.3f} ' + unit_choice + '<extra></extra>'
))
fig_pred.update_layout(
    xaxis_title="Competitor Concentration (mM)",
    yaxis_title=f"Predicted Modulus ({unit_choice})",
    title="Predicted Modulus vs. Competitor Concentration",
    hovermode="x",
    template="simple_white"
)
st.plotly_chart(fig_pred, use_container_width=True)

# --- Downloadable Data ---
df = pd.DataFrame({
    'Competitor Concentration (mM)': crange_mM,
    'Normalized g₀': normalized_g0,
    f'Predicted g₀ ({unit_choice})': pred_g0 / unit_options[unit_choice]
})

st.subheader("Download Output Data")
st.download_button("Download CSV", df.to_csv(index=False), "predicted_g0_mM.csv")
