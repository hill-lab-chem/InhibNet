import streamlit as st
import numpy as np
import plotly.graph_objects as go

# Compute the modulus surface
def compute_modulus_surface(conc_cross, Kab, Kac_range, comp_conc_range, g0_init_user):
    Kac_vals = np.linspace(Kac_range[0], Kac_range[1], 100)
    comp_concs = np.linspace(comp_conc_range[0], comp_conc_range[1], 100)
    CCOMP, KAC = np.meshgrid(comp_concs, Kac_vals)

    Kab_app = Kab / (1 + KAC * (CCOMP / 1000))
    inner = 1 / (2 * conc_cross / 1000 * Kab_app)
    p = (1 + inner) - np.sqrt((1 + inner) ** 2 - 1)

    sqrt_term = (4 / p) - 3
    sqrt_term[sqrt_term < 0] = np.nan
    root = np.sqrt(sqrt_term)
    g0 = (conc_cross / 1000 / 16) * (3 - root) ** 3 * (root + 1)
    g0[g0 < 0] = np.nan

    g0_init = g0[0, 0]
    normalized_g0 = g0 / g0_init
    real_g0 = normalized_g0 * g0_init_user
    return comp_concs, Kac_vals, real_g0, CCOMP, KAC

# Streamlit UI
st.title("3D Competitive Inhibition Modulus Prediction Tool")

st.sidebar.header("Model Parameters")
Kab = st.sidebar.number_input("Kab (Keq of crosslink)", min_value=0.0, value=2000.0)
conc_cross = st.sidebar.number_input("Crosslink Concentration (mM)", min_value=0.0, value=80.0)
g0_init_user = st.sidebar.number_input("Initial Modulus (kPa)", min_value=0.0, value=20.0)

st.sidebar.markdown("### Competitor Parameters")
kac_min = st.sidebar.number_input("Kac min", min_value=0.0, value=0.0)
kac_max = st.sidebar.number_input("Kac max", min_value=0.0, value=2000.0)

comp_min = st.sidebar.number_input("Competitor conc min (mM)", min_value=0.0, value=0.0)
comp_max = st.sidebar.number_input("Competitor conc max (mM)", min_value=0.0, value=100.0)

st.sidebar.markdown("### Target Modulus Highlight")
target_modulus = st.sidebar.number_input("Target Modulus (kPa)", min_value=0.0, value=7.0)
tolerance = st.sidebar.slider("Tolerance (±)", min_value=0.01, max_value=5.0, value=0.5)

# Compute surface
comp_concs, kac_vals, modulus, CCOMP, KAC = compute_modulus_surface(
    conc_cross, Kab, (kac_min, kac_max), (comp_min, comp_max), g0_init_user
)

# Mask for values near the target modulus
mask = np.abs(modulus - target_modulus) <= tolerance
highlight_x = CCOMP[mask]
highlight_y = KAC[mask]
highlight_z = modulus[mask]

# Set up camera (zoomed out)
camera = dict(eye=dict(x=2, y=2, z=1.5))

# Create Plotly figure
fig = go.Figure(data=[
    go.Surface(
        z=modulus,
        x=comp_concs,
        y=kac_vals,
        colorscale='Blues',
        colorbar=dict(title='Modulus (kPa)'),
        name='Modulus Surface',
        hovertemplate="Conc: %{x:.2f} mM<br>Kac: %{y:.2f} M^-1<br>Modulus: %{z:.2f} kPa<extra></extra>"
    ),
    go.Scatter3d(
        x=highlight_x,
        y=highlight_y,
        z=highlight_z,
        mode='markers',
        marker=dict(size=4, color='red'),
        name=f'Modulus ≈ {target_modulus}±{tolerance}',
        hovertemplate="Conc: %{x:.2f} mM<br>Kac: %{y:.2f} M^-1<br>Modulus: %{z:.2f} kPa<extra></extra>"
    )
])


fig.update_layout(
    scene=dict(
        xaxis_title='Competitor Conc (mM)',
        yaxis_title='Competitor Keq (Kac)',
        zaxis_title='Predicted Modulus (kPa)',
        camera=camera
    ),
    margin=dict(l=0, r=0, b=0, t=40),
    height=700,
    title="Modulus Surface with Competitive Inhibition"
)

# Display the interactive chart
st.plotly_chart(fig, use_container_width=True)

# Additional info
st.markdown("### Instructions")
st.markdown("""
- Adjust the parameters in the sidebar to explore how competitive inhibition affects modulus.
- Use the "Target Modulus" input to highlight regions near that modulus value.
- Hover over the 3D surface to see precise values.
""")
