![Python](https://img.shields.io/badge/python-3.9%2B-blue)
![Streamlit](https://img.shields.io/badge/Streamlit-1.x-red)
![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)

# InhibNet: A Tool for Predicting Hydrogel Mechanics from Competitive Inhibition of Dynamic Crosslinks.

Welcome! This app allows you to explore how competitive inhibition affects the mechanical properties of polymer networks.
Associated publication: TBD
Preprint: https://doi.org/10.26434/chemrxiv-2025-97wbj

For the cloud version of the app please visit: https://inhibnet.streamlit.app/


# Background on Work

This app allows you to explore how competitive inhibition affects the mechanical properties of polymer networks.

---

The purpose of this work is to quantify how a polymer network's material properties can be altered by adding a small molecule that can compete with a crosslink. It has long been shown in the literature that the material properties (i.e., modulus and relaxation time) of dynamic polymer networks are determined by the dynamics of the crosslink.  

We propose that adding a small molecule competitor can alter the dynamics of the crosslink through the strength of the competing interaction. In turn, the altered dynamics lead to a predictable change in material properties.  
This framework has been validated across boronate ester and hydrazone dynamic crosslinked networks, demonstrating predictive accuracy for both shear modulus and relaxation time.

We hypothesized that principles from enzyme competitive inhibition could be adapted to dynamic hydrogels to provide a similarly simple framework for predicting how key network properties change in the presence of competing species.  

In this analogy, the apparent equilibrium constant ($K_{a,app}$), commonly used in enzyme kinetics to describe effective binding affinities under competition, can be translated to dynamic networks as a predictor of mechanical response.  Due to the crosslinks existing in a ternary equilibrium among unbound strands, crosslinked pairs, and crosslink–competitor complexes, we reasoned that the association constants for crosslink formation ($K_{a,XL}$) and competitor binding ($K_{a,C}$) could be used within the $K_{a,app}$ framework to capture the effective network behavior.

**Figure:** Ternary equilibrium, allowing for the $K_{a,app}$ assumption  
![Ternary Equilibrium](pics/TernaryEq.jpg)

---

## Equations

### 1. Crosslink Conversion

$$
p = \left(1 + \frac{1}{2 N_a K_{a,app}}\right) - \sqrt{\left(1 + \frac{1}{2 N_a K_{a,app}}\right)^2 - 1}
$$

where:

- $p$ — crosslink conversion  
- $N_a$ — concentration of crosslinks in solution  

---

### 2. Shear Modulus (Affine Network)

Define the effective probability term:

$$
P_{out} = \sqrt{\frac{1}{p} - \frac{3}{4}} - \frac{1}{2}
$$

Then the probability that 3 crosslinks will form ($P_3$):

$$
P_3 = 4 P_{out} (1 - P_{out})^3
$$

Then the probability that 4 crosslinks will form ($P_4$):

$$
P_4 = (1 - P_{out})^4
$$

The elastically active chains ($v_e$) are:

$$
v_e = \frac{N_a}{4} \left(\frac{3}{2}P_3 + 2P_4 \right)
$$

The modulus is given by:

$$
\frac{G_p}{k_b T} = v_e
$$

---

### 3. Shear Modulus (Phantom Network)

For more dilute systems (close to the overlap concentration), it is more appropriate to use a phantom model for network elasticity, which accounts for the movement of strand junctions:

$$
\frac{G_p}{k_b T} = v_e - \mu_e
$$

where $\mu_e$ represents the number density of elastically active junctions:

$$
\mu_e = \frac{N_a}{4} (P_3 + P_4)
$$

Simplified modulus equation:

$$
G_p = \frac{N_a}{16} \left(3 - \sqrt{\frac{4}{p} - 3}\right)^3 \left(\sqrt{\frac{4}{p} - 3} + 1\right)
$$

---

### 4. Relaxation Time Prediction

The relaxation time ($\tau$) decreases with increasing competitor concentration due to disruption of elastically active crosslinks.  
We model this using a Langmuir-type decay relation scaled by the fraction of active strands ($v_e$):

$$
\tau(C) = \tau_0 - (\tau_0 - \tau_{min}) \left(\frac{2 v_e}{N_a}\right)
$$

where:  
- $\tau_0$ — relaxation time in the absence of competitor  
- $\tau_{min}$ — minimum $\tau$ value as competitor concentration approaches infinity  
- $v_e$ and $N_a$ — as defined above
## Tools Available in App

**Five main tools are included:**

1. **Modulus vs. Competitor Concentration**  
   Generate a 2D plot showing how the shear modulus changes with increasing competitor concentration.  
   Input features of the crosslink, competitor, and the initial stiffness of your gel.  
   The app outputs predicted modulus values across a range of competitor concentrations and allows you to download a CSV of the results.

2. **Predict $K_{a,C}$ or $K_{a,XL}$ from Experimental Data**  
   Input experimental modulus and concentration data to estimate either the competitor association constant ($K_{a,C}$) or the crosslink association constant ($K_{a,XL}$).  
   This is useful for fitting experimental data to the inhibition framework.

3. **Relaxation Time as a Function of Competitor**  
   Predict how the relaxation time ($\tau$) changes with increasing competitor concentration.  
   Simply input the uninhibited relaxation time ($\tau_0$) and thermodynamic information about the crosslink and competitor, and the app will generate a predicted relaxation curve.

4. **Modulus Surface Visualization**  
   Explore the full design space of the model in 3D.  
   Input features such as the crosslink association constant ($K_{a,XL}$), crosslink concentration ($N_a$), competitor association constant ($K_{a,C}$), and competitor concentration ([C]).  
   The app generates an interactive 3D plot of modulus vs. $K_{a,C}$ vs. [C], allowing you to visualize how competitor strength and concentration jointly affect the network stiffness.

5. **Network Visualization**  
   Visualize how elastically active strands vary under different parameters.  
   This tool simulates a 2D grid of network states based on the mean-field model, calculating the probability of each node being elastically active under given thermodynamic conditions.

---

### Network Models

You can choose between two network elasticity models:

- **Phantom Network** – accounts for network fluctuations and is most appropriate for dilute polymer networks (near the overlap concentration).  
- **Affine Network** – assumes a fully connected network with affine deformations, best for more concentrated systems (typically >4× the overlap concentration).

---

### Instructions for App Use

- Use the **sidebar** to navigate between tools.  
- Adjust model parameters interactively on each page.  
- For the **Modulus Surface Visualization** tool, hover over the 3D plot to view modulus values at each coordinate.  
- Use the **model selector** to toggle between Phantom and Affine network calculations.


       
## Detailed Local Install Instructions

### 1) Install Python
- Install **Python 3.9+** from https://www.python.org/downloads/  
- On Windows, check **“Add Python to PATH”** during installation.
### 2) Download this repository and navigate to the path
- This repo can be downloaded from the green code button
   - Download zip
   - Then you can expand the zip file
   - It should be labeled as a folder "InhibNet-Main"
- Now we need to navigate to the folder "InhibNet-Main"
   - This can be done through Terminal, command line, etc.
- Navigate to the folder containing InhibNet:

```bash
cd /path/to/InhibNet-main
#note that this won't be your exact path, find where you saved the Inhib-Net folder in file explorer, or finder and you can copy the file path
```

### 3) Create a virtual environment named "InhibNetEnv"
- you can call this environment whatever you want here I call it InhibNetEnv
- note: some people may need to write "py" instead of python because of how your computer saved python
```bash
python -m venv InhibNetEnv
```

### 4) Activate your environment

- Windows (PowerShell)
```bash
InhibNetEnv\Scripts\Activate.ps1
```
- Windows (cmd)
```bash
InhibNetEnv\Scripts\activate.bat
```
- macOS / Linux
```bash
source InhibNetEnv/bin/activate 
```
- After activating your environment instead of something like (base) you will see (InhibNetEnv) at the start of each line of code
### 5) Install Dependencies
```
# Make sure your pip is up to date
python -m pip install --upgrade pip
```
```
# use the requirements file in the repository to make sure you have the required packages to run the app
pip install -r requirements.txt
```
### 6) Run the app
```
streamlit run INHIBNET_App.py
```
### 7) Running the app locally in the future
- Now that you have all the requirements ready, you can run the app with much fewer steps
- You simply need to activate your enviroment (step 4)
- Then run the app (step 6)
