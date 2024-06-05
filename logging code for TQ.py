import matplotlib.pyplot as plt
import CoolProp.CoolProp as cp
import numpy as np 
import statistics as st
from tabulate import tabulate
# Defining a function to produce a TQ plot
def plot_TQ_diagram(ax, P2, P3, H2, H3, T_water_inlet, T_water_outlet, pinch_diff, fluid, m_dot_water, m_dot_refrigerant):
    # Constants for water
    c_p_water = 4.18  # kJ/kg.K, specific heat capacity of water

    # Discretize the heat exchanger process
    num_points = 100
    Q_refrigerant = np.linspace(0, ((H2 - H3) * m_dot_refrigerant / 1000), num_points)  # Heat transferred by refrigerant in kJ
    Q_water = np.linspace(0, m_dot_water * c_p_water * (T_water_outlet - T_water_inlet), num_points)  # Heat absorbed by water in kJ

    # Temperature profiles
    Hr_profile = np.linspace(H3, H2, num_points)  # Makes a linear enthalpy profile to calculate temperature profile for refrigerant
    Tr_profile = cp.PropsSI('T', 'H', Hr_profile, 'P', np.linspace(P3, P2, num_points), fluid)  # Temperature profile of the refrigerant
    T_refrigerant = Tr_profile  # T2 decreasing to T3
    T_water = np.linspace(T_water_inlet, T_water_outlet, num_points)  # T_water_in increasing to T_water_out

    # Plotting the T-Q diagram
    ax.plot(Q_refrigerant, T_refrigerant - 273.15, label='Refrigerant')  # Convert K to °C
    ax.plot(Q_water, T_water - 273.15, label='Water')  # Convert K to °C

# Initiate lists for values, number is circuit locations, w is water, a is air, r is refrigerant
# p1 is compressor inlet, 2 is outlet, see diagram in handout for reference
tim = []    
T1w = []
T2w = []
T1a = []
T2a = []
T1r = []
T2r = []
T3r = []
T4r = []
p1=[]
p2 = []
I = []
qw= []

# Pick the text file to use
file_name = '/Users/levendavies/Documents/Heat Pump/Heat-Pump--1/Heat-Pump--5/run3'
# Add values from files to list 

try:
    with open(file_name, "r") as data_file:  # Open file using "with" statement for automatic closure
        title = data_file.readline()  # Read first header line
        props = data_file.readline()  # Read second header line

        # Now read in data line by line
        for line in data_file:
            vals = line.split() 
            try:
                # Try to process the line as numbers
                numbers = [float(num) for num in line.split()]
                # Process the numbers here # Split each line into "words"
                if len(vals) >= 2:  # Ensure at least two values are present
                    tim.append(float(vals[0]))  # Convert values to float and assign to corresponding list 
                    T1w.append(float(vals[1]) + 273.15)  # Convert the second value to floating point and append to T1w[]
                    T2w.append(float(vals[2]) + 273.15)  # Convert temp to K
                    T1a.append(float(vals[3]) + 273.15)
                    T2a.append(float(vals[4]) + 273.15) 
                    T1r.append(float(vals[5]) + 273.15)
                    T2r.append(float(vals[6]) + 273.15)
                    T3r.append(float(vals[7]) + 273.15)
                    T4r.append(float(vals[8]) + 273.15)
                    p1.append(float(vals[9]) * 10**5)  # Convert pressures to Pa
                    p2.append(float(vals[10]) * 10**5)
                    I.append(float(vals[11]))
                    qw.append(float(vals[12]))
                else:
                    print("Skipping line with insufficient data:", line)
            except ValueError: 
                continue 
except FileNotFoundError:
    print("File not found:", file_name)
except Exception as e:
    print("An error occurred:", e)

# Set up lists for values, h is corresponding enthalpy, m is mass flow, p is power 

h1w = [] 
I_no_comp = 3.2740721739130434  # Calculated current when the fan is on full and compressor is not running.
h2w = []
hdifw = []
h1r = []
h2r = []
h3r = []
h3r_dry01 =[]
h3r_no_p_loss =[]
h4r = []
mw  = []
mr = [] 
mr_dry01 =[]
mr_no_p_loss =[]
ma = []
power_w = []
pw2 =[]
power_r =[]
power_r_dry01 =[]
copw = []
copr = []
wc = []
power_draw = []  # Power drawn by system
power_draw_compressor = []  # Power drawn by compressor
power_draw_actual_compressor = []
copr2 = []
copr_dry01 = []
copr_no_p_loss =[]
s4r = []
s3r = []
s2r = []
s1r = []
irr_gen_comp = []
irr_gen_heat_ex =[]
irr_gen_throttle = []
irr_gen_evap = []
p3r = []
p3r_loss = []
h2r_s =[]
eta_c =[]  # Isentropic compressor efficiency
s1 = []
s2 =[]
p4 = []

for i in range(len(T1w)):
    # Calculate enthalpies using CoolProp
    hdifw.append(4200 * (T2w[i] - T1w[i]))
    
    h1r.append(cp.PropsSI('H', 'P|gas', p1[i], 'T', T1r[i], "R134a"))  # Using T and P here gives severely oscillating values of h1 
    s1.append(cp.PropsSI('S', 'P|gas', p1[i], 'T', T1r[i], "R134a"))
    h2r.append(cp.PropsSI('H', 'T', T2r[i], 'P', p2[i], "R134a"))
    h2r_s.append(cp.PropsSI('H', 'P', p2[i], 'S', s1[i], "R134a"))  # T2 strongly correlated with h2 so need to use pressure for isentropic efficiency
    s2.append(cp.PropsSI('S', 'T', T2r[i], 'P', p2[i], "R134a"))
    eta_c.append((h2r_s[i] - h1r[i]) / (h2r[i] - h1r[i]))
    h3r.append(cp.PropsSI('H', 'T', T3r[i], 'Q', 0.0, "R134a"))  # Assuming wet saturated on exit from the condenser
    h3r_dry01.append(cp.PropsSI('H', 'T', T3r[i], 'Q', 0.1, "R134a"))  # Assuming dryness fraction of 0.1 at condenser exit
    h3r_no_p_loss.append(cp.PropsSI('H', 'T', T3r[i], 'P', p2[i], "R134a"))  # Assuming no pressure loss across condenser
    h4r.append(h3r[i])
    
    power_draw.append(I[i] * 245 * 0.98)  # Now only considering compressor work and using power factor

    mw.append(qw[i] / 60)  # Calculate water mass flow (assume rho = 1000)
    mr.append((mw[i] * hdifw[i]) / ((h2r[i] - h3r[i])))  # Assuming adiabatic condenser 
    mr_dry01.append((mw[i] * hdifw[i]) / ((h2r[i] - h3r_dry01[i])))  # Assuming adiabatic condenser 
    power_draw_actual_compressor.append((I[i] - I_no_comp) * 249 * 0.98)
    power_draw_compressor.append((h2r[i] - h1r[i]) * mr[i]) 
    power_w.append(mw[i] * (hdifw[i]))  # Calculate heat output based on enthalpy gained by water
    copw.append(power_w[i] / power_draw[i])  # Calculate COP based on water
    p4.append(cp.PropsSI('P', 'T', T4r[i], 'Q', 0.0, "R134a"))  # Calculate p4, will be in 2-phase region 

    copr2.append((h2r[i] - h3r[i]) / (h2r[i] - h1r[i]))  # Calculating COP purely thermodynamically
    copr_dry01.append((h2r[i] - h3r_dry01[i]) / (h2r[i] - h1r[i])) 
    copr_no_p_loss.append((h2r[i] - h3r_no_p_loss[i]) / (h2r[i] - h1r[i])) 
    # Calculating air mass flow rate 
    ma.append(mr[i] * (h4r[i] - h1r[i]) / (1005 * (T2a[i] - T1a[i])))

evap_p_loss = (st.mean(p4) / 10**5 - st.mean(p1) / 10**5) / (st.mean(p4) / 10**5)
T3_compare = []
# Looking at irreversible entropy generation and pressure losses.
for i in range(len(T1w)):
    s3r.append(cp.PropsSI('S', 'T', T3r[i], 'Q', 0, "R134a"))  
    s2r.append(cp.PropsSI('S', 'T', T2r[i], 'P', p2[i], "R134a"))
    s1r.append(cp.PropsSI('S', 'P|gas', p1[i], 'T', T1r[i], "R134a"))
    s4r.append(cp.PropsSI('S', 'P', p4[i], 'H', h4r[i], "R134a"))
    p3r.append(cp.PropsSI('P', 'Q', 0.0, 'T', T3r[i], "R134a"))
    p3r_loss.append(p2[i] - p3r[i])
    irr_gen_comp.append(s2r[i] - s1r[i])
    irr_gen_throttle.append(s4r[i] - s3r[i])

cond_p_loss = st.mean(p3r_loss) / st.mean(p2)
ma_av = st.mean(ma)
irr_gen_comp_av = st.mean(irr_gen_comp)
irr_gen_throttle_av = st.mean(irr_gen_throttle)
p3r_loss_av = st.mean(p3r_loss)
copw_av = st.mean(copw)
copr_dry01_av = st.mean(copr_dry01)
copr_no_p_loss_av = st.mean(copr_no_p_loss)
copr2_av = st.mean(copr2)
eta_c_av = st.mean(eta_c)
power_draw_compressor_av = st.mean(power_draw_compressor)
power_draw_actual_compressor_av = st.mean(power_draw_actual_compressor)

T2w = np.array(T2w)
T1w = np.array(T1w)
h2r = np.array(h2r)
h1r = np.array(h1r)
p1 = np.array(p1)
power_comp_frac = st.mean(power_draw_compressor) / st.mean(power_draw)

# Values for TQ plot:
T3r_av = st.mean(T3r)
p2r_av = st.mean(p2)
p3r_av = st.mean(p3r)  # Pressure at state 3 after loss
p3r = np.array(p3r)
s2r_av = st.mean(s2r)
s3r_av = st.mean(s3r)
h2r_av = st.mean(h2r)
h3r_av = st.mean(h3r)  # Assuming wet saturated
T1w_av = st.mean(T1w)
T2w_av = st.mean(T2w)
mr_av = st.mean(mr)
mw_av = st.mean(mw)
T3_compare = cp.PropsSI('T','H',h3r_av,'P',p3r_av,"R134a")
T3_compare = np.array(T3_compare)
T3r = np.array(T3r)
print(T3_compare)

# Plotting TQ diagram and doing exergy analysis over the heat exchanger
pinch_diff = 10  # Defining FIXED pinch point temperature difference
# Create the figure and the initial plot
fig, ax = plt.subplots()
plt.subplots_adjust(left=0.25, bottom=0.35, top=0.85)  # Adjust top to make room for title
ax.set_xlabel('Heat Transfer (Q) [kJ]')
ax.set_ylabel('Temperature (°C)')
ax.grid(True)
fluid = "R134a"

P_water = 1e5
S1_water = cp.PropsSI('S', 'P', P_water, 'T', T1w_av, 'water')  # This calculates the entropy of the water at inlet
S2_water = cp.PropsSI('S', 'P', P_water, 'T', T2w_av, 'water')  # This calculates the entropy of water at the outlet
entropy_water = mw_av * (S2_water - S1_water)
entropy_ref = mr_av * (s3r_av - s2r_av)
entropy = entropy_ref + entropy_water  # Calculates overall entropy increase over the heat exchanger
entropy_vs_heat = entropy / ((h2r_av - h3r_av) * mr_av)  # This calculates entropy change per unit heat transfer

# Calculate the enthalpy and temperature profiles for the T-Q diagram
Hr_profile = np.linspace(h3r_av, h2r_av, 50)
Tr_profile = cp.PropsSI('T', 'H', Hr_profile, 'P', np.linspace(p3r_av, p2r_av, 50), fluid)  # Using consistent enthalpy profile approach

# Calculate the temperature difference between water and refrigerant
T_profile_water = np.linspace(T1w_av, T2w_av, 50)
T_diff = Tr_profile - T_profile_water
T_coeff = np.mean(T_diff)  # Find mean temperature difference between water and refrigerant over the cycle

# Recalculate the Heat Transfer Coefficient
Heat_coeff = ((h2r_av - h3r_av) * mr_av) / T_coeff  # Use mr_av instead of mr

# Plotting the T-Q diagram
fig, ax = plt.subplots()
plt.subplots_adjust(left=0.25, bottom=0.35, top=0.8)
ax.set_xlabel('Heat Transfer (Q) [kJ]')
ax.set_ylabel('Temperature (°C)')
ax.set_title(f'T-Q Diagram of the Heat Exchanger \n Entropy Generation = {entropy: .2f} [J/Ks] \n Entropy Generation per unit heat transferred = {entropy_vs_heat: .6f} [1/K]\nHeat Transfer Coefficient = {Heat_coeff: .2f} [J/Ks]')
ax.grid(True)

# Update plot with new values
T3 = cp.PropsSI('T','H',h3r_av,'P',p3r_av,"R134a")
print(T3)
plot_TQ_diagram(ax, p2r_av, p3r_av, h2r_av, h3r_av, T1w_av, T2w_av, pinch_diff, fluid, mw_av, mr_av)
ax.legend()
plt.show()


# Data for the table
table_data = [
    ["Calculated mass flow rate of air through evaporator (kg/s)", ma_av],
    ["Percentge pressure loss of refrigerant lost in the condensor(%)", cond_p_loss*100],
    ["Irreversible entropy generation due to compressor (J/kgK)", irr_gen_comp_av],
    ["Irreversible entropy generation due to throttle (J/kgK)", irr_gen_throttle_av],
    ["COP Based on water", copw_av],
    ["COP Based on Refrigerantassuming sat liquid @3", copr2_av],
    ["COP Based on Refrigerant assuming dryness 0.1@3", copr_dry01_av],
    ["COP Based on Refrigerant assuming no pressure loss", copr_no_p_loss_av],
    #['Power ouptut based on water (kW)' , st.mean(power_w)],
    #['Power output based on refrigerant (kW)', st.mean(power_r)],
    #['Percentage of Power Draw by Compressor (%)', power_comp_frac*100],
    ['Compressor Isentropic Efficiency', eta_c_av],
    ['Refrigerant mass flow assuming sat liquid @3', st.mean(mr)],
    ['Refrigerant mass flow assuming dryness 0.1@3', st.mean(mr_dry01)],
    ['Water volumetric flow rate', st.mean(qw)],
    ['Water mass flow rate', st.mean(mw)],
    ['Average T2w', st.mean(T2w)-273.15],
    ['Percentage pressure loss in the evaporator (%)', evap_p_loss*100],
    ['Predicted power drawn from compressor, mr (W) = ', power_draw_compressor_av],
    ['Measured power drawn from compressor, mr (W) = ', power_draw_actual_compressor_av],
    ['T3 of refrigerant (K) = ', T3r_av],
]

# Print the table using tabulate
print(tabulate(table_data, headers=["Metric", "Value"], tablefmt="grid"))
