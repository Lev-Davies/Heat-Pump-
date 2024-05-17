import matplotlib.pyplot as plt
import CoolProp.CoolProp as cp
import numpy as np 
import statistics as st
from tabulate import tabulate

# Initiate lists for values, number is circuit locations,w is wataer, a is air , r is refrigerant
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
#pick the text file to use
file_name = "HP_May10_05"

#add values from files to list 

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
                    tim.append(float(vals[0]))  # Convert values to floart and assign to corresponding list 
                    T1w.append(float(vals[1])+273.15)  # Convert the second value to floating point and append to T1w[]
                    T2w.append(float(vals[2])+273.15)  # convert temp to K
                    T1a.append(float(vals[3])+273.15)
                    T2a.append(float(vals[4])+273.15)
                    T1r.append(float(vals[5])+273.15)
                    T2r.append(float(vals[6])+273.15)
                    T3r.append(float(vals[7])+273.15)
                    T4r.append(float(vals[8])+273.15)
                    p1.append(float(vals[9])*10**5) # convert pressures to Pa
                    p2.append(float(vals[10])*10**5)
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

#set up lists for values, h is corresponding enthalpy, m is mass flow, p is power 

h1w = [] 
I_no_comp = 3.2740721739130434 # calculated current when the fan is on full and compressor is not running.
h2w = []
hdifw = []
h1r = []
h2r = []
h3r= []
h4r = []
mw  = []
mr = [] 
ma = []
pw = []
pw2 =[]
pr =[]
copw = []
copw2 =[] 
copr = []
wc = []
pdraw = [] # power drawn by system
pdraw_compressor = [] # power drawn by compressor
copr2 = []
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


for i in range(700):
    #calculate enthalpies using coolprop
    h1w.append(cp.PropsSI ('H','T',T1w[i],'Q',0.0,"Water")) 
    h2w.append(cp.PropsSI ('H','T',T2w[i],'Q',0.0,"Water"))
    hdifw.append(4200*(T2w[i]-T1w[i]))

    h1r.append(cp.PropsSI ('H','P|gas', p1[i] ,'T',T1r[i],"R134a")) # using t and p here gives sevrely oscillating values of h1 
    h2r.append(cp.PropsSI ('H','T',T2r[i],'P',p2[i],"R134a"))
    h3r.append(cp.PropsSI ('H','T',T3r[i],'Q',0,"R134a")) # assuming wet saturated on exit from the condenser
    h4r.append(h3r[i])
    
    pdraw.append(I[i]*240*0.98)
    pdraw_compressor.append((I[i]-I_no_comp)*240*0.98) # now only considering compressor work and using power factor

    mw.append(qw[i]/60) # calculate water mass flow (assume rho = 1000)
    mr.append(pdraw_compressor[i]/(h2r[i]-h1r[i])) 

    pw.append(mw[i]*(hdifw[i])) # calculate heat output base on enthalpy gained by water
    
    copw.append(pw[i]/pdraw[i]) #calculate COP based on water

    #copr.append(pr[i]/pdraw[i])
    copr2.append((h2r[i] - h3r[i])/(h2r[i] - h1r[i])) #calculating COP purely thermodynamically

    # calculating air mass flow rate 
    ma.append(mr[i]*(h4r[i] - h1r[i])/(1005*(T2a[i] - T1a[i])))


#looking at irreversible entropy generation and pressure losses.
for i in range(700):
    s3r.append(cp.PropsSI ('S','T',T3r[i],'Q',0,"R134a"))
    s2r.append(cp.PropsSI ('S','T',T2r[i],'P',p2[i],"R134a"))
    s1r.append(cp.PropsSI ('S','P|gas', p1[i] ,'T',T1r[i],"R134a"))
    s4r.append(cp.PropsSI ('S','P', p1[i] ,'T',T4r[i],"R134a"))
    p3r.append(cp.PropsSI ('P','Q',0.0 ,'T',T3r[i],"R134a"))
    p3r_loss.append(p2[i] - p3r[i])
    irr_gen_comp.append(s2r[i] - s1r[i])
    irr_gen_throttle.append(s4r[i] - s3r[i])


ma_av = st.mean(ma)
irr_gen_comp_av = st.mean(irr_gen_comp)
irr_gen_throttle_av = st.mean(irr_gen_throttle)
p3r_loss_av = st.mean(p3r_loss)
copw_av = st.mean(copw)
#copr_av = st.mean(copr)
copr2_av = st.mean(copr2)


T2w = np.array(T2w)
T1w = np.array(T1w)
h2r = np.array(h2r)
h1r = np.array(h1r)
p1 = np.array(p1)

# Data for the table
table_data = [
    ["Calculated mass flow rate of air through evaporator (kg/s)", ma_av],
    ["Pressure of refrigerant lost in the heat exchanger (bar)", p3r_loss_av / 1e5],
    ["Irreversible entropy generation due to compressor (J/kgK)", irr_gen_comp_av],
    ["Irreversible entropy generation due to throttle (J/kgK)", irr_gen_throttle_av],
    ["COP Based on water", copw_av],
    ["COP Based on Refrigerant", copr2_av]
]

# Print the table using tabulate
print(tabulate(table_data, headers=["Metric", "Value"], tablefmt="grid"))