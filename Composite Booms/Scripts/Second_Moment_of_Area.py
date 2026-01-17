import numpy as np
import math
import scipy.integrate as integrate
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True

h = 3e-2 #Flattened Height
w = 2e-3 #Web Height

Angle_Array_deg = np.linspace(5, 90, 200)
Angle_Array_rad = np.deg2rad(Angle_Array_deg)
Radius_Array = (h-w)/Angle_Array_rad

def yc_x (x,R,w):
    return (R**2 - (x-R)**2)**0.5 + w

def xc_y (y,R,w):
    return R - (R**2 - (y-w)**2)**0.5


Ixx_Array = np.zeros_like(Radius_Array)
Iyy_Array = np.zeros_like(Radius_Array)
Izz_Array = np.zeros_like(Radius_Array)

for i in range (len(Radius_Array)):
    R            = Radius_Array[i]
    theta        = Angle_Array_rad[i]
    Ixx = integrate.quad(lambda x: ((w + R*np.sin(theta))**3 - (yc_x(x, R, w)**3))/3, 0, R - R*np.cos(theta))
    Iyy = integrate.quad(lambda y: (xc_y(y, R, w)**3)/3, w, w + R * np.sin(theta))
    Izz = 2*Ixx[0] + 2*Iyy[0]
    Ixx_Array[i] = 2*Ixx[0]
    Iyy_Array[i] = 2*Iyy[0]
    Izz_Array[i] = Izz

# Find maximum Ixx
max_ixx_idx = np.argmax(Ixx_Array)
max_ixx_value = Ixx_Array[max_ixx_idx]
max_ixx_radius = Radius_Array[max_ixx_idx]
max_ixx_angle = Angle_Array_deg[max_ixx_idx]

print(f"Maximum Ixx Value is= {max_ixx_value} with Radius= {max_ixx_radius} and Subtended Angle= {max_ixx_angle}")

# Find maximum Iyy
max_iyy_idx = np.argmax(Iyy_Array)
max_iyy_value = Iyy_Array[max_iyy_idx]
max_iyy_radius = Radius_Array[max_iyy_idx]
max_iyy_angle = Angle_Array_deg[max_iyy_idx]

print(f"Maximum Iyy Value is= {max_iyy_value} with Radius= {max_iyy_radius} and Subtended Angle= {max_iyy_angle}")

# Check if they occur at the same radius and angle
if max_ixx_idx == max_iyy_idx:
    print("The maximum Ixx and Iyy values occur at the same radius and subtended angle.")
else:
    print("The maximum Ixx and Iyy values do not occur at the same radius and subtended angle.")


theta_plot = np.deg2rad(max_iyy_angle)

x_end = max_iyy_radius * (1 - np.cos(theta_plot))
x_array = np.linspace(0, x_end, 400)
y_array = yc_x(x_array, max_iyy_radius, w)

# convert to cm for plotting
x_cm = x_array * 1e2
y_cm = y_array * 1e2
x_cm = np.append([0], x_cm)
y_cm = np.append([0], y_cm)

plt.figure()
plt.plot( x_cm, y_cm, label="profile")
plt.plot(-x_cm, y_cm, label="mirrored")

plt.gca().set_aspect("equal", adjustable="box")
plt.grid(True)
plt.xlabel(r"$x$ [cm]")
plt.ylabel(r"$y$ [cm]")
plt.title(r"Geometry Maximising $I_{yy}$")

textbox_text = (
    rf"Radius: {max_iyy_radius*1e2:.4f} cm\\"
    rf"Subtended Angle: {max_iyy_angle:.2f}$^\circ$"
)
plt.text(
    0.05, 0.95, textbox_text,
    transform=plt.gca().transAxes,
    fontsize=10, va="top",
    bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8)
)

plt.show()

