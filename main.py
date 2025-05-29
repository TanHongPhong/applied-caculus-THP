import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar
import numpy as np

def lens_sag_numeric(y, R, K):
    sqrt_term = np.sqrt(1 - (1 + K) * y**2 / R**2)
    sag = y**2 / (R * (1 + sqrt_term))
    return sag

def lens_sag_inverted(y, R, K ):
    return -lens_sag_numeric(y, R, K) + D

def lens_line_intersection(flip, a, b, K, R):
    def f(x):
        y = a * x + b
        radical = 1 - ((1 + K) * y**2) / R**2
        if radical < 0:
            return np.nan
        sag_x = y**2 / (R * (1 + np.sqrt(radical)))
        return sag_x - x
    def flip_f(x):
        y = a * x + b
        radical = 1 - ((1 + K) * y**2) / R**2
        if radical < 0:
            return np.nan
        sag_x = y**2 / (R * (1 + np.sqrt(radical)))
        return  - sag_x - x + 0.2
    sol = root_scalar(
        f if flip==0 else flip_f,
        bracket=[0, 0.2],
        method='brentq')
    if sol.converged:
        x_intersect = sol.root
        y_intersect = a * x_intersect + b
    return (x_intersect,y_intersect)

def line_from_point_and_angle(point, angle):
    x0, y0 = point
    
    # Handle special case: angle 90° or 270° (vertical line)
    if np.isclose(np.cos(angle), 0):
        return None, None  # Vertical line does not have the form y = ax + b

    a = np.tan(angle)
    b = y0 - a * x0
    return a, b

def lens_sag_derivative_numeric(y, R, K):
    sqrt_term = np.sqrt(1 - (1 + K) * y**2 / R**2)
    numerator = 2 * y * R * (1 + sqrt_term) + (1 + K) * y**3 / (R * sqrt_term)
    denominator = (R * (1 + sqrt_term))**2
    dz_dy = numerator / denominator
    theta_rad = np.arctan(dz_dy)
    return abs(theta_rad)

def caculate_refraction(denta, i, r, n1, n2):
    if r== None:
        r = (np.arcsin((n1/n2)*np.sin(i+denta)) - denta)
        return r if theta > 0 else -r
    elif i== None:
        i = -(-np.arcsin((n2/n1)*np.sin(-r+denta)) + denta)
        return i if theta > 0 else -i

# Function to draw line from point A to B
def draw_line(A=(0, 0), B=(0,0),name="AB",color='b'):
    x_vals = [A[0], B[0]]
    y_vals = [A[1], B[1]]
    plt.plot(x_vals, y_vals, color)

P = (-10,0)
plt.plot(P[0], P[1], 'ro')  # 'r' is red color, 'o' is circle marker

O1 = (-5,0)
O2 = (4.8,0)
D = 0.2

R = R1 = R2 = 10
k=0
theta_deg = 7
theta = np.deg2rad(theta_deg) 
n1=1
n2=1.5

a,b = line_from_point_and_angle(P,theta)

A = lens_line_intersection( flip=0, a=a, b=b, K = k, R = R)

denta1 = lens_sag_derivative_numeric(A[1],R,k)

r = caculate_refraction(denta1, i=theta, r=None, n1=1, n2=1.5)

c,d = line_from_point_and_angle(A,r)

B = lens_line_intersection( flip=1, a=c, b=d, K = k, R = R)

denta2 = lens_sag_derivative_numeric(B[1],R,k)

i = caculate_refraction(denta2 , i= None, r = r, n1=1.5, n2=1)

e,f = line_from_point_and_angle(B,i)

P1 = (-f/e,0)

draw_line(P,A,"PA")
draw_line(A,B,"AB")
draw_line(B,P1,"BP1")
draw_line(P,P1)

y_vals = np.linspace(-1.40854, 1.4054, 500)
# Calculate sag (z) for both surfaces
x_vals = lens_sag_numeric(y_vals, R, k)
x_inv_vals = lens_sag_inverted(y_vals, R, k)
# Plot curves: swap x and y -> x=z, y=y
plt.plot(x_vals, y_vals, label='Original lens surface', color='orange')
plt.plot(x_inv_vals, y_vals, label='Inverted lens surface', color='orange')



theta_list_deg = np.linspace(1, 7, 15)  # from 1 to 7 degrees
theta_list = np.deg2rad(theta_list_deg)

def simulate(K):
    x_positions = []
    for theta in theta_list:
        try:
            # Calculate line from P at angle theta
            a, b = line_from_point_and_angle(P, theta)

            # Intersection with first spherical surface
            A = lens_line_intersection(flip=0, a=a, b=b, K=K, R=R)
            if A is None: continue
            denta1 = lens_sag_derivative_numeric(A[1], R, K)
            r1 = caculate_refraction(denta1, i=theta, r=None, n1=n1, n2=n2)

            # Refracted ray line
            c, d = line_from_point_and_angle(A, r1)

            # Intersection with second spherical surface
            B = lens_line_intersection(flip=1, a=c, b=d, K=K, R=R)
            if B is None: continue
            denta2 = lens_sag_derivative_numeric(B[1], R, K)
            i2 = caculate_refraction(denta2, i=None, r=r1, n1=n2, n2=n1)

            # Final exiting ray
            e, f = line_from_point_and_angle(B, i2)
            P1 = (-f/e,0)
            x_positions.append(P1[0])
            print(P1)
            plt.plot(P1[0], P1[1], 'ro')  # 'r' is red color, 'o' is circle marker

            # Draw lines for ray path
            draw_line(P,A,"PA")
            draw_line(A,B,"AB")
            draw_line(B,P1,"BP1")
            draw_line(P,P1)
        except Exception as e:
            continue  # skip errors during simulation steps
    # if len(x_positions) < 2:
    #     return None
    # return abs(max(x_positions) - min(x_positions))

# K_values = np.linspace(6, 7, 100)
# dist_list = [simulate(K) for K in K_values]
# valid_data = [(i, dist) for i, dist in enumerate(dist_list) if dist is not None]
# min_val = min(valid_data, key=lambda x: x[1])
# print(K_values)
# print(f"Min distance: {min_val[1]:.6f} at index {min_val[0]}")
simulate(k)

plt.axis('equal')  # Correct aspect ratio between axes
plt.grid(True)
plt.title("Simulation of refraction through two spherical surfaces")
plt.legend()
plt.show()
