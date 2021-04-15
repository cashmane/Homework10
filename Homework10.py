import numpy as np
import matplotlib.pyplot as plt

#Fix this
G = 6.6738e-11 #Newtons KG**-2 m**2
M = 1.9891e30 #Mass of sun in KG
AU = 1.496e11 #1 astronomical unit in meters

def f(r):
    """derivative function to pass to rk4
    pass state vector r = [x, y, vx, vy]"""

    x, y, vx, vy = r
    rcubed = np.sqrt(x**2 + y**2)**3

    fx = vx #return 1st parameter
    fy = vy #return 2nd parameter
    fvx = -1*G*M*(x/rcubed) #return 3rd parameter
    fvy = -1*G*M*(y/rcubed) #return last parameter
    return np.array([fx, fy, fvx, fvy])

def rk4_step(r=None, h=None, f=None):
    '''returns one step of rk4 method, used in both fixed and adaptive functions.
        r = initial state vector
        h = step size
        f = derivative of function'''
    k1 = h*f(r)
    k2 = h*f(r+0.5*k1)
    k3 = h*f(r+0.5*k2)
    k4 = h*f(r+k3)
    return ((k1+2*k2+2*k3+k4)/6)

def run_rk4_fixed(initial_state=None, initial_h=None, tmax=None):
    '''Runs RK4 with a fixed step size.
        initial_state = vector with positionX, positionY, velocityX, velocityY
        initial_h = step size for rk4
        tmax = the amount of iterations to calculate'''
    r = initial_state
    h = initial_h
    xpoints = []
    ypoints = []
    t=0
    while t<tmax:
        r = r + rk4_step(r=r, h=h, f=f)
        t+=h
        xpoints.append(r[0])
        ypoints.append(r[1])
    return np.array([xpoints, ypoints])
    

def run_rk4_adaptive(initial_state=None, initial_h=None, tmax=None):
    '''Rk4 with an adaptive step size.
        initial_state = vector with positionX, positionY, velocityX, velocityY
        initial_h = initial step size for rk4
        tmax = the amount of iterations to calculate'''
    r = initial_state
    h = initial_h
    t=0
    xpoints = []
    ypoints = []
    iterations = 0
    while t<tmax:
        #Do one large step
        r1 = r + rk4_step(r=r, h=2*h, f=f)
        #Do two small steps
        r2 = r + rk4_step(r=r, h=h, f=f)
        r2 = r2 + rk4_step(r=r2,h=h, f=f)

        #calculate value of rho
        ex = (1/30)*(r2[0]-r1[0])
        ey = (1/30)*(r2[1]-r1[1])
        rho = (30*h*delta)/(np.sqrt(ex**2+ey**2)+1e-5)

        #calculate new values of t, h, r
        #update points if appropriate
        if rho>=1.0:
            t = t+h
            h = h*rho**(1/4)
            r = r + rk4_step(r=r, h=h, f=f)
            xpoints.append(r[0])
            ypoints.append(r[1])
            iterations += 1
            pass
        else:
            h = h*rho**(1/4)
            pass
    return np.array([xpoints, ypoints])

if __name__ == "__main__":
    h0 = 1.0e4 #initial step size
    tmax = 3.0e9 #total time

    delta = 1e3/(365.25*24*3600) #meters accuracy per second

    x0, y0 = 4e12, 0 #starting pos, 4 billion kilometers
    vx0, vy0 = 0, 500 #starting velocity, m/s
    r0 = np.array([x0, y0, vx0, vy0])

    xpos, ypos = run_rk4_fixed(initial_state=r0, initial_h=h0, tmax=tmax)
    #Make the plot of fixed step size
    plt.plot(xpos/AU, ypos/AU, alpha = 0.5)
    plt.plot(xpos/AU, ypos/AU, 'k.')
    plt.title("Fixed Step Size")
    plt.show()
    plt.clf()
    print('For the fixed step size, I used a step size of 1e4, and a tmax of 3e9 to get an accurate solution.')
    print('It took between 5 and 10 seconds to run with these settings')
    
    xpos, ypos = run_rk4_adaptive(initial_state=r0, initial_h=h0, tmax=tmax)
    #Make the plot of adaptive step size
    plt.plot(xpos/AU, ypos/AU, alpha = 0.5)
    plt.plot(xpos/AU, ypos/AU, 'k.')
    plt.title("Adaptive Step Size")
    plt.show()
    print('The adaptive method was much faster, taking less than a second. The step size does noticeably change along the curve.')
    
