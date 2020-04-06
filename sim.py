# COVID19 simulations

import numpy as np
from scipy.integrate import solve_ivp

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib.style
import matplotlib as mpl
mpl.style.use('default')


# PLOTTING OPTIONS
colors = {'S': '#1b9e77', 'I': '#d95f02', 'R': '#7570b3',
          'Slight' : '#66c2a5', 'Ilight' : '#fc8d62', 'Rlight': '#8da0cb'}
lw = 2

def sir_model(t, y, beta, gamma, N):
    """
    Defines the differential equations for the SIR model system.

    Arguments:
        t :  time
        y :  vector of state variables:
                    y = [S, I, R]
        beta, gamma, N : model parameters
    """

    S, I, R = y

    # Output array of form f = (S', I', R')
    f = [- (beta * I * S)/N,
           (beta * I * S)/N - gamma * I,
           gamma * I]

    return f

def draw_plot(S,I,R,time):

    fig = plt.figure(constrained_layout=True, figsize=(10, 6))
    spec = gridspec.GridSpec(ncols=2, nrows=2, figure=fig)
    ax1 = fig.add_subplot(spec[0, 0])
    ax2 = fig.add_subplot(spec[0, 1])
    ax3 = fig.add_subplot(spec[1, 0])
    ax4 = fig.add_subplot(spec[1, 1])

    xaxis1 = time
    ax1.plot(xaxis1, S/N, color = colors['S'], label = 'Healthy', lw=lw)
    ax1.plot(xaxis1, I/N, color = colors['I'], label = 'Infectious', lw=lw)
    ax1.plot(xaxis1, R/N, color = colors['R'], label = 'Removed', lw=lw)
    ax1.set_xlabel('Time (days)')
    ax1.set_ylabel('Proportion of population')
    ax1.legend(frameon=False)
    # Hide the right and top spines
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)

    ax2.fill_between(xaxis1, 0, R, color = colors['Rlight'])
    ax2.hlines(np.max(R), 0, np.max(xaxis1), colors='black', linestyles='dashed')
    ax2.text(0, 0.8*np.max(R), str(int(np.max(R))) + ' total infections \n('+str(np.round(np.max(R)/N * 100, 2))+' % of population)')
    ax2.set_xlabel('Time (days)')
    ax2.set_ylabel('Cumulative cases')
    # Hide the right and top spines
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)


    # calculate deltaI
    deltaI = [I[i]- I[i-1] for i in np.arange(1,np.size(I))]

    ax4.fill_between(xaxis1[1:], 0, deltaI, color = colors['Ilight'])
    ax4.set_xlabel('Time (days)')
    ax4.set_ylabel(r'$\Delta$ Infectious population')
    # Hide the right and top spines
    ax4.spines['right'].set_visible(False)
    ax4.spines['top'].set_visible(False)

    return fig


# TO DO
# Add death rate, hospital bed threshold, interventions

if __name__ == "__main__":

    # Starting values and Parameters
    N = 82790000
    R_init = 30284
    I_init = 69848
    S_init = N - I_init - R_init
    gamma = 0.025  # rate of recovery / deat
    beta = 0.06  # rate of infection
    T = 500  # Number of days to simulate

    t_span = (0.0, T)

    sol = solve_ivp(sir_model, t_span, [S_init, I_init, R_init],
                    args=(beta, gamma, N),
                    first_step=1, max_step=1)

    St, It, Rt = sol['y']
    time = sol['t']

    draw_plot(St, It, Rt, time)
    plt.show()

