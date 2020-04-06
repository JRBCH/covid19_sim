def sim_fweuler(S_init, I_init, R_init, gamma, beta, T):


    N = S_init + I_init + R_init
    trange = np.arange(0,T)
    S = np.zeros(trange.size)
    I = np.zeros(trange.size)
    R = np.zeros(trange.size)

    S[0] = S_init
    I[0] = I_init
    R[0] = R_init

    for t in trange[1:]:
        S[t] = S[t-1] - (beta * I[t-1] * S[t-1])/N
        I[t] = I[t-1] + (beta * I[t-1] * S[t-1])/N - gamma * I[t - 1]
        R[t] = R[t-1] + gamma * I[t-1]

    return S,I,R
