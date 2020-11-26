from labutil.src.plugins.ising import *
import matplotlib.pyplot as plt

def problem1a():
    N_list = [8, 16, 32]
    T_list = [2, 2.5]
    slope_window = 1000

    n_eq = 100
    n_mc = 10

    all_E_data = []
    all_M_data = []

    for N in N_list:
        singleN_E_data = []
        singleN_M_data = []
        for T in T_list:
            E, M, E_eq, M_eq, imagelist = mc_run(N, n_eq, n_mc, T)
            singleN_E_data.append(E_eq)
            singleN_M_data.append(M_eq)

        all_E_data.append(singleN_E_data)
        all_M_data.append(singleN_M_data)

    N_len = len(N_list)
    T_len = len(T_list)

    fig1, axs1 = plt.subplots(N_len, T_len)

    for i in range(N_len):
        for j in range(T_len):
            plot_data(all_E_data[i][j], axs1[i, j])
            plot_rolling_slope(all_E_data[i][j], slope_window, axs1[i, j].twinx())
            axs1[i, j].set_xlabel('Steps')
            axs1[i, j].set_ylabel('Energy', color='tab:blue')
            axs1[i, j].set_title('N = ' + str(N_list[i]) + ', T = ' + str(T_list[j]))
    plt.subplots_adjust(hspace=0.5)
    plt.subplots_adjust(wspace=0.3)

def plot_data(data, ax):
    color='tab:blue'
    ax.plot(range(len(data)), data, color=color, linewidth=0.75)
    ax.tick_params(axis='y', labelcolor=color)


def plot_rolling_slope(data, window, ax):
    color='tab:red'

    fit_data = []
    for i in range(len(data) - window + 1):
        fit = np.polyfit(range(window), data[i:i + window], 1) # fit a line
        fit_data.append(fit[0])
        #fit_data.append(np.std(data[i:i + window])) # fluctuations?

    ax.plot([np.floor(window/2) + i for i in range(len(fit_data))], fit_data, color=color)
    ax.ticklabel_format(style='sci')
    ax.tick_params(axis='y', labelcolor=color)
    ax.axhline(0, color='black')
    ax.set_ylabel('Slope', color=color)
    #ax.set_ylabel('Std Dev', color=color)


def problem1b():
    Tstart = 1.0
    Tstop = 3.5
    numTs = 15

    N = 16
    n_eq = 75
    n_mc = 1000

    T_list = np.linspace(Tstart, Tstop, numTs)
    M_list = []
    for T in T_list:
        E, M, E_eq, M_eq, imagelist = mc_run(N, n_eq, n_mc, T)
        M_list.append(abs(np.mean(M)))

    fig, ax = plt.subplots(1, 1)
    ax.plot(T_list, M_list, marker='o')


def problem1b_multiple_runs():
    # throwing out more data but hopefully random initializations help
    # initial N=12 didn't look to neat so expanded to other sizes, still not neat, explain
    Tstart = 1.0
    Tstop = 3.5
    numTs = 15
    num_runs = 1
    N_list = [3, 12]

    n_eq = 75
    n_mc = 1000

    fig, ax = plt.subplots(1, 1)

    N_list = [3, 12]
    T_list = np.linspace(Tstart, Tstop, numTs)
    last_images = []
    for N in N_list:
        M_list = []
        for T in T_list:
            M_runs = []
            for i in range(num_runs):
                E, M, E_eq, M_eq, imagelist = mc_run(N, n_eq, n_mc, T)
                M_runs.append(abs(np.mean(M)))
            M_list.append(np.mean(M_runs))

        last_images.append(imagelist[-1]) # saves last image on last temp, need to save first and last temp
        ax.plot(T_list, M_list, marker='o', label=str(N))

    ax.legend()

    fig1, axs1 = plt.subplots(1, len(N_list))
    for i in range(len(N_list)):
        axs1[i].imshow(last_images[i])


def problem1b_n_mc_convergence():
    Tstart = 1.0
    #Tstart = 2.0
    Tstop = 3.5
    #Tstop = 2.5
    numTs = 15

    N = 16
    n_eq = 100
    n_mc_list = [1000]

    T_list = np.linspace(Tstart, Tstop, numTs)
    M_list_list = []
    images = []
    for n_mc in n_mc_list:
        print(n_mc)
        M_list = []
        for T in T_list:
            print(T)
            E, M, E_eq, M_eq, imagelist = mc_run(N, n_eq, n_mc, T)
            M_list.append(abs(np.mean(M)))
            if n_mc == n_mc_list[-1]:
                if T == T_list[4]:
                    images.append(imagelist[-1])
                elif T == T_list[-5]:
                    images.append(imagelist[-1])


        M_list_list.append(M_list)

    print(T_list)
    print(M_list)
    fig, ax = plt.subplots(1, 1)
    for i in range(len(n_mc_list)):
        ax.plot(T_list, M_list_list[i], marker='o', label=str(n_mc_list[i]))
    ax.legend()
    ax.set_xlabel('T')
    ax.set_ylabel('<M>')
    ax.set_title('Magnetization vs Temperature')

    fig_image, (ax_left, ax_right) = plt.subplots(1, 2)
    ax_left.imshow(images[0])
    ax_left.set_title('T = ' + str(T_list[4]))
    ax_right.imshow(images[-1])
    ax_right.set_title('T = ' + str(T_list[-5]))


def problem2a():
    Trange = 3
    numTs = 20
    T_list = np.linspace(2.27 - Trange/2, 2.27 + Trange/2, numTs)
    N = 16
    num_runs = 1

    n_eq = 75
    n_mc = 1000

    E_list = []

    for T in T_list:
        print(T)
        E_runs = []
        for i in range(num_runs):
            print(i)
            E, M, E_eq, M_eq, imagelist = mc_run(N, n_eq, n_mc, T)
            E_runs.append(np.mean(E))
        E_list.append(np.mean(E_runs))

    fig, ax = plt.subplots(1, 1)
    ax2 = plot_and_derivative(T_list, E_list, ax)
    ax.set_xlabel('T')
    ax.set_ylabel('E', color='tab:blue')
    ax.set_title('Energy and Heat Capacity vs Temperature')

    ax2.set_ylabel('c', color='tab:red')


def plot_and_derivative(x, y, ax):
    n = len(x)
    center_x = []
    dydx = []
    for i in range(n - 1):
        center_x.append((x[i] + x[i + 1])/2)
        dydx.append((y[i + 1] - y[i])/(x[i + 1] - x[i]))

    color1 = 'tab:blue'
    ax.plot(x, y, color=color1, marker='o')
    #ax.plot(x, y, color=color1, marker='o', label='Energy')
    ax.tick_params(axis='y', labelcolor=color1)

    color2 = 'tab:red'
    ax2 = ax.twinx()
    ax2.plot(center_x, dydx, color=color2, marker='s')
    #ax2.plot(center_x, dydx, color=color2, marker='s', label='c, derivative')
    ax2.tick_params(axis='y', labelcolor=color2)

    return ax2

def problem2a_and_b():
    Trange = 3
    numTs = 20
    T_list = np.linspace(2.27 - Trange/2, 2.27 + Trange/2, numTs)
    N = 16
    num_runs = 5

    n_eq = 75
    n_mc = 1000

    E_list = []
    c_list = []

    testEnergies = [] #testing

    for T in T_list:
        print(T)
        E_runs = []
        c_runs = []
        for i in range(num_runs):
            print(i)
            E, M, E_eq, M_eq, imagelist = mc_run(N, n_eq, n_mc, T)
            E_runs.append(np.mean(E))
            c_runs.append(np.power(T, -2)*(np.mean([np.power(element, 2) for element in E]) - np.power(np.mean(E), 2)))
        testEnergies.append(E)  # testing
        E_list.append(np.mean(E_runs))
        c_list.append(np.power(N, 2)*np.mean(c_runs))

    fig, ax = plt.subplots(1, 1)
    ax2 = plot_and_derivative(T_list, E_list, ax)
    color3 = 'tab:red'
    ax2.plot(T_list, c_list, color=color3, marker='^', label='c, fluctuations')

    ax.set_xlabel('T')
    ax.set_ylabel('E')
    ax.set_title('E and c vs T')
    ax2.set_ylabel('c')
    ax.legend()
    ax2.legend()

    testfig, testax = plt.subplots(1,1) #testing
    for i in range(len(T_list)):
        testax.plot(range(len(testEnergies[i])), testEnergies[i], label=str(T_list[i]))
    testax.legend()
    testax.set_xlabel('steps')
    testax.set_ylabel('E')


def problem2c():
    Trange = 2.5
    numTs = 20
    T_list = np.linspace(2.27 - Trange / 2, 2.27 + Trange / 2, numTs)

    N_list = [4, 8, 16, 32]
    n_eq = 75
    n_mc = 1000

    E_list_list = []
    c_list_list = []

    for N in N_list:
        print(N)
        E_list = []
        c_list = []
        for T in T_list:
            print(T)
            E, M, E_eq, M_eq, imagelist = mc_run(N, n_eq, n_mc, T)
            E_list.append(np.mean(E))
            c_list.append(np.power(N, 2)*np.power(T, -2) * (np.mean([np.power(element, 2) for element in E]) - np.power(np.mean(E), 2)))
        E_list_list.append(E_list)
        c_list_list.append(c_list)

    fig, (ax_left, ax_right) = plt.subplots(1, 2, sharey=True)

    for i in range(len(N_list)):
        plot_only_derivative(T_list, E_list_list[i], ax_left)
        ax_right.plot(T_list, c_list_list[i], marker='o', label=str(N_list[i])+'x'+str(N_list[i]))

    ax_right.legend()
    ax_right.set_xlabel('T')
    ax_right.set_title('c vs T, fluctuation')
    ax_left.set_xlabel('T')
    ax_left.set_ylabel('c')
    ax_left.set_title('c vs T, derivative')


def plot_only_derivative(x, y, ax):
    n = len(x)
    center_x = []
    dydx = []
    for i in range(n - 1):
        center_x.append((x[i] + x[i + 1])/2)
        dydx.append((y[i + 1] - y[i])/(x[i + 1] - x[i]))
    ax.plot(center_x, dydx, marker='o')

def problem3():
    Trange = 2.5
    numTs = 5
    T_list = np.linspace(2.27 - Trange / 2, 2.27 + Trange / 2, numTs)

    N_list = [8, 16]
    n_eq = 75
    n_mc = 300

    X_list_list = []

    for N in N_list:
        print(N)
        X_list = []
        for T in T_list:
            print(T)
            E, M, E_eq, M_eq, imagelist = mc_run(N, n_eq, n_mc, T)
            #test_fig = plt.figure()
            #test_ax = test_fig.add_subplot(1, 1, 1)
            #test_ax.plot(range(len(M)), M)
            X_list.append(np.power(N, 2)*np.power(T, -1)*(np.mean([np.power(element, 2) for element in M]) - np.power(np.mean(M), 2)))
        X_list_list.append(X_list)

    fig, ax = plt.subplots(1, 1)

    for i in range(len(N_list)):
        ax.plot(T_list, X_list_list[i], marker='o', label=str(N_list[i])+'x'+str(N_list[i]))

    ax.set_yscale('log')
    ax.legend()


def problem3_multirun():
    Tstart = 2.0
    Tstop = 3.5
    numTs = 10
    T_list = np.linspace(Tstart, Tstop, numTs)
    num_runs = 5

    N_list = [4, 8, 16, 32]
    n_eq = 75
    n_mc = 1000

    X_list_list = []

    for N in N_list:
        print(N)
        X_list = []
        for T in T_list:
            print(T)
            X_runs = []
            for i in range(num_runs):
                print(i)
                E, M, E_eq, M_eq, imagelist = mc_run(N, n_eq, n_mc, T)
                X_runs.append(np.power(N, 2)*np.power(T, -1)*(np.mean([np.power(element, 2) for element in M]) - np.power(np.mean(M), 2)))
            X_list.append(np.mean(X_runs))
        X_list_list.append(X_list)

    fig, ax = plt.subplots(1, 1)
    for i in range(len(N_list)):
        ax.plot(T_list, X_list_list[i], marker='o', label=str(N_list[i])+'x'+str(N_list[i]))
    ax.axvline(2.27, color='lightgray', linestyle='--')

    ax.set_yscale('log')
    ax.legend()
    ax.set_xlabel('T')
    ax.set_ylabel('X')
    ax.set_title('Susceptibility vs Temperature for different lattice sizes')

if __name__ == "__main__":
    #problem1a() # use for 1a - under 1 min
    #problem1b() # looks good, do multiple n_mcs though
    #problem1b_multiple_runs() #ignore
    #problem1b_n_mc_convergence() # use for 1b - under 2 mins
    #problem2a() # use for 2a - under 10 mins
    problem2a_and_b() # use for 2b - same as prev
    #problem2c()
    #problem3_multirun()
    plt.show()