'''
@author: Deli, BYK
@contact: gulen.ilker@hotmail.com, madbyk@gmail.com
@summary:
    Provides solving and post-processing functions.
    Solving function uses SciPy's sparse solving function from sparse linear algebra module(sparse.spsolve).
    Post process shapes the coordinates and results into a proper form to be plotted,
    refines the plotting range and approximates mid values using node neighbors
    and then produces a smooth contour plot.
@version: 1.6
'''

def solve_system(K, F):
    """
    Solves the K * x = F system using the linear sparse matrix solver.
    The spsolve function used here can be replaced with qmr or some other
    non-linear solver function for non-linear problems such as N-S problems.
    """
    from scipy.sparse import linalg

    print("Solving system...")
    return linalg.spsolve(K, F)

def save_solution(file_name, solution):
    from json import dump

    print(" * Writing output file...")
    output_file = open(file_name, "w")
    dump({"T": solution.tolist()}, output_file)
    output_file.close()

def plot_solution(problem_data, solution):
    from numpy import linspace
    from matplotlib.pylab import contour, colorbar, contourf, xlabel, ylabel, title, show
    from matplotlib.mlab import griddata

    print(" * Preparing for plotting...")
    NN = problem_data["NN"]

    #Extract node coordinates seperately for plotting
    x, y = [0] * NN, [0] * NN
    for i, node in enumerate(problem_data["nodes"]):
        x[i] = node[0]
        y[i] = node[1]

    #refine the contour plot mesh for a "smoother" image, generate a 200*200 grid
    xi = linspace(min(x), max(x), 200)
    yi = linspace(min(y), max(y), 200)
    #approximate the mid values from neighbors
    zi = griddata(x, y, solution, xi, yi)

    print(" * Plotting...")
    #plot the contour lines with black
    contour(xi, yi, zi, 15, linewidths = 0.5, colors = 'k')
    #plot the filled contour plot
    plot = contourf(xi, yi, zi, 15, antialiased = True)

    colorbar(plot, format = "%.3f").set_label("T")
    xlabel('X')
    ylabel('Y')
    title("Contour plot of T values for {0}".format(problem_data["title"]))

    show()

def post_process(problem_data, solution, arguments):
    """
    Performs the necessary post processing operations on the "solution"
    using the problem_data such as contour plotting, SV calculating etc.
    """
    print("Post processing...")

    if not arguments.dontsave:
        save_solution(problem_data["output"], solution)

    if not arguments.dontplot:
        plot_solution(problem_data, solution)
