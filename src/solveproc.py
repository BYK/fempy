'''
@author: Deli, BYK
@contact: gulen.ilker@hotmail.com, madbyk@gmail.com
@summary:
	Provides solving and post-processing functions.
	Solving function uses SciPy's sparse solving function from sparse linear algebra module(sparse.spsolve).
	Post process shapes the coordinates and results into a proper form to be plotted,
	refines the plotting range and approximates mid values using node neighbors
	and then produces a smooth contour plot.
@version: 1.5
'''
from numpy import linspace
import matplotlib.pylab as pylab
import matplotlib.cm as cm
from matplotlib.mlab import griddata
from scipy.sparse import linalg
from json import dump

def solve_system(K, F):
	"""
	Solves the K * x = F system using the linear sparse matrix solver.
	The spsolve function used here can be replaced with qmr or some other
	non-linear solver function for non-linear problems such as N-S problems.
	"""
	print "Solving system..."
	return linalg.spsolve(K, F)

def post_process(problem_data, solution):
	"""
	Performs the necessary post processing operations on the "solution"
	using the problem_data such as contour plotting, SV calculating etc.
	"""
	output_file = open(problem_data["filename"] + "_output.json", "w")
	dump({"T": solution.tolist()}, output_file, indent = 2)
	output_file.close()

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

	#plot the contour lines with black
	pylab.contour(xi, yi, zi, 15, linewidths = 0.5, colors = 'k')
	#plot the filled contour plot
	plot = pylab.contourf(xi, yi, zi, 15, antialiased = True)

	pylab.colorbar(plot, format = "%.3f").set_label("T")
	pylab.xlabel('X')
	pylab.ylabel('Y')
	pylab.title("Contour plot of T values for {0}".format(problem_data["title"]))

	pylab.show()
