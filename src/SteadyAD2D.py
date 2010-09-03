'''
@author: BYK, Deli
@contact: gulen.ilker@hotmail.com, madbyk@gmail.com
@summary: 
	A steady 2D advection-diffusion FEM solver in Python 2.6
	using NumPy, SciPy and MatPlotLib
@version: 1.5
'''

if __name__ == "__main__":
	from psetup import get_problem_data
	from gsystem import calc_global
	from solveproc import post_process, solve_system
	from time import time
	problem_data = get_problem_data()
	t = time()
	K, F = calc_global(problem_data)
	solution = solve_system(K, F)
	t = time() - t
	#print solution
	print ("Total run time: {0} seconds.".format(t))
	post_process(problem_data, solution)
