'''
@author: BYK, Deli
@contact: gulen.ilker@hotmail.com, madbyk@gmail.com
@summary: 
	A steady 2D advection-diffusion FEM solver in Python 2.6
	using NumPy, SciPy and MatPlotLib
@version: 1.5
'''

if __name__ == "__main__":
	import argparse

	from psetup import get_problem_data
	from gsystem import calc_global
	from solveproc import post_process, solve_system
	from time import time

	parser = argparse.ArgumentParser(description = 'Solves steady and 2D advection/diffusion problems using finite elements method.')
	parser.add_argument('-i', '--input', default = '', help = 'Input file path.')
	parser.add_argument('-o', '--output', default = '', help = 'Output file path.')
	arguments = parser.parse_args()

	problem_data = get_problem_data(arguments.input, arguments.output)

	#Exclude input reading time from total time
	t = time()

	#Calculate the system
	K, F = calc_global(problem_data)

	#Solve the system
	solution = solve_system(K, F)

	#Calculate the total running time
	t = time() - t

	print("Total run time: {0} seconds.".format(t))

	#Exclude the post processing time from total time
	post_process(problem_data, solution)
