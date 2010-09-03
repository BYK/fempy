'''
@author: Deli, BYK
@contact: gulen.ilker@hotmail.com, madbyk@gmail.com
@summary:
	Provides the routines to calculate elemental systems, assembling them and applying boundary conditions.
@version: 1.1
'''

from numpy import array, zeros, linalg, matrix
from scipy import sparse
from math import sqrt

global NEN, NEN_range, functions, a, V1, V2, c, f, shape_funcs

def get_element_coords(problem_data, e_nodes):
	"""
	Elemental coordinates are needed in both elemental and global system calculations.
	"""
	def get_node_coords_(node_no):
		return array(problem_data["nodes"][node_no])
	return array(map(get_node_coords_, e_nodes))

def calc_elem(problem_data, e_nodes):
	"""
	Elemental system calculation function to be used in global system calculations.
	"""


	#Initializing elemental values
	Fe = zeros((NEN, 1))
	Ke = zeros((NEN, NEN))
	Se = [0] * NEN
	DKe = [0] * NEN
	DEe = [0] * NEN

	e_coord = get_element_coords(problem_data, e_nodes)
	for GQ_info in problem_data["GQ"]:
		ksi = GQ_info["coord"][0]
		eta = GQ_info["coord"][1]
		for k, shape_func in enumerate(shape_funcs):
			Se[k] = shape_func["main"](ksi, eta)
			DKe[k] = shape_func["dKsi"](ksi, eta)
			DEe[k] = shape_func["dEta"](ksi, eta)

		DS = matrix((DKe, DEe))
		Jacob = DS * matrix(e_coord)
		invJacob = linalg.inv(Jacob)
		detJ = linalg.det(Jacob)
		gDS = invJacob * DS

		#Global coordinate calculation
		x, y = 0, 0
		for i, coord in enumerate(e_coord):
			x += Se[i] * coord[0]
			y += Se[i] * coord[1]

		#Main loop for elemental K calculation
		for i in NEN_range:
			weight = GQ_info["weight"]
			Fe[i] = Fe[i] + Se[i] * f(x, y) * detJ * weight
			for j in NEN_range:
				Ke[i][j] += (a(x, y) * (gDS[0, i] * gDS[0, j] + gDS[1, i] * gDS[1, j]) + Se[i] * (V1(x, y) * gDS[0, j] + V2(x, y) * gDS[1, j]) + c(x, y) * Se[i] * Se[j]) * detJ * weight

	return Ke, Fe


def calc_global(problem_data):
	"""
	Calculates global stiffness matrix, assembly of elemental systems are included here
	instead of defining an extra function for assembly
	"""
	print("Calculating global system...")
	
	global NEN, NEN_range, functions, a, V1, V2, c, f, shape_funcs
	
	#Defining global variables
	NEN = problem_data["NEN"]
	NEN_range = range(NEN)

	#Taking coefficient functions of DE out of problem data
	functions = problem_data["functions"]
	a = functions["a"]
	V1 = functions["V1"]
	V2 = functions["V2"]
	c = functions["c"]
	f = functions["f"]
	
	#Defining shape functions
	shape_funcs = problem_data["shapefunc"]

	print(" * Creating matrixes...")
	NN = problem_data["NN"]
	K = sparse.lil_matrix((NN, NN))
	F = zeros((NN, 1))

	print(" * Calculating K and F matrixes...")
	for e_nodes in problem_data["LtoG"]:
		Ke, Fe = calc_elem(problem_data, e_nodes)
		for i, node_i in enumerate(e_nodes):
			F[node_i] += Fe[i]
			for j, node_j in enumerate(e_nodes):
				K[node_i, node_j] += Ke[i][j]

	print(" * Freeing up memory (1/2)...")
	del problem_data["GQ"]
	del problem_data["UV"]
	del problem_data["functions"]

	K, F = apply_bc(problem_data, K, F)
	print(" * Freeing up memory (2/2)...")
	del problem_data["LtoG"]
	del problem_data["BCs"]

	print(" * Converting LIL to CSR format...")
	K = K.tocsr()
	return K, F


def apply_bc(problem_data, K, F):
	""" 
	Applies all boundary conditions, according to input
	"""
	print(" * Applying boundary conditions...")

	print("  * Applying EBCs...")
	for BC in problem_data["BCs"]["EBC"]:
		node = BC["node"]
		data = BC["data"][0]
		F[node] = data
		K[node, :] = 0.0
		K[node, node] = 1.0

	print("  * Applying NBCs...")
	NEN = problem_data["NEN"]
	for BC in problem_data["BCs"]["NBC"]:
		node1 = BC["face"]
		node2 = (node1 + 1) % NEN
		SV = BC["data"][0]

		e_nodes = problem_data["LtoG"][BC["element"]]
		e_coord = get_element_coords(problem_data, e_nodes)
		length = sqrt(((e_coord[node1] - e_coord[node2]) ** 2).sum())

		SV *= .5 * length
		F[e_nodes[node1]] += SV
		F[e_nodes[node2]] += SV

	print("  * Applying MBCs...")
	for BC in problem_data["BCs"]["MBC"]:
		element = BC["element"]
		node1 = BC["face"]
		node2 = (node1 + 1) % NEN
		alpha = BC["data"][0]
		beta = BC["data"][1]

		e_nodes = problem_data["LtoG"][element]
		e_coord = get_element_coords(problem_data, e_nodes)
		length = sqrt(((e_coord[node1] - e_coord[node2]) ** 2).sum())

		node1 = e_nodes[node1]
		node2 = e_nodes[node2]
		F[node1] += 0.5 * beta * length
		F[node2] += 0.5 * beta * length
		K[node1, node1] -= (alpha * length) / 3.
		K[node2, node2] -= (alpha * length) / 6.

	return K, F
