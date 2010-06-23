'''
@author: BYK
@contact: madbyk@gmail.com
@summary:
	Provides the routines to read an ".inp" file whose format is defined by Dr. Sert into
	an easy to use data structure in the program. When run on its own, it converts ".inp"
	files into ".json" files.
@version: 1.2
'''
import re
from json import dump as json_dump
from os.path import splitext as os_path_splitext

def imap_dict(func, dict_obj):
	"""
	Applies "func" to all of the values of dict and replaces the original value with the result of "func".
	"""
	for key in dict_obj:
		dict_obj[key] = func(dict_obj[key])
	return dict_obj

def locate_value_line(file_object, info_pattern):
	"""
	Reads lines from file_object until it reaches a line which matches the RegEx pattern given in "info_pattern".
	Used to position file cursor to the correct place for reading in arbitrary order.
	"""
	info_matcher = re.compile(info_pattern)
	info_line = ' ';
	while not info_matcher.search(info_line) and info_line != "":
		info_line = file_object.readline()

def read_fundamental_variables(file_object):
	"""
	Searches for "eType NE NN NEN NGP" line to the end of file from current position
	then constructs the dictionary object which contains the values as integers under the same names.
	"""
	locate_value_line(file_object, r"^eType\s+NE\s+NN\s+NEN\s+NGP")
	eType_TRANSLATIONS = {1: "quad", 2: "tri"}
	fund_vars = imap_dict(int, re.search(r"(?P<eType>[12])\s+(?P<NE>\d+)\s+(?P<NN>\d+)\s+(?P<NEN>\d+)\s+(?P<NGP>\d+)", file_object.readline()).groupdict())
	fund_vars['eType'] = eType_TRANSLATIONS[fund_vars['eType']]
	return fund_vars

def read_problem_functions(file_object):
	"""
	Searches for "a V1 V2 c f exactSoln" line to the end of file from current position
	then constructs the dictionary object which contains the values as strings under the same names.
	"""
	locate_value_line(file_object, r"^a\s+V1\s+V2\s+c\s+f\s+exactSoln")
	return {
					"a": file_object.readline().strip(),
					"V1": file_object.readline().strip(),
					"V2": file_object.readline().strip(),
					"c": file_object.readline().strip(),
					"f": file_object.readline().strip(),
					"exactSoln": file_object.readline().strip()
				 }

def read_mesh(file_object, node_count):
	"""
	Searches for "Node# x y" or "Node No x y" line to the end of the file from current position
	then reads node_count lines from the file and constructs the nodes array.
	The nodes can be given in any order in the file provided that they have node numbers in the first column.
	"""
	value_matcher = re.compile(r"(?P<node>\d+)\s+(?P<x>\S+)\s+(?P<y>\S+)")
	locate_value_line(file_object, r"^(Node#|Node No)\s+x\s+y")
	nodes = [[]] * node_count
	for i in range(node_count):
		values = value_matcher.search(file_object.readline()).groupdict()
		nodes[int(values['node']) - 1] = ([float(values['x']), float(values['y'])])
	return nodes

def one_less_int_(val):
	return int(val) - 1

def read_LtoG(file_object, element_count):
	"""
	Searches for the line starting with "Elem# node1 node2 node3" or "Elem No node1 node2 node3" to
	the end of the file from current position then reads element_count lines from the file and
	constructs the LtoG matrix. Elements can be given in any order in the file provided that
	they have element numbers in the first column.
	"""
	value_matcher = re.compile(r"\d+")
	locate_value_line(file_object, r"^(Elem#|Elem No)\s+node1\s+node2\s+node3")
	elements = [[]] * element_count
	for i in range(element_count):
		values = map(one_less_int_, value_matcher.findall(file_object.readline()))
		elements[values[0]] = values[1:]
	return elements

def read_boundary_conditions(file_object):
	"""
	Searches for the line starting with "nBCdata" to determine the BC count then
	searches for the line starting with "nEBCnodes nNBCfaces nMBCfaces" to determine
	the nodes and faces which are subject to provided BCs and puts this information
	into problem_data for a more easy-to-use manner.
	"""
	value_matcher = re.compile(r"\S+")
	locate_value_line(file_object, r"^nBCdata")
	BC_data_count = int(file_object.readline())
	BC_datas = []

	for i in range(BC_data_count):
		values = value_matcher.findall(file_object.readline())
		BC_datas.append(map(float, values[1:]))

	locate_value_line(file_object, r"^nEBCnodes\s+nNBCfaces\s+nMBCfaces")
	BC_applicant_count = imap_dict(int, re.search(r"(?P<EBC>\d+)\s+(?P<NBC>\d+)\s+(?P<MBC>\d+)", file_object.readline()).groupdict())
	BCs = {"EBC": [], "NBC": [], "MBC": []}

	locate_value_line(file_object, r"^EBC Data\s+\(Node\s+BCno\)")
	for i in range(BC_applicant_count['EBC']):
		values = map(one_less_int_, value_matcher.findall(file_object.readline()))
		BCs["EBC"].append({"node": values[0], "data": BC_datas[values[-1]]})

	locate_value_line(file_object, r"^NBC Data\s+\(Elem\s+Face\s+BCno\)")
	for i in range(BC_applicant_count['NBC']):
		values = map(one_less_int_, value_matcher.findall(file_object.readline()))
		BCs["NBC"].append({"element": values[0], "face": values[1], "data": BC_datas[values[-1]]})

	locate_value_line(file_object, r"^MBC Data\s+\(Elem\s+Face\s+BCno\)")
	for i in range(BC_applicant_count['MBC']):
		values = map(one_less_int_, value_matcher.findall(file_object.readline()))
		BCs["MBC"].append({"element": values[0], "face": values[1], "data": BC_datas[values[-1]]})

	return BCs

def read_UV_values(file_object, node_count):
	"""
	Searches for the line starting with "Node No U V" to determine then reads
	"node_count" number of lines and parses the values inside them to the UV list.
	"""
	value_matcher = re.compile(r"\S+")
	locate_value_line(file_object, r"^Node No\s+U\s+V")
	UV = [[0] * node_count, [0] * node_count]
	for i in range(node_count):
		values = value_matcher.findall(file_object.readline())
		if values.__len__() == 0:
			return None
		values[1:3] = map(float, values[1:3])
		UV[0][one_less_int_(values[0])] = values[1]
		UV[1][one_less_int_(values[0])] = values[2]
	return UV

def read_input_data(file_object):
	"""
	Reads the file located at file_name using the functions above and constructs the "data"
	dictionary which contains all the problem information in a structured manner.
	"file_name" should be a full file name including the file extension.
	"""
	print "Parsing input data..."
	data = read_fundamental_variables(file_object)
	data.update(
		{
			"functions": read_problem_functions(file_object),
			"nodes": read_mesh(file_object, data['NN']),
			"LtoG": read_LtoG(file_object, data['NE']),
			"BCs": read_boundary_conditions(file_object),
			"UV": read_UV_values(file_object, data['NN'])
		}
	)
	return data

if __name__ == "__main__":
	file_name = raw_input('Enter input file name: ')
	file_name_parts = os_path_splitext(file_name)

	if file_name_parts[1] != '.inp' and file_name_parts[1] != '':
		print "There is nothing I can do with this file, sorry."
		exit()
	elif file_name_parts[1] == '':
		file_name += '.inp'

	input_file = open(file_name, 'r')
	json_file = open(file_name_parts[0] + '.json', 'w')
	data = read_input_data(input_file)
	print "Writing JSON file..."
	json_dump(data, json_file, indent = 2)
	json_file.close()
	input_file.close()
	print "JSON file created successfully."
