"""
@author: Deli, BYK
@contact: gulen.ilker@hotmail.com, madbyk@gmail.com
@summary:
    Provides the routines to read and process problem data to be used in the
    calculation of the global system. Also includes Gauss Quadrature and shape
    function definitions. Calling get_problem_data() from an external module is
    enough the get the problem data easily.
@version: 1.2
"""

from math import sqrt
from json import load as json_load
import os

GQ = {
    #----------------------------------------------------
    # Definition of Gauss Quadrature Points
    # Format for use: GQ[eType][NGP]["coords" or "weight"]
    #----------------------------------------------------

    "quad":
    {
        1:
        [
            {"coord": (0., 0.), "weight": 4.}
        ],
        4:
        [
            {"coord": (-sqrt(1. / 3.), -sqrt(1. / 3.)), "weight": 1.},
            {"coord": (sqrt(1. / 3.), -sqrt(1. / 3.)), "weight": 1.},
            {"coord": (-sqrt(1. / 3.), sqrt(1. / 3.)), "weight": 1.},
            {"coord": (sqrt(1. / 3.), sqrt(1. / 3.)), "weight": 1.}
        ],
        9:
        [
            {"coord": (-sqrt(3. / 5.), -sqrt(3. / 5.)), "weight": 25. / 81.},
            {"coord": (0., -sqrt(3. / 5.)), "weight": 40. / 81.},
            {"coord": (sqrt(3. / 5.), -sqrt(3. / 5.)), "weight": 25. / 81.},
            {"coord": (-sqrt(3. / 5.), 0.), "weight": 40. / 81.},
            {"coord": (0., 0.), "weight": 64. / 81.},
            {"coord": (sqrt(3. / 5.), 0.), "weight": 40. / 81.},
            {"coord": (-sqrt(3. / 5.), sqrt(3. / 5.)), "weight": 25. / 81.},
            {"coord": (0., sqrt(3. / 5.)), "weight": 40. / 81.},
            {"coord": (sqrt(3. / 5.), sqrt(3. / 5.)), "weight": 25. / 81.}
        ],
        16:
        [
            {"coord": (-0.8611363116, -0.8611363116), "weight": 0.3478548451 * 0.3478548451},
            {"coord": (-0.3399810435, -0.8611363116), "weight": 0.3478548451 * 0.6521451548},
            {"coord": (0.3399810435, -0.8611363116), "weight": 0.3478548451 * 0.6521451548},
            {"coord": (0.8611363116, -0.8611363116), "weight": 0.3478548451 * 0.3478548451},
            {"coord": (-0.8611363116, -0.3399810435), "weight": 0.6521451548 * 0.3478548451},
            {"coord": (-0.3399810435, -0.3399810435), "weight": 0.6521451548 * 0.6521451548},
            {"coord": (0.3399810435, -0.3399810435), "weight": 0.6521451548 * 0.6521451548},
            {"coord": (0.8611363116, -0.3399810435), "weight": 0.6521451548 * 0.3478548451},
            {"coord": (-0.8611363116, 0.3399810435), "weight": 0.6521451548 * 0.3478548451},
            {"coord": (-0.3399810435, 0.3399810435), "weight": 0.6521451548 * 0.6521451548},
            {"coord": (0.3399810435, 0.3399810435), "weight": 0.6521451548 * 0.6521451548},
            {"coord": (0.8611363116, 0.3399810435), "weight": 0.6521451548 * 0.3478548451},
            {"coord": (-0.8611363116, 0.8611363116), "weight": 0.3478548451 * 0.3478548451},
            {"coord": (-0.3399810435, 0.8611363116), "weight": 0.3478548451 * 0.6521451548},
            {"coord": (0.3399810435, 0.8611363116), "weight": 0.3478548451 * 0.6521451548},
            {"coord": (0.8611363116, 0.8611363116), "weight": 0.3478548451 * 0.3478548451}
        ]
    },
    "tri":
    {
        1:
        [
            {"coord": (1. / 3., 1. / 3.), "weight": 0.5}
        ],
        3:
        [
            {"coord": (0.5, 0.), "weight": 1. / 6.},
            {"coord": (0., 0.5), "weight": 1. / 6.},
            {"coord": (0.5, 0.5), "weight": 1. / 6.}
        ],
        4:
        [
            {"coord": (1. / 3., 1. / 3.), "weight":-27. / 96.},
            {"coord": (0.6, 0.2), "weight": 25. / 96.},
            {"coord": (0.2, 0.6), "weight": 25. / 96.},
            {"coord": (0.2, 0.2), "weight": 25. / 96.}
        ],
        7:
        [
            {"coord": (1. / 3., 1. / 3.), "weight": 0.225 / 2.},
            {"coord": (0.059715871789770, 0.470142064105115), "weight":	0.132394152788 / 2.},
            {"coord": (0.470142064105115, 0.059715871789770), "weight":	0.132394152788 / 2.},
            {"coord": (0.470142064105115, 0.470142064105115), "weight":	0.132394152788 / 2.},
            {"coord": (0.101286507323456, 0.797426985353087), "weight":	0.125939180544 / 2.},
            {"coord": (0.101286507323456, 0.101286507323456), "weight":	0.125939180544 / 2.},
            {"coord": (0.797426985353087, 0.101286507323456), "weight":	0.125939180544 / 2.}
        ]
    }
}

Shape = {
    #==========================================================================
    # Definition of Shape functions and their derivatives
    # Format: Shape[eType]["linear" or "quadratic"][main/dKsi/dEta](ksi, eta)
    #==========================================================================

    "tri":
    {
        "linear":
        [
            {"main": lambda ksi, eta: 1 - ksi - eta, "dKsi": lambda ksi, eta:-1, "dEta": lambda ksi, eta:-1},
            {"main": lambda ksi, eta: ksi, "dKsi": lambda ksi, eta: 1, "dEta": lambda ksi, eta: 0},
            {"main": lambda ksi, eta: eta, "dKsi": lambda ksi, eta: 0, "dEta": lambda ksi, eta: 1}
        ],
        "quadratic":
        [
            {"main": lambda ksi, eta: 2 * (1 - ksi - eta) * (.5 - ksi - eta), "dKsi": lambda ksi, eta :-3 + 4 * eta + 2 * ksi, "dEta": lambda ksi, eta:-3 + 4 * ksi + 2 * eta},
            {"main": lambda ksi, eta: 2 * ksi * (ksi - .5), "dKsi": lambda ksi, eta : 4 * ksi - 1, "dEta": lambda ksi, eta:0},
            {"main": lambda ksi, eta: 2 * eta * (eta - .5), "dKsi": lambda ksi, eta :0, "dEta": lambda ksi, eta: 4 * ksi - 1},
            {"main": lambda ksi, eta: 4 * (1 - ksi - eta) * ksi, "dKsi": lambda ksi, eta : 4 * (1 - 2 * ksi - eta), "dEta": lambda ksi, eta:-4 * ksi},
            {"main": lambda ksi, eta: 4 * ksi * eta, "dKsi": lambda ksi, eta :4 * eta, "dEta": lambda ksi, eta:4 * ksi},
            {"main": lambda ksi, eta: 4 * (1 - ksi - eta) * eta, "dKsi": lambda ksi, eta :-4 * eta, "dEta": lambda ksi, eta:4 * (1 - 2 * eta - ksi)}
        ]
    },
    "quad":
    {
        "linear":
        [
            {"main": lambda ksi, eta: .25 * (1 - ksi) * (1 - eta), "dKsi": lambda ksi, eta:-.25 * (1 - eta), "dEta": lambda ksi, eta:-.25 * (1 - ksi)},
            {"main": lambda ksi, eta: .25 * (1 + ksi) * (1 - eta), "dKsi": lambda ksi, eta: .25 * (1 - eta), "dEta": lambda ksi, eta:-.25 * (1 + ksi)},
            {"main": lambda ksi, eta: .25 * (1 + ksi) * (1 + eta), "dKsi": lambda ksi, eta: .25 * (1 + eta), "dEta": lambda ksi, eta: .25 * (1 + ksi)},
            {"main": lambda ksi, eta: .25 * (1 - ksi) * (1 + eta), "dKsi": lambda ksi, eta:-.25 * (1 + eta), "dEta": lambda ksi, eta: .25 * (1 - ksi)}
        ],
        "quadratic":
        [
            {"main": lambda ksi, eta: .25 * (ksi ** 2 - ksi) * (eta ** 2 - eta), "dKsi": lambda ksi, eta: .25 * (2 * ksi - 1) * (eta ** 2 - eta), "dEta": lambda ksi, eta: .25 * (ksi ** 2 - ksi) * (2 * eta - 1)},
            {"main": lambda ksi, eta: .25 * (ksi ** 2 + ksi) * (eta ** 2 - eta), "dKsi": lambda ksi, eta: .25 * (2 * ksi + 1) * (eta ** 2 - eta), "dEta": lambda ksi, eta: .25 * (ksi ** 2 + ksi) * (2 * eta - 1)},
            {"main": lambda ksi, eta: .25 * (ksi ** 2 + ksi) * (eta ** 2 + eta), "dKsi": lambda ksi, eta: .25 * (2 * ksi + 1) * (eta ** 2 + eta), "dEta": lambda ksi, eta: .25 * (ksi ** 2 + ksi) * (2 * eta + 1)},
            {"main": lambda ksi, eta: .25 * (ksi ** 2 - ksi) * (eta ** 2 + eta), "dKsi": lambda ksi, eta: .25 * (2 * ksi - 1) * (eta ** 2 + eta), "dEta": lambda ksi, eta: .25 * (ksi ** 2 - ksi) * (2 * eta + 1)},
            {"main": lambda ksi, eta: .5 * (1 - ksi ** 2) * (eta ** 2 - eta), "dKsi": lambda ksi, eta: .5 * -2 * ksi * (eta ** 2 - eta), "dEta": lambda ksi, eta: .5 * (1 - ksi ** 2) * (2 * eta - 1)},
            {"main": lambda ksi, eta: .5 * (ksi ** 2 + ksi) * (1 - eta ** 2), "dKsi": lambda ksi, eta: .5 * (2 * ksi + 1) * (1 - eta ** 2), "dEta": lambda ksi, eta: .5 * (ksi ** 2 + ksi) * -2 * eta},
            {"main": lambda ksi, eta: .5 * (1 - ksi ** 2) * (eta ** 2 + eta), "dKsi": lambda ksi, eta: .5 * -2 * ksi * (eta ** 2 + eta), "dEta": lambda ksi, eta: .5 * (1 - ksi ** 2) * (2 * eta + 1)},
            {"main": lambda ksi, eta: .5 * (ksi ** 2 - ksi) * (1 - eta ** 2), "dKsi": lambda ksi, eta: .5 * (2 * ksi - 1) * (1 - eta ** 2), "dEta": lambda ksi, eta: .5 * (ksi ** 2 - ksi) * -2 * eta},
            {"main": lambda ksi, eta: (1 - ksi ** 2) * (1 - eta ** 2), "dKsi": lambda ksi, eta:-2 * ksi * (1 - eta ** 2), "dEta": lambda ksi, eta: (1 - ksi ** 2) * -2 * eta}
        ]
    }
}

Order = {
    "tri": {3: "linear", 6: "quadratic"},
    "quad": {4: "linear", 9: "quadratic"}
}


def process_functions(functions, UV_data, nodes):
    """
    Processes coefficient functions of the DE/problem to create directly
    callable functions from Python.
    """

    default_lambda = "lambda x,y:"
    global_vars = None
    for name in functions:
        if functions[name] == '?':
            functions[name] = '0'
        elif functions[name] == "x" or functions[name] == "y":
            #If it is indicated that the provided U & V values to be used
            if not global_vars:
                x, y = [0] * nodes.__len__(), [0] * nodes.__len__()
                for i, node in enumerate(nodes):
                    x[i] = node[0]
                    y[i] = node[1]

                from scipy.interpolate import bisplrep, bisplev
                # Fit a bivariate B-spline to U and V values to calculate
                # values that are not on the nodes  This "global_vars"
                # dictionary is provided to eval for the lambda's to work
                # properly
                global_vars = {
                    "x_tck": bisplrep(x, y, UV_data[0]),
                    "y_tck": bisplrep(x, y, UV_data[1]),
                    "bisplev": bisplev
                }

            functions[name] = eval(
                "lambda x,y: bisplev(x, y, {0}_tck)".format(functions[name]),
                global_vars
            )
            continue

        functions[name] = default_lambda + functions[name]
        functions[name] = eval(functions[name])
    return functions


def process_problem_data(problem_data):
    """
    Takes the raw problem data then converts the string functions into usable
    functions with process_functions, determines necessary shape funcstions and
    embeds them to problem_data with necessary GQ info.
    """

    eType = problem_data["eType"]
    eOrder = Order[eType][problem_data["NEN"]]
    problem_data["GQ"] = GQ[eType][problem_data["NGP"]]
    problem_data["shapefunc"] = Shape[eType][eOrder]

    if not "UV" in problem_data:
        problem_data["UV"] = None

    if not "title" in problem_data:
        problem_data["title"] = "Untitled Problem"

    process_functions(
        problem_data["functions"],
        problem_data["UV"],
        problem_data["nodes"]
    )

    return problem_data


def read_problem_data(input_name='', output_name=''):
    """
    Reads the problem data from the user provided file name.
    The file can either be an .inp file or a .json file.
    """

    if not input_name:
        input_name = raw_input('Enter input file name: ')

    file_name, file_ext = os.path.splitext(input_name)

    if file_ext == "":
        if os.path.exists(file_name + ".json"):
            file_ext = ".json"
        else:
            print("Cannot find valid input file. Expecting a .json file.")
            exit(127)

    with open(file_name + file_ext, "r") as input_file:
        problem_data = json_load(input_file)

    if output_name:
        problem_data["output"] = output_name

    if ("output" not in problem_data) or (not problem_data["output"]):
        problem_data["output"] = file_name + "_output.json"

    return problem_data


def get_problem_data(input_name='', output_name=''):
    return process_problem_data(read_problem_data(input_name, output_name))
