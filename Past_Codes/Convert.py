import sys

parameters_to_convert = ["ux", "epsWCA", "Dmorse"]
path = "Data/" + sys.argv[1]
with open(path, 'r') as file:
    lines = file.readlines()

for i, line in enumerate(lines):
    line_split = line.split()
    # find the parameter needed to be converted
    for parameter in parameters_to_convert:
        try:
            order_parameter = line_split.index(parameter)
            index_line = i
        except:
            pass
    # convert the value of the parameter
    try:
        if i == (index_line + 1):
            value = 1
            for digit in line_split[order_parameter].split('*'):
                value = value*float(digit)
            line_split[order_parameter] = str(value)
            lines[i] = '   '.join(line_split) + "\n"
    except:
        pass
# write the converted value
with open(path, 'w') as file:
    file.writelines( lines )