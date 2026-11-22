#!/usr/bin/env python3

# to analyze individual-based sims, but then parameters are given upfront

import os, re, sys

first = True


def analyze_parameters(lines,first=False):

    pars = {}

    for line in lines:
        mobj = re.match("(.*);(.*)",line)
        if mobj != None:
            pars[mobj.group(1)] = mobj.group(2)

    return(pars)

def analyze_data(lines):

    data = [ [ float(celli) ] for celli in lines[0].split(";")[0:-1] ]

    # loop through the lines to collect the data
    for line in lines[1:]:
        splitline = line.split(";")[0:-1]

        for i in range(0,len(splitline)):
            data[i].append(float(splitline[i]))

    # now take averages

    avgs = []
    for i in range(0,len(data)):
        avgs.append(sum(data[i])/len(data[i]))

    return(avgs)

def analyze_file(filename):

    global first;

    # open file; read data
    f = open(filename)
    fl = f.readlines()
    f.close

    endline = 0

    # first see until what line the parameters stretch
    for idx, line_i in enumerate(fl):
        if re.match("^time",line_i) != None:
            endline = idx
            break

    # something went wrong. Could not find 
    # the header. return
    if endline == len(fl):
        return

    flhead = fl[endline]

    parameters = analyze_parameters(fl[0:endline])

    # check whether simulation actually returned output
    # otherwise skip it
    if re.match("^\d.*",fl[-1]) == None:
        return
    
    if first:
        print(";".join(parameters.keys()) + ";" + flhead.strip() + "file")
        first = False

    print(";".join(parameters.values()) + ";" + fl[-1].strip() + filename)


for root, dirs, files in os.walk(sys.argv[1]):
    for name in files:
        if re.match("(sim|iter).*",name) != None:
#            print(os.path.join(root,name))
            data = analyze_file(os.path.join(root,name))


