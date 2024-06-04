from os import getcwd

def read_input():
    cwd = getcwd()
    main_infile = open('%s/asc.in'%cwd, 'r')
    line = main_infile.readline()
    global indict
    indict = {}
    while line:
        line = main_infile.readline()
        llist = line.split('=')
        
        if llist != ['\n'] and llist != ['']:
            if llist[0][0] != '#':
                inputlist = [i.strip().split() for i in llist]
                if inputlist[1] == []:
                    with open('../log.asc', 'a') as logfile:
                        print >>logfile, "Please give the value(s) for: %s" % inputlist[0][0]
                        
                else:
                    indict[inputlist[0][0]] = inputlist[1]
                    
    return indict


indict = read_input()
