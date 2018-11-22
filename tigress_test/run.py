#!/usr/bin/env python

import numpy as np
import os
import subprocess
import time

#read KLM.in file

input_file = "cooling.in"
dir_out = os.getcwd() + "/"

#read list of problems in input file
fin = open(input_file, "r")
fprob = []
exe = []
dout_prob = []
for l in fin.readlines():
    if l[0] != "#" and (l.strip()):
        lsp = l.split()
        fprob.append(lsp[0])
        exe.append(lsp[1] + ".exe")
        dout_prob.append(dir_out + lsp[2])
fin.close()

nprob = len(fprob)

#loop over each problem, do
#(1) mkdir/cleanup output directory
#(2) rewrite makefile
#(3) call make to compile
#(4) run in background
for i in xrange(nprob):
    #(1) mkdir/cleanup output directory
    if os.path.isdir(dout_prob[i]):
        command_clean = "rm {}/*".format(dout_prob[i])
        print command_clean
        proc1 = subprocess.Popen(command_clean, shell=True).wait()
        #pid1, sts1 = os.waitpid(proc1.pid, 0)
    else:
        proc1 = subprocess.Popen("mkdir {}".format(dout_prob[i]),
                shell=True).wait()
        #pid1, sts1 = os.waitpid(proc1.pid, 0)

    #(2) rewrite makefile
    fmakein = open("Makefile.in", "r")
    lines_makein = fmakein.readlines()
    fmakein.close()
    fmake = open ("Makefile", "w")
    fmake.write("prob = {}\n".format(fprob[i]))
    fmake.write("exe = {}\n".format(exe[i]))
    for l_makein in lines_makein:
        fmake.write(l_makein)
    fmake.close()

    #(3) call make to compile
    proc3 = subprocess.Popen("make clean", shell=True).wait()
    proc3 = subprocess.Popen("make", shell=True).wait()
    #pid3, sts3 = os.waitpid(proc3.pid, 0)

    #(4) run in background
    run_command = "./{}".format(exe[i])
    print run_command
    proc4 = subprocess.Popen(run_command)
    proc4.wait()
    #pid4, sts4 = os.waitpid(proc4.pid, 0)
