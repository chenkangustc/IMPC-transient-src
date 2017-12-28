#!/usr/bin/env python3
# -*- coding:utf-8 -*-
'''
parse src files and generate dependence relationship for makefile
'''

import os, sys

if __name__ == "__main__":
    startpath = "."
    if len(sys.argv) > 1:
        startpath = sys.argv[1]
    
    outpfile = "dependence.inc"
    extfile = "dependence.ext"
    extname = (".f90", ".F90")
    
    # get module names in every src file
    file2mod = dict()
    for root, dir, files in os.walk(startpath):
        for ff in files:
            fullpath = os.path.join(root, ff)
            if ff.endswith(extname):
                inp = open(fullpath, "rt", errors="ignore")
                mods = list()
                for line in inp:
                    if line.isspace():
                        continue
                    if line.strip().startswith("!"):
                        continue
                    if line.strip().lower().split()[0] == "module" and (not "procedure" in line.lower()):
                        k = line.strip().lower().split()[1]
                        mods.append(k)
                fbase, _ = os.path.splitext(ff)
                file2mod[fbase] = mods
                inp.close()
            
    # for every module in src file, get the file name that contains it
    csum = 0
    extlib = dict()
    outp = open(outpfile, "wt")
    for root, dir, files in os.walk(startpath):
        for ff in files:
            fullpath = os.path.join(root, ff)
            if ff.endswith(extname):
                inp = open(fullpath, "rt", errors="ignore")
                csum += 1
                print("----> {0: >3d}: {1:s}".format(csum, ff))
                is_nonnull = False
                mods = list()
                for line in inp:
                    if line.isspace():
                        continue
                    if line.strip().startswith("!"):
                        continue
                    if (line.strip().lower().split()[0] == "use") and (not "intrinsic" in line.lower().split()):
                        mod = line.strip().lower().split()[1]
                        if "," in mod:
                            cloc = mod.index(",")
                            mod = mod[ :cloc]
                        
                        is_find = False
                        for k, v in zip(file2mod.keys(), file2mod.values()):
                            if mod == v or (mod in v):
                                is_find = True
                                if not mod in mods:
                                    is_nonnull = True
                                    mods.append(mod)
                                    src = k
                                    dst, _ = os.path.splitext(ff)
                                    if dst == src:
                                        continue
                                    cnt = "{0}.o: {1}.o".format(dst, src)
                                    outp.write(cnt + "\n")
                                    break
                        if not is_find:
                            extlib[ff] = mod
                if is_nonnull:
                    outp.write("\n")
                inp.close()
    outp.close()
                   
    # if a module is not contained by any other src file, it is an external module
    if len(extlib.keys()) > 0:
        outp = open(extfile, "wt")
        for f, m in extlib.items():
            cnt = "{0: >30s}: {1: >20s}".format(f, m)
            outp.write(cnt + "\n")
        outp.close()
    