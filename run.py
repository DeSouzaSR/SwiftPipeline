#!/usr/bin/env python
#-*- coding:utf-8 -*-
from configobj import ConfigObj
import os, sys
import shutil

def create_file(file):
    simulation_name = file.split(".")[0]
    config = ConfigObj(file)
    suffix_simu = config["suffix_simu"]
    number_clones = int(config["number_clones"])
    
    if os.path.isdir(simulation_name):
        shutil.rmtree(simulation_name)
        os.mkdir(simulation_name)
    else:
        os.mkdir(simulation_name)

    for i in suffix_simu:
        os.makedirs(simulation_name + "/" + simulation_name + i)
        for j in range(number_clones):
            os.makedirs(simulation_name + "/" + simulation_name + i + "/" + \
                        simulation_name + i + "{:03d}".format(j))
  
def main():
    create_file(sys.argv[1])
        
if __name__ == '__main__':
    main()
    
