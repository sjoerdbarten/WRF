#!/usr/bin/env python
# MyPythonToolbox.py

"""
Author : van der Molen 

Revision History:
File created on 14 Feb 2012.
"""

def msg(string,type=2):
    import sys

    if type == 1: print(string),
    if type == 2: print(string)
    
    sys.stdout.flush()
    return    