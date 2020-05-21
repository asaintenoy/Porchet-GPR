#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 21 17:16:54 2020

@author: el
"""
import numpy as np

from joblib import Parallel, delayed


def ola(a,b):
    return np.sqrt(a+b)


Parallel(n_jobs=2, prefer="threads")(delayed(ola)(i ** 2,18) for i in range(10))