import sys

sys.path.insert(1, '../src/python/lib/')

import pytest
import pandas as pd
import math
from sum_of_squares import calc_sse
from fitness_score import calc_fitness

def test_calc_sse():

    assert calc_sse(pd.DataFrame([[1, 2, 3], [4, 5, 6], [7, 8, 9]]), pd.DataFrame([[9, 8, 7], [6, 5, 4], [3, 2, 1]]), 3) == 240
    assert calc_sse(pd.DataFrame([[10, 20, 30], [40, 50, 60], [70, 80, 90]]), pd.DataFrame([[90, 80, 70], [60, 50, 40], [30, 20, 10]]), 3) == 24000
    assert calc_sse(pd.DataFrame([[10, 20, 30], [40, 50, 60], [70, 80, 90]]), pd.DataFrame([[90, 80, 70], [60, 50, 40], [30, 20, 10]]), 3) != 20
    assert calc_sse(pd.DataFrame([[35, 99, 106], [32, 76, 83], [22, 5, 68]]), pd.DataFrame([[49, 53, 110], [12, 88, 91], [70, 80, 55]]), 3) == 11034

def test_calc_fitness():

    assert calc_fitness(10000, 30000, 1.0, 0.001) == 1.0
    assert calc_fitness(50000, 10000, 1.0, 0.001) != 1.0
    assert calc_fitness(30000, 5000, 1.0, 0.001) == 1.9548290415713273e-22
    assert calc_fitness(60000, 45000, 1.0, 0.001) != 1e-8
