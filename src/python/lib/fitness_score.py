import math
import sys

#Calculates the fitness of the new mutation
def calc_fitness(variant_fit, orig_fit, N, beta_val):
    """
    Calculates the fitness of the new mutation and compares it to the fitness of the old mutation.
    The 'orig_fit' is the old fitness value previously calculated.
    The 'variant_fit' is the new fitness value calculated to determine if mutation is accepted or rejected.
    The 'beta' value is determined in relationship to scale of the sum of squares values.
    """

    Ne = N
    beta = beta_val
    thresholds = 0


    xi = calc_x(orig_fit, beta, thresholds)
    xj = calc_x(variant_fit, beta, thresholds)


    if xj >= xi:
        return((1.0))
    else:
        exponent = -2 * float(Ne) * (xi - xj)
        return(safe_calc(exponent))


def calc_x(data, beta, threshold):
    total = 0
    exponent = float(beta) * (float(data) - float(threshold))
    total += -math.log(safe_calc(exponent) + 1)
    return(total)


def safe_calc(exponent):
    if exponent > 700:
        print("system maxed")
        return(sys.float_info.max)
    else:
        return(math.exp(exponent))
