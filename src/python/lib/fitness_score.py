import math
import sys

#1e-5

#Calculates the fitness of the new mutation
def calc_fitness(variant_fit, orig_fit, generations, count):
    """
    Calculates the fitness of the new mutation and compares it to the fitness of the old mutation.
    Code adapted from Ashley Teufel and Claus Wilke.
    Input(s):
    variant_fit is the new fitness value calculated to determine if mutation is accepted or rejected.
    orig_fit is the old fitness value previously calculated.
    N is the effective population size that determines how strict selection should be.
    beta is the value determined in relationship to the scale of the sum of squares values.
    Output(s):
    Returns the comparison between fitness values.
    """

    Ne = 1000
    #Determines the value that controls variation within the simulation
    if count <= 0.1 * generations:
        beta = 1e-7
    else:
        #y = mx + b: linearly increases
        slope = (1e-4 - 1e-7) / generations
        beta = (slope * (count - (0.1 * generations))) + 1e-7
    thresholds = 0

    #Fitness values are calculated based on the new and current sum of squared values
    xi = calc_x(orig_fit, beta, thresholds)
    xj = calc_x(variant_fit, beta, thresholds)

    #Fitness values are compared to determine if a mutation should be accepted
    if xj >= xi:
        return 1.0
    #Deleterious mutations are accepted exponentially
    else:
        exponent = -2 * float(Ne) * (xi - xj)
        return safe_calc(exponent)


def calc_x(data, beta, threshold):
    """
    Calculates the fitness values based on the sum of squared error value.
    Code adapted from Ashley Teufel and Claus Wilke.
    Input(s):
    data is the sum of squared value.
    beta is the value determined in relationship to the scale of the sum of squares values.
    threshold is the ...
    Output(s):
    Returns the calculated fitness value for comparison.
    """

    total = 0
    exponent = float(beta) * (float(data) - float(threshold))
    total += -math.log(safe_calc(exponent) + 1)
    return total


def safe_calc(exponent):
    """
    Verify the value is less than 700 as any exponent greater would be too large and unnecessary to calculate.
    Code adapted from Ashley Teufel and Claus Wilke.
    Input(s):
    exponent is the value pertaining to the exponent needed to convert a sum of squared value to a fitness.
    Output(s):
    Returns an exponent that is within a "safe" range that is not too big to handle.
    """

    if exponent > 700:
        return sys.float_info.max
    else:
        return math.exp(exponent)
