import numpy as np

# OLS regression functions
def ols_slope(x, y):
    """
    Calculates the slope of the linear regression line.
    Formula: slope = Σ((x - mean_x) * (y - mean_y)) / Σ((x - mean_x)^2)
    x: independent variable (e.g., T_B).
    y: dependent variable (e.g., H_v).
    Returns the slope of the regression line.
    """
    x_mean = np.mean(x)
    y_mean = np.mean(y)
    numerator = np.sum((x - x_mean) * (y - y_mean))
    denominator = np.sum((x - x_mean) ** 2)
    return numerator / denominator

def ols_intercept(x, y):
    """
    Calculates the intercept of the linear regression line.
    Formula: intercept = mean_y - (slope * mean_x)
    x: independent variable.
    y: dependent variable.
    Returns the intercept of the regression line.
    """
    x_mean = np.mean(x)
    y_mean = np.mean(y)
    slope = ols_slope(x, y)
    return y_mean - slope * x_mean

def ols(x, y):
    """
    Combines both slope and intercept calculations.
    x: independent variable.
    y: dependent variable.
    Returns a tuple (slope, intercept).
    """
    slope = ols_slope(x, y)
    intercept = ols_intercept(x, y)
    return slope, intercept