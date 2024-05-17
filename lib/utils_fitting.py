import os, glob, pickle, pprint, copy
import numpy as np

from skimage.measure import label
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import r2_score
from scipy.optimize import curve_fit
import warnings


def efron_rsquare(y, y_pred):
	n = float(len(y))
	t1 = np.sum((y - y_pred) * (y - y_pred))
	t2 = np.sum((y - np.average(y)) * (y - np.average(y)))
	return 1.0 - (t1 / t2)

def get_count_most_freq_outcome(y):
	num_0 = 0
	num_1 = 0
	for p in y:
		if p == 1.0:
			num_1 += 1
		else:
			num_0 += 1
	return float(max(num_0, num_1))

def count_adjusted_rsquare(y, y_pred, t=0.5):
    correct = get_num_correct(y, y_pred, t)
    total = float(len(y))
    n = get_count_most_freq_outcome(y)
    return (correct - n) / (total - n)


def get_num_correct(y, y_pred, t=0.5):
	y_correct = np.array([0.0 if p < t else 1.0 for p in y_pred])
	return sum([1.0 for p, p_pred in zip(y, y_correct) if p == p_pred])

def count_rsquare(y, y_pred, t=0.5):
	n = float(len(y))
	num_correct = get_num_correct(y, y_pred, t)
	return num_correct / n


def hill(x, a, b, c): # Hill sigmoidal equation from zunzun.com
	return  a * np.power(x, b) / (np.power(c, b) + np.power(x, b))




