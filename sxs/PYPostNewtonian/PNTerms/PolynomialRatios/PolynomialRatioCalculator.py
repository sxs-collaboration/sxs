#! /usr/bin/env ipython
# coding=utf-8
import sys
import pickle
import sympy

ip = get_ipython()
for order in [int(o) for o in range(10)]:
    print("Calculating with polynomials to order {0}".format(order))
    Num = sympy.var('Num:{0}'.format(order), real=True)
    Den = sympy.var('Den:{0}'.format(order), real=True)
    PolynomialVariable = sympy.Symbol('PolynomialVariable', real=True)
    p_Num = sum(Num[i]*PolynomialVariable**i for i in range(order))
    p_Den = sum(Den[i]*PolynomialVariable**i for i in range(order))
    print(p_Num, p_Den, Num, Den)
    try:
        p_Ratio = sympy.series(p_Num/p_Den,x=PolynomialVariable,x0=0,n=order)
    except:
        p_Ratio = 1
    #ip.magic("time p_Ratio = sympy.series(p_Num/p_Den,x=PolynomialVariable,x0=0,n=order)")
    with open('PolynomialRatioSeries_Order{0}.dat'.format(order), 'wb') as fff:
        pickle.dump(p_Ratio, fff)

# Load later, in a different python session, with:
# p_Ratio = pickle.load(file('PolynomialRatioSeries_Order{0}.dat'.format(order)))
