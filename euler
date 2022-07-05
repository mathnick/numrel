import numpy as np

#exemplo para y'= 3 - y/x

def euler(f,x0,y0,h,n):

    for k in range(n):
        yk = y0 + h*f(x0,y0)
        xk = x0 + h
        print(yk,xk)
        y0=yk
        x0=xk



def f(x,y):
    return 3-x/y
x0=2
y0=2
h=0.1
n=10

euler(f,x0,y0,h,n)

