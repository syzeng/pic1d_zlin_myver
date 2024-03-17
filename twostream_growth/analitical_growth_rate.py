import cmath  

deltax=0.245437
ngrid=64
xsize=ngrid*deltax

km=2
vb=1.5
x=km*vb*2*cmath.pi/(xsize)
print(x)
x2=x**2
y=x2+0.5-0.5*(8*x2+1)**0.5
print(y)
gamma=cmath.sqrt(y).imag
print(gamma)