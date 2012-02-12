from __future__ import division
#Copyright 2011 Thomas A Caswell
#tcaswell@uchicago.edu
#http://jfi.uchicago.edu/~tcaswell
#
#This program is free software; you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation; either version 3 of the License, or (at
#your option) any later version.
#
#This program is distributed in the hope that it will be useful, but
#WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
#General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program; if not, see <http://www.gnu.org/licenses>.

#from 
# http://www.tandfonline.com/doi/abs/10.1080/00268977000101421
# 1 W. Smith and D. Henderson, Molecular Physics 19, 411-415 (1970).

#

def S(t,eta):
	#equation 3
	return (1-eta)**2 * t**3 + 6*eta*(1-eta) * t**2 + 18*(eta**2)* t - 12*eta * (1+2*eta)

def S1(t,eta):
	# derivative of (3)
	return 3*(1-eta)**2 * t**2 + 12*eta*(1-eta) * t + 18*(eta**2)
def S2(t,eta):
	# second derivative of (3)
	return 6*(1-eta)**2 * t + 12*eta*(1-eta)

def L(t,eta):
	#equation (4)
	return (1+0.5*eta)*t + 1 + 2 * eta

def L1(t,eta):
	# derivative of (4)
	return (1+0.5*eta)

def A1(t,eta):
	# equation (9)
	return (L(t,eta)/S1(t,eta))**2

def A2(t,eta):
	# equation (9)
	return 15*S2(t,eta) - 4*S1(t,eta)*S3(t,eta)
def A3(t,eta):
	# equation (9)
	return	L(t,eta) + 4*t*L1(t,eta)
def A4(t,eta):
	# equation (9)
	return L(t,eta)*S2(t,eta)/S1(t,eta)
def A5(t,eta):
	# equation (9)
	return	2*L(t,eta) + 3*t*L1(t,eta)
	
def ti(i,eta):
    # equation (5)
    J = np.exp(np.pi*2j/3)
    return (-2*eta + ((2*eta*f(eta))**(1/3))*(y_p(eta)*(J**i) + y_m(eta)*(J**-i)))/(1-eta)


def y_p(eta):
    # equation (6)
    return  (1+np.sqrt(1+2*( (eta**2) / f(eta) )**2) )**(1/3)
    
def y_m(eta):
	# equation (6)
	tmp = (1-np.sqrt(1+2*( (eta**2) / f(eta) )**2) )
	sg = np.sign(tmp)
	return sg*(np.abs(tmp)**(1/3))

def f(eta):
	# equation (6)
	return 3 +3*eta - eta**2


# below all from (8)
def a0(i,eta):
	_ti = ti(i,eta)
	return	_ti*L(_ti,eta)/S1(_ti,eta)

def b0(i,eta):
	_ti = ti(i,eta)
	return -12*eta*L(_ti,eta)/(S1(_ti,eta)**2)

def b1(i,eta):
	_ti = ti(i,eta)
	return (1-_ti*S2(_ti,eta)/S1(_ti,eta))*L(_ti,eta) + 2*_ti*L1(_ti,eta)

def b2(i,eta):
	_ti = ti(i,eta)
	return _ti*L(_ti,eta)

def c0(i,eta):
	_ti = ti(i,eta)
	return	72*(eta**2)*L(_ti,eta)/(S1(_ti,eta)**3)

def c1(i,eta):
	_ti = ti(i,eta)
	return A1(_ti,eta)*t * (3*S2(_ti,eta)**2 - S1(_ti,eta)*S3(_ti,eta)) - 3*L(_ti,eta)*S2(_ti,eta)/(S1(_ti,eta)*(L(_ti,eta) + 3*_ti*L1(_ti,eta))) + 6*L1(_ti,eta)*(L(_ti,eta) + _ti*L1(_ti,eta))

def c2(i,eta):
	_ti = ti(i,eta)
	return (6*_ti*L1(_ti,eta) + L(_ti,eta)*(2-3*_ti*S2(_ti,eta)/S1(_ti,eta)))*L(_ti,eta)

def c3(i,eta):
	_ti = ti(i,eta)
	return _ti*L(_ti,eta)**2

def d0(i,eta):
	_ti = ti(i,eta)
	return -288*(eta**3)*L(_ti,eta)/(S1(_ti,eta)**4)

def d1(i,eta):
	_ti = ti(i,eta)
	return	5*_ti*((L(_ti,eta)/S1(_ti,eta))**3) * S2(_ti,eta) * (2*S1(_ti,eta)*S3(_ti,eta) - 3*S2(_ti,eta)**2) + A1(_ti,eta)*A2(_ti,eta)*A3(_ti,eta) - 24*L1(_ti,eta)*A4(_ti,eta)*A5(_ti,eta) + 12*(L1(_ti,eta)**2)*(3*L(_ti,eta) + 2*_ti*L1(_ti,eta))

def d2(i,eta):
	_ti = ti(i,eta)
	return (_ti*A1(_ti,eta)*A2(_ti,eta) - 12*(A3(_ti,eta)*A4(_ti,eta) - L1(_ti,eta)*A5(_ti,eta)))*L(_ti,eta)

def d3(i,eta):
	_ti = ti(i,eta)
	return (3*A3(_ti,eta) - 6*_ti*A4(_ti,eta))*L(_ti,eta)**2

def d4(i,eta):
	_ti = ti(i,eta)
	return _ti*L(_ti,eta)**3


def g1_part(i,eta,x):
	return a0(i,eta) *np.exp(ti(i,eta) *(x-1))

def g1(eta,x):
	return (g1_part(0,eta,x) + g1_part(1,eta,x)+ g1_part(2,eta,x))/x

def g2_part(i,eta,x):
	return b0(i,eta) * np.exp(ti(i,eta)*(x-2))*(b1(i,eta) + b2(i,eta)*(x-2))

def g2(eta,x):
	return (g2_part(0,eta,x) + g2_part(1,eta,x)+ g2_part(2,eta,x))/x
