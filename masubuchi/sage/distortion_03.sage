import math

def bkls_miller(P, Q, n):

  if Q.is_zero():
      raise ValueError("Q must be nonzero.")
  if n.is_zero():
      raise ValueError("n must be nonzero.")
  n_is_negative = False
  if n < 0:
      n = n.abs()
      n_is_negative = True
  one = P.curve().base_field().one()
  t = one ##1が入っている
  V = P
  nbin = n.bits()
  i = n.nbits() - 2
  count = 1
  Q= distortion_map(Q)
  print('count, t', [count, t])
  while i > -1:
      ell = V._line_(V, Q)
      t = (t**2)*ell
      count = 2*count
      print('count, t', [count, t])
      V = 2*V
      if nbin[i] == 1:
          ell = V._line_(P, Q)
          t = t*ell
          count += 1
          print('count, t', [count, t])
          V = V+P
      i = i-1
  print('count, t', [count, t])
  if n_is_negative:
      vee = V._line_(-V, Q)
      t = 1/(t*vee)
  return t

def window_bkls_miller(P, Q, n, Pi, fi, w):

  if Q.is_zero():
      raise ValueError("Q must be nonzero.")
  if n.is_zero():
      raise ValueError("n must be nonzero.")
  n_is_negative = False
  if n < 0:
      n = n.abs()
      n_is_negative = True
  one = P.curve().base_field().one()
  t = one ##1が入っている
  V = P
  nbin = n.bits()
  i = n.nbits() - 2
  count = 1
  Q= distortion_map(Q)
  print('count, t', [count, t])
  while i > -1:
      j = w
      while j > 0:
          ell = V._line_(V, Q)
          t = (t**2)*ell
          V = 2*V
          j -= 1
          count = 2*count
          print('count, t', [count, t])

      j = i-w+1
      m = 0
      while j <= i:
          m += nbin[j]*2^(j-i+w-1)
          j+=1

      if m != 0:
          ell = V._line_(Pi[m], Q)
          t = t*fi[m]*ell
          V = V + Pi[m]
          count = count + m
          print('count, t', [count, t])

      i = i-w
  if n_is_negative:
      vee = V._line_(-V, Q)
      t = 1/(t*vee)
  print('count, t', [count, t])
  return t

def precomptation_scalar(P, w):
  one = P.curve().base_field().one()
  t = one
  p = [0]
  p.insert(1, P)
  for i in range(2, (2^w)):
   p.insert(i, p[i-1] + p[1])
  return p

def precomptation_bkls_func(P, Q, w):
  p = [0]
  p.insert(1, P)
  t = [0]
  t.insert(1,1)
  # print('count, t', [count, t])
  Q = distortion_map(Q)
  for i in range(2, (2^w)):
    p.insert(i, p[i-1] + p[1])
    ell = p[1]._line_(p[i-1], Q)
    t.insert(i, t[i-1]*ell)

  return t

def distortion_map(P):
  global s
  global t
  global Ex
  global Fx
  global phi
  x = phi(P.xy()[0])
  y = phi(P.xy()[1])
  x_2 = x+ s^2
  y_2 = y+s*x+t
  # 2^m 拡大体だと、errorが出てしまった
  return Ex(x_2,y_2)

# ペアリング演算の高速化に関する研究（定義などの参考).pdf を参考にパラメータを選定

m=5
q=2^m
F.<a>=GF(2^m) # dunderlying field
E = EllipticCurve(F,[0, 0, 1, 1, 0])
Fx.<b>=GF(2^(4*m))
Ex = EllipticCurve(Fx,[0, 0, 1, 1, 0])
t=Fx(math.sqrt(q))
s=Fx(-(t^2 + t))
print(t)
print(s)
phi=Hom(F,Fx)(F.gen().minpoly().roots(Fx)[0][0])

P=E(a^3, a^2 + a)
Q=E(a^3 + 1, a^3 + a^2 + a)
n=25
w=2
print('bkls miller')
bkls_miller=bkls_miller(P, Q, n)

Pi = precomptation_scalar(P, w)
fi = precomptation_bkls_func(P, Q, w)

print('window bkls miller')
window_bkls_miller = window_bkls_miller(P, Q, n, Pi, fi, w)

print(bkls_miller==window_bkls_miller)
r"""
sage: E_2.points()
[(0 : 0 : 1), (0 : 1 : 0), (0 : 1 : 1), (1 : 0 : 1), (1 : 1 : 1)]
sage: arr=E_2.points()
sage: map(lambda x:x.order(), arr)
[5, 1, 5, 5, 5]
"""
