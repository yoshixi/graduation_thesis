import time

phi=0
Ex =0

def original_miller(P, Q, n):
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
  S = 2*V
  nbin = n.bits()
  i = n.nbits() - 2
  count = 1
  while i > -1  :
      S = 2*V
      ell = V._line_(V, Q)
      vee = S._line_(-S, Q)
      t = (t**2)*(ell/vee)
      count = 2*count
      V = S
      if nbin[i] == 1:
          S = V+P
          ell = V._line_(P, Q)
          vee = S._line_(-S, Q)
          t = t*(ell/vee)
          count += 1
          V = S
      i = i-1
  if n_is_negative:
      vee = V._line_(-V, Q)
      t = 1/(t*vee)
  return t

def bkls_miller(P, Q, n):
  global phi
  global Ex
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
  S = 2*V
  nbin = n.bits()
  i = n.nbits() - 2
  count = 1
  Q = Ex(phi(Q.xy()[0]),phi(Q.xy()[1]))
  while i > -1  :
      S = 2*V
      ell = V._line_(V, Q)
      t = (t**2)*ell
      count = 2*count
      V = S
      if nbin[i] == 1:
          S = V+P
          ell = V._line_(P, Q)
          t = t*ell
          count += 1
          V = S
      i = i-1
  if n_is_negative:
      vee = V._line_(-V, Q)
      t = 1/(t*vee)
  return t


def window_bkls_miller(P, Q, n, pi, fi, w):
  global phi
  global Ex
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
  S = 2*V
  nbin = n.bits()
  i = n.nbits() - 2
  count = 1
  Q = Ex(phi(Q.xy()[0]),phi(Q.xy()[1]))
  while i > -1:
      j = w
      while j > 0:
          S = 2*V
          ell = V._line_(V, Q)
          vee = S._line_(-S, Q)
          t = (t**2)*(ell/vee)
          V = S
          j -= 1

      j = i-w+1
      m = 0
      while j <= i:
          m += nbin[j]*2^(j-i+w-1)
          j+=1

      if m != 0:
          S = V+pi[m]
          ell = V._line_(pi[m], Q)
          vee = S._line_(-S, Q)
          t = t*fi[m]*(ell/vee)
          V = V + pi[m]
      i = i-w
  if n_is_negative:
      vee = V._line_(-V, Q)
      t = 1/(t*vee)
  return t

def precomptation_scalar(P, w):
  one = P.curve().base_field().one()
  t = one
  p = [0]
  p.insert(1, P)
  for i in range(2, (2^w)):
   p.insert(i, p[i-1] + p[1])
  return p

def precomptation_func(P, Q, w):
  p = [0]
  p.insert(1, P)
  t = [0]
  t.insert(1,1)
  count = 1
  for i in range(2, (2^w)):
    p.insert(i, p[i-1] + p[1])
    ell = p[1]._line_(p[i-1], Q)
    vee = p[i]._line_(-p[i], Q)
    t.insert(i, t[i-1]*(ell/vee))
    count += 1

  return t
# 超楕円
# http://doc.sagemath.org/html/en/reference/curves/sage/schemes/hyperelliptic_curves/hyperelliptic_generic.html


k=2 #埋め込み次数
# z=5
# p = 36*(z^4) + 36*(z^3) +    24*(z^2) + 6*(z) + 1
# n = 36*(z^4) + 36*(z^3) + 18*(z^2) + 6*(z) + 1
p = 11
n = 6
F.<a>=GF(p) # Finite Field of size 27631
E=EllipticCurve(F,[0,0,0,0,1]) #楕円曲線を定義
Fx.<b>=GF(p^(k))  # Fの拡大を定義
Ex=EllipticCurve(Fx,[0,0,0,0,1]) # 拡大体上の楕円曲線を定義
phi=Hom(F,Fx)(F.gen().minpoly().roots(Fx)[0][0])

P = E(2,3) # 6 torsion points
Qx = Ex.random_point()
Q = E(2,8)

start = time.time()
miller = original_miller(P, Qx, n)
process_time = time.time() - start

start2 = time.time()
bkls_miller = bkls_miller(P, Q, n)
process_time2 = time.time() - start2

print(miller)
print(bkls_miller)

print(process_time)
print(process_time2)
