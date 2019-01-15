var phi=0
var Ex =0

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
      ell = V._line_(V, Ex(phi(Q.xy()[0]),phi(Q.xy()[1])))
      t = (t**2)*ell
      count = 2*count
      V = S
      if nbin[i] == 1:
          S = V+P
          ell = V._line_(P, Ex(phi(Q.xy()[0]),phi(Q.xy()[1])))
          t = t*ell
          count += 1
          V = S
      i = i-1
  if n_is_negative:
      vee = V._line_(-V, Q)
      t = 1/(t*vee)
  return t

# BN curve
# https://www.sci.kanagawa-u.ac.jp/info/matsuo/theses/pdf/matsuo_lab_20150203_kobayashi_abst.pdf
global phi
global Ex

z=5
k=2 #埋め込み次数
p = 36*(z^4) + 36*(z^3) + 24*(z^2) + 6*(z) + 1
n = 36*(z^4) + 36*(z^3) + 18*(z^2) + 6*(z) + 1 #torsion
F.<a>=GF(p) # Finite Field of size 27631
E=EllipticCurve(F,[0,0,0,0,28]) #楕円曲線を定義
Fx.<b>=GF(p^(k))  # Fの拡大を定義
Ex=EllipticCurve(Fx,[0,0,0,0,28]) # 拡大体上の楕円曲線を定義
phi=Hom(F,Fx)(F.gen().minpoly().roots(Fx)[0][0])
P = E.points()[1] # torsion points
Qx = Ex.random_point()
Q = E.points()[3]
miller = original_miller(P,Qx, n)
miller = original_miller(P,Qx, n)

print(miller)
