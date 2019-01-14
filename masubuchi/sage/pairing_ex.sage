# p = 103; A = 1; B = 18; E = EllipticCurve(GF(p), [A, B])
# P = E(33, 91); n = P.order(); n
# k = GF(n)(p).multiplicative_order(); k
#
# e = Integer((p^k-1)/n); e
# print(P.weil_pairing(Q, n)^e)
# print(P._miller_(Q, n)^e)
# P.tate_pairing(Q, n, k) == P._miller_(Q, n)^e
# Q.tate_pairing(P, n, k) == Q._miller_(P, n)^e
# P.tate_pairing(Q, n, k)/Q.tate_pairing(P, n, k)

# sample
# t = cputime()
# p = 103; A = 1; B = 18; E = EllipticCurve(GF(p), [A, B])
# P = E(33, 91); n = P.order(); n
# k = GF(n)(p).multiplicative_order(); k
# e = (p^k-1)/n; e
# print(p^k-1)
# print(n)
# print(e)
# print(P.tate_pairing(P, n, k))
# Q = E(87, 51)
# P.tate_pairing(Q, n, k)
# set_random_seed(35)
# P.tate_pairing(P,n,k)
# print(cputime(t))

# https://ask.sagemath.org/question/34642/how-to-implement-weil-pairing/
# t = cputime()
# F.<a>=GF(2^5)
# E=EllipticCurve(F,[0,0,1,1,1])
# P = E(a^4 + 1, a^3)
# Fx.<b>=GF(2^(4*5))
# Ex=EllipticCurve(Fx,[0,0,1,1,1])
# phi=Hom(F,Fx)(F.gen().minpoly().roots(Fx)[0][0])
# Px=Ex(phi(P.xy()[0]),phi(P.xy()[1]))
# O = Ex(0)
# Qx = Ex(b^19 + b^18 + b^16 + b^12 + b^10 + b^9 + b^8 + b^5 + b^3 + 1, b^18 + b^13 + b^10 + b^8 + b^5 + b^4 + b^3 + b)
# Px.weil_pairing(Qx,41) == b^19 + b^15 + b^9 + b^8 + b^6 + b^4 + b^3 + b^2 + 1
# Px.weil_pairing(17*Px,41) == Fx(1)
# Px.weil_pairing(O,41) == Fx(1)
# print(cputime(t))

def miller(P, Q, n):
  if Q.is_zero():
      raise ValueError("Q must be nonzero.")
  if n.is_zero():
      raise ValueError("n must be nonzero.")
  n_is_negative = False
  if n < 0:
      n = n.abs()
      n_is_negative = True

  one = P.curve().base_field().one()
  t = one
  V = P
  S = 2*V
  nbin = n.bits()
  i = n.nbits() - 2
  while i > -1:
      S = 2*V
      ell = V._line_(V, Q)
      vee = S._line_(-S, Q)
      t = (t**2)*(ell/vee)
      V = S
      if nbin[i] == 1:
          S = V+P
          ell = V._line_(P, Q)
          vee = S._line_(-S, Q)
          t = t*(ell/vee)
          V = S
      i = i-1
  if n_is_negative:
      vee = V._line_(-V, Q)
      t = 1/(t*vee)
  return t

# Finite Fields
# http://www.sagemath.org/files/thesis/hansen-thesis-2009.pdf

# http://doc.sagemath.org/pdf/ja/tutorial/tutorial-jp.pdf
# <a>は変数aを定義している。以下、29ページに記載されている。
# 楕円曲線　EllipticCurve([a1, a2, a3, a4, a6 ])の引数は52pに記載されている
F.<a>=GF(2^5)
E=EllipticCurve(F,[0,0,1,1,1])
P = E(a^4 + 1, a^3)
Fx.<b>=GF(2^(4*5))
Ex=EllipticCurve(Fx,[0,0,1,1,1])
phi=Hom(F,Fx)(F.gen().minpoly().roots(Fx)[0][0])
# print('phi: ',phi)
# print('F:',F)
# print('F.gen:',F.gen())
# print('F.minipoly:',F.gen().minpoly())
# print('F.roots:',F.gen().minpoly().roots(Fx))
# print('F:',F.gen().minpoly().roots(Fx)[0][0])
Px=Ex(phi(P.xy()[0]),phi(P.xy()[1]))
print('Px: ',Px.xy()[0])
print('Px: ',Px.xy()[1])
Qx = Ex(b^19 + b^18 + b^16 + b^12 + b^10 + b^9 + b^8 + b^5 + b^3 + 1, b^18 + b^13 + b^10 + b^8 + b^5 + b^4 + b^3 + b)
print('Px-miller: ', Px._miller_(Qx,41))
print('Qx-miller: ', Qx._miller_(Px,41))
print(Px._miller_(Qx,41) == b^17 + b^13 + b^12 + b^9 + b^8 + b^6 + b^4 + 1)
print(Qx._miller_(Px,41) == b^13 + b^10 + b^8 + b^7 + b^6 + b^5)
print(window_miller(Px,Qx,41)==Px._miller_(Qx,41))

nagative_number = -12
# print(typeof(Px))
# Px.miller(Qx,41)==b^17+b^13+b^12+b^9+b^8+b^6+b^4+1
# Qx.miller(Px,41)==b^13+b^10+b^8+b^7+b^6+b^5


# ('phi: ', Ring morphism:
#   From: Finite Field in a of size 2^5
#   To:   Finite Field in b of size 2^20
#   Defn: a |--> b^17 + b^16 + b^15 + b^13 + b^10 + b^9 + b^7 + b^6 + b)

E = EllipticCurve('11a'); E
P = E(0); P
P.division_points(5)
