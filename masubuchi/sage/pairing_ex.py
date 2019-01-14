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
t = cputime()
F.<a>=GF(2^5)
E=EllipticCurve(F,[0,0,1,1,1])
P = E(a^4 + 1, a^3)
Fx.<b>=GF(2^(4*5))
Ex=EllipticCurve(Fx,[0,0,1,1,1])
phi=Hom(F,Fx)(F.gen().minpoly().roots(Fx)[0][0])
Px=Ex(phi(P.xy()[0]),phi(P.xy()[1]))
O = Ex(0)
Qx = Ex(b^19 + b^18 + b^16 + b^12 + b^10 + b^9 + b^8 + b^5 + b^3 + 1, b^18 + b^13 + b^10 + b^8 + b^5 + b^4 + b^3 + b)
Px.weil_pairing(Qx,41) == b^19 + b^15 + b^9 + b^8 + b^6 + b^4 + b^3 + b^2 + 1
Px.weil_pairing(17*Px,41) == Fx(1)
Px.weil_pairing(O,41) == Fx(1)
print(cputime(t))
