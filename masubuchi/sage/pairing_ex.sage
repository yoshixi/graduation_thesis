p = 103; A = 1; B = 18; E = EllipticCurve(GF(p), [A, B])
P = E(33, 91); n = P.order(); n
k = GF(n)(p).multiplicative_order(); k

e = Integer((p^k-1)/n); e
print(P.weil_pairing(Q, n)^e)
print(P._miller_(Q, n)^e)
P.tate_pairing(Q, n, k) == P._miller_(Q, n)^e
Q.tate_pairing(P, n, k) == Q._miller_(P, n)^e
P.tate_pairing(Q, n, k)/Q.tate_pairing(P, n, k)
