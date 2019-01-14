def base_10_to_n(X, n):
    if (int(X/n)):
        return str(X%n)+','+base_10_to_n(int(X/n), n)
    return str(X%n)

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
  print(count)
  if n_is_negative:
      vee = V._line_(-V, Q)
      t = 1/(t*vee)
  return t

def window_miller(P, Q, n, pi, fi, w):
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
  print(t)
  V = P
  S = 2*V
  nbin = n.bits()
  i = n.nbits() - 2
  count = 1
  while i > -1  :
      j = w
      while j > 0:
          S = 2*V
          ell = V._line_(V, Q)
          vee = S._line_(-S, Q)
          t = (t**2)*(ell/vee)
          count = 2*count
          V = S
          j -= 1

    #  j = i-w+1
    #  m = 0
    #  while j <= i:
    #      m += nbin[j]*2^(j-i+w-1)
    #      j+=1

      if i+1 < n.nbits():
          n_2 = nbin[i+1]
      else:
          n_2 = 0
      m = nbin[i] + n_2*2
      print('m: ', m)
      if m != 0:
          S = V+pi[m]
          ell = V._line_(pi[m], Q)
          vee = S._line_(-S, Q)
          t = t*fi[m]*(ell/vee)
          count += m
          V = V + pi[m]
      i = i-w
  if n_is_negative:
      vee = V._line_(-V, Q)
      t = 1/(t*vee)
  print(count)
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
# print('Px: ',Px.xy()[0])
# print('Px: ',Px.xy()[1])
Qx = Ex(b^19 + b^18 + b^16 + b^12 + b^10 + b^9 + b^8 + b^5 + b^3 + 1, b^18 + b^13 + b^10 + b^8 + b^5 + b^4 + b^3 + b)
#print('Px-miller: ', Px._miller_(Qx,41))
#print('Qx-miller: ', Qx._miller_(Px,41))
#print(Px._miller_(Qx,41) == b^17 + b^13 + b^12 + b^9 + b^8 + b^6 + b^4 + 1)
#print(Qx._miller_(Px,41) == b^13 + b^10 + b^8 + b^7 + b^6 + b^5)

# precomptation_funcがちゃんと計算できているか
n = 41
miller = original_miller(Px,Qx, n)

w = 2
pi = precomptation_scalar(Px, w)
fi = precomptation_func(Px,Qx, w)
w_miller = window_miller(Px,Qx, n, pi, fi, w)
print(miller==w_miller)