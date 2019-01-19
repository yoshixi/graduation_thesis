Exx =0
s=0
t=0

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
  print('count, t', [count, t])
  Q= distortion_map(Q)
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
  x = Fx(P.xy()[0])
  y = Fx(P.xy()[1])
  x_2 = x+ s^2
  y_2 = y+s*x+t
  # 2^m 拡大体だと、errorが出てしまった
  return Ex(x_2,y_2)

m=283
k=4 # 埋め込み字数
F.<a>=GF(2^m) # dunderlying field
E = EllipticCurve(F,[0, 0, 1, 1, 0])
Fx.<a>=GF((2^m)^k) #distortion_mapで飛ばすところ
Ex = EllipticCurve(Fx,[0, 0, 1, 1, 0])
Fxx.<b>=GF(2^(4*m))  # sの条件式を考える体
Exx = EllipticCurve(Fxx,[0, 0, 1, 1, 0])
s=Fxx(0)
t=Fxx(1)

P = E(a^282 + a^278 + a^276 + a^275 + a^272 + a^271 + a^270 + a^269 + a^268 + a^266 + a^264 + a^263 + a^262 + a^259 + a^258 + a^251 + a^242 + a^240 + a^239 + a^238 + a^236 + a^235 + a^231 + a^228 + a^226 + a^224 + a^222 + a^215 + a^214 + a^213 + a^207 + a^205 + a^204 + a^202 + a^201 + a^199 + a^198 + a^194 + a^190 + a^187 + a^186 + a^185 + a^184 + a^181 + a^179 + a^178 + a^177 + a^175 + a^174 + a^173 + a^172 + a^169 + a^168 + a^164 + a^163 + a^162 + a^161 + a^160 + a^159 + a^156 + a^154 + a^152 + a^149 + a^147 + a^146 + a^143 + a^140 + a^137 + a^136 + a^135 + a^133 + a^132 + a^131 + a^129 + a^127 + a^126 + a^125 + a^121 + a^119 + a^118 + a^117 + a^116 + a^112 + a^108 + a^107 + a^106 + a^103 + a^102 + a^101 + a^100 + a^99 + a^96 + a^95 + a^94 + a^92 + a^90 + a^89 + a^86 + a^85 + a^84 + a^83 + a^78 + a^77 + a^76 + a^75 + a^72 + a^71 + a^70 + a^68 + a^67 + a^63 + a^62 + a^61 + a^51 + a^47 + a^45 + a^44 + a^43 + a^40 + a^39 + a^38 + a^35 + a^33 + a^32 + a^28 + a^27 + a^25 + a^24 + a^23 + a^20 + a^19 + a^18 + a^15 + a^14 + a^9 + a^8 + a^7 + a^6 + a^5 + a^3, a^282 + a^280 + a^279 + a^274 + a^271 + a^269 + a^267 + a^264 + a^261 + a^258 + a^255 + a^249 + a^248 + a^245 + a^244 + a^239 + a^237 + a^236 + a^235 + a^233 + a^228 + a^222 + a^218 + a^217 + a^216 + a^214 + a^212 + a^211 + a^209 + a^207 + a^204 + a^198 + a^197 + a^194 + a^193 + a^192 + a^190 + a^189 + a^185 + a^183 + a^182 + a^180 + a^179 + a^178 + a^177 + a^175 + a^172 + a^170 + a^168 + a^166 + a^164 + a^163 + a^162 + a^161 + a^158 + a^155 + a^153 + a^152 + a^151 + a^149 + a^144 + a^142 + a^141 + a^139 + a^136 + a^135 + a^134 + a^132 + a^131 + a^130 + a^129 + a^124 + a^123 + a^121 + a^115 + a^114 + a^112 + a^111 + a^108 + a^107 + a^106 + a^104 + a^103 + a^100 + a^98 + a^96 + a^95 + a^94 + a^93 + a^92 + a^91 + a^83 + a^82 + a^81 + a^80 + a^78 + a^76 + a^71 + a^68 + a^66 + a^65 + a^62 + a^58 + a^57 + a^53 + a^51 + a^50 + a^49 + a^48 + a^45 + a^43 + a^40 + a^38 + a^33 + a^32 + a^30 + a^29 + a^27 + a^24 + a^22 + a^20 + a^18 + a^16 + a^15 + a^12 + a^11 + a^10 + a^9 + a^8 + a^7 + a^6 + 1)

Q = E(a^281 + a^279 + a^277 + a^275 + a^271 + a^266 + a^265 + a^264 + a^259 + a^258 + a^254 + a^253 + a^251 + a^250 + a^249 + a^246 + a^245 + a^244 + a^243 + a^241 + a^240 + a^238 + a^233 + a^232 + a^231 + a^230 + a^228 + a^226 + a^224 + a^220 + a^217 + a^216 + a^214 + a^213 + a^210 + a^209 + a^205 + a^204 + a^201 + a^200 + a^199 + a^198 + a^195 + a^194 + a^193 + a^190 + a^189 + a^187 + a^186 + a^180 + a^178 + a^175 + a^174 + a^172 + a^164 + a^162 + a^159 + a^158 + a^157 + a^156 + a^154 + a^151 + a^149 + a^148 + a^147 + a^145 + a^143 + a^141 + a^140 + a^139 + a^138 + a^137 + a^136 + a^135 + a^132 + a^129 + a^127 + a^126 + a^125 + a^123 + a^121 + a^118 + a^117 + a^114 + a^113 + a^112 + a^110 + a^109 + a^108 + a^107 + a^106 + a^105 + a^104 + a^102 + a^99 + a^98 + a^96 + a^94 + a^92 + a^89 + a^87 + a^86 + a^85 + a^84 + a^83 + a^80 + a^75 + a^73 + a^72 + a^70 + a^69 + a^67 + a^63 + a^62 + a^60 + a^59 + a^56 + a^55 + a^52 + a^50 + a^49 + a^47 + a^44 + a^39 + a^38 + a^37 + a^34 + a^33 + a^32 + a^31 + a^29 + a^24 + a^23 + a^22 + a^20 + a^18 + a^17 + a^15 + a^14 + a^12 + a^11 + a^10 + a^9 + a^7 + a^5 + a^3 + a^2 + a + 1, a^278 + a^276 + a^275 + a^273 + a^272 + a^271 + a^268 + a^265 + a^264 + a^262 + a^260 + a^259 + a^258 + a^252 + a^251 + a^250 + a^249 + a^243 + a^242 + a^241 + a^239 + a^238 + a^237 + a^236 + a^234 + a^232 + a^231 + a^229 + a^228 + a^227 + a^226 + a^221 + a^218 + a^216 + a^215 + a^210 + a^209 + a^208 + a^206 + a^203 + a^201 + a^199 + a^198 + a^196 + a^195 + a^194 + a^192 + a^189 + a^188 + a^187 + a^186 + a^184 + a^183 + a^182 + a^181 + a^179 + a^177 + a^176 + a^175 + a^174 + a^169 + a^168 + a^167 + a^166 + a^164 + a^163 + a^162 + a^161 + a^158 + a^153 + a^151 + a^149 + a^147 + a^146 + a^143 + a^141 + a^138 + a^137 + a^136 + a^134 + a^133 + a^132 + a^131 + a^129 + a^128 + a^127 + a^126 + a^125 + a^120 + a^119 + a^118 + a^117 + a^115 + a^113 + a^112 + a^109 + a^108 + a^107 + a^106 + a^105 + a^103 + a^101 + a^99 + a^96 + a^94 + a^92 + a^91 + a^90 + a^89 + a^87 + a^86 + a^83 + a^81 + a^79 + a^78 + a^77 + a^75 + a^73 + a^71 + a^70 + a^68 + a^63 + a^60 + a^59 + a^58 + a^57 + a^54 + a^52 + a^48 + a^47 + a^42 + a^40 + a^39 + a^37 + a^36 + a^34 + a^31 + a^29 + a^28 + a^25 + a^19 + a^18 + a^16 + a^15 + a^13 + a^12 + a^8 + a^7 + a^6 + a^4 + a^2 + a + 1)
n = 15541351137805832567355695254588151253139249137230816537358713893981666119551291490305
w = 2

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
