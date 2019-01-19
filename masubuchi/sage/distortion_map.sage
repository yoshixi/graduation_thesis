
m=2
F.<a>=GF(2)
E_2 = EllipticCurve(F,[0, 0, 1, 1, 1])
Fx.<a>=GF(2^m)
Ex = EllipticCurve(Fx,[0, 0, 1, 1, 1])
Fst=GF(2^(4*m))
s=Fst(0)
t=F(-1)


  def distortion_map(P):
    x = P.xy()[0]+ s^2
    y = P.xy()[1]+s*x+t
    return Ex(x,y)




s, t = var('s, t')
solve([s+s^4==0, t+t^2+s^6+s^2==0], s, t)
