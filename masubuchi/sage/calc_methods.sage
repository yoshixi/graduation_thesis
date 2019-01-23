def tadd(k):
  return 4*k**2 + 6*k + 14

def tdbl(k):
  return 6*k**2 + 6*k + 14

def b_tadd(k):
  return k**2 + 3*k + 13

def b_tdbl(k):
  return 2*k**2 + 3*k + 13

def miller(k, l, w):
  return round(l * tdbl(k) + (l/2) * tadd(k), 2)

def w_miller(k, l, w):
  return round(l * tdbl(k) + (l/w + 2**w - 2) * tadd(k),2)

def bkls(k, l, w):
  return round(l * b_tdbl(k) + (l/2) * b_tadd(k), 2)

def w_bkls(k, l, w):
  return round(l * b_tdbl(k) + (l/w + 2**w - 2) * b_tadd(k), 2)

def main(k, l, w):
  var_bkls = bkls(k, l, w)
  var_w_bkls = w_bkls(k, l, w)
  print('k: ', k)
  print('l: ', l)
  print('m: ', w)
  print('miller :', miller(k, l, w))
  print('w_miller :', w_miller(k, l, w))
  print('bkls :', var_bkls)
  print('w_bkls :', var_w_bkls)
  print('w_bkls / bkls ', round(1 - (var_w_bkls / var_bkls), 3))
