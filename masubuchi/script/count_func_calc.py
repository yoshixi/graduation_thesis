# -*- coding: utf-8 -*-

# input ex = 29 = [1, 1, 1, 0, 1]
# output 29

def default_miller():
  a = 1
  f = 1
  # n = [1, 0, 1, 1, 1]
  n = [1, 1, 1, 1, 1, 0, 1]
  w =  3
  l = len(n)
  for i in reversed(range(l-1)):
    #print(i)
    f = f * 2
    a = a * 2
    if n[i] == 1:
      f = f + 1
      a = a + 1
    print(a)

def window_miller():
  # m = 2^w と定義している
  a = 1
  f = 1
  m =  4
  num = base_10_to_n(95, m)
  arr_num = num.split(',')
  n = list(map(lambda x: int(x), arr_num))
  l = len(n)
  print(n)

  for i in reversed(range(l)):
    print('i:', i)
    f = f * 2
    a = a * 2

    m = 0
    for j in range(i-w+1, i):
      m += n[j] * (2 ** (j-i+w-1))
      print('m: ', m)
    if m != 1:
      f = f + m
      a = a + m
    print(a)

def window_method():
  # d: bit列
  # for

  # 95 = (1011111)_2

  m = 4
  num = base_10_to_n(95, 4)
  arr_num = num.split(',')
  n = list(map(lambda x: int(x), arr_num))
  print(n)
  l = len(n)
  q = 0

  for j in reversed(range(l)):
    print('j: ', j)
    q = (m)*q
    q = q + n[j]
  print(q)

def binary_method():
  q = 0
  n = [1, 1, 1, 1, 1, 0, 1]
  l = len(n)
  for j in reversed(range(l)):
    q = 2 * q
    if n[j] == 1:
      q = q + 1
  print(q)

def base_10_to_n(X, n):
    if (int(X/n)):
        return str(X%n)+','+base_10_to_n(int(X/n), n)
    return str(X%n)

if __name__ == "__main__":
  # default_miller()
  # window_miller()
  window_method()
  # binary_method()

  # num = base_10_to_n(95, 2)
  # arr_num = num.split(',')
  # arr = list(map(lambda x: int(x), arr_num))
  # print(arr)
