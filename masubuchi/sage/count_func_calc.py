# -*- coding: utf-8 -*-
import math
# input ex = 29 = [1, 1, 1, 0, 1]
# output 29

def default_miller():
  a = 1
  f = 1
  # n = [1, 0, 1, 1, 1]
  num = base_10_to_n(16, 2)
  arr_num = num.split(',')
  n = list(map(lambda x: int(x), arr_num))
  print(n)
  w =  3
  l = len(n)
  for i in reversed(range(l-1)):
    #print(i)
    f = f * 2
    a = a * 2
    if n[i] == 1:
      f = f + 1
      a = a + 1
  print('default-miller: ', a)

def proposed_window_miller():
    p = 1
    w = 2
    m =  2**w
    # num = base_10_to_n(22, 2)
    # arr_num = num.split(',')
    # n = list(map(lambda x: int(x), arr_num))
    n = [0 ,0, 0, 1]
    l = len(n)
    print('n: ', n)
    print('l: ', l)
    i = l - 1
    while i > -1:
        print('==================roop: i', i)
        step = w
        while step > 0:
            p = 2*p
            step -= 1
        m = 0
        j = w - 1
        while j > -1:
            if i - j < 0:
                break
            m += n[i - j]
            print('i-j:', i-j)
            j -= 1
        if m !=0:
            p = p + m
            print('p: ', p)
        i = i - w
    print(p)


def window_miller():
  # m = 2^w と定義している
  a = 1
  f = 1
  w = 2
  m =  2**w
  num = base_10_to_n(22, 2)
  arr_num = num.split(',')
  n = list(map(lambda x: int(x), arr_num))
  l = len(n)
  print('n: ', n)
  step = math.floor(l/w)
  print('step: ', step)

  max_i = w*step
  l_max_i = l - max_i
  print('max_i: ', max_i)

  i = max_i
  while i > -1:
    # if i == max_i:
    #     for k in range(l_max_i):
    #         f = f * 2
    #         a = a * 2
    # else:
    #     for k in range(w -1):
    #         f = f * 2
    #         a = a * 2
    m = 0
    if i == max_i:
        width = l_max_i
    else:
        width = w

    for j in reversed(range(i,i+width)):
      m += n[j] * (2 ** (j))
    if m != 0:
      f = f + m
      a = a + m
    print(a)
    if step == 0:
      break
    step -= 1
    i = i - w
    print('loop end')
    print('====================')
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

def binary_window_method():
  # d: bit列
  # for

  # 95 = (1011111)_2

  m = 4
  w = 2
  num = base_10_to_n(41, 2)
  arr_num = num.split(',')
  n = list(map(lambda x: int(x), arr_num))
  l = len(n)
  q = 0
  print('n: ', n)
  for j in reversed(range(0,l,w)):
    q = (m)*q
    if j+1 < l:
        n_2 = n[j+1]
    else:
        n_2 = 0
    q = q + n[j] + n_2*2
  print('binary_window_msthod', q)


def binary_window_method_while_loop():
  # d: bit列
  # for
  # 95 = (1011111)_2

  w = 2
  m = 2**w
  num = base_10_to_n(41, 2)
  arr_num = num.split(',')
  n = list(map(lambda x: int(x), arr_num))
  l = len(n)
  q = 0
  j = l-2
  while j > -1:
    # q = (m)*q
    print('j: ', j)
    print('q: ', q)
    k = w
    while k > 0:
        q = 2*q
        k -= 1

    if j+1 < l:
        n_2 = n[j+1]
    else:
        n_2 = 0
    m = n[j] + n_2*2
    print('m: ', m)
    q = q + m
    j -= w
  print('while loop q:', q)


def binary_scalar_add_method():
  num = base_10_to_n(41, 2)
  arr_num = num.split(',')
  n = list(map(lambda x: int(x), arr_num))
  q = 0
  l = len(n)
  for j in reversed(range(l)):
    print('j: ', j)
    print('q: ', q)
    q = 2 * q
    if n[j] == 1:
      q = q + 1
  print('binary_scalar_add_method', q)

def base_10_to_n(X, n):
    if (int(X/n)):
        return str(X%n)+','+base_10_to_n(int(X/n), n)
    return str(X%n)

if __name__ == "__main__":
  # default_miller()
  # window_miller()
  # window_method()
  binary_scalar_add_method()
  # binary_window_method()
  binary_window_method_while_loop()
  # binary_method()
  # proposed_window_miller()
  # num = base_10_to_n(95, 2)
  # arr_num = num.split(',')
  # arr = list(map(lambda x: int(x), arr_num))
  # print(arr)
  #  a = base_10_to_n(19, 2)
  #  print('2: ', a)
  #  b = base_10_to_n(19, 4)
  #  print('4: ', b)
