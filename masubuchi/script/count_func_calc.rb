
# input ex = 29 = [1, 1, 1, 0, 1]
# output 29


a = 1
f = 1
n = [1, 1, 1, 0, 1]
w =  3

n.each_with_index do |el, idx|
 f = f * 2
 a = a * 2

 m = 0
 for j in  (i -w +1)...i do
   m += n[j]
 end
end

    
