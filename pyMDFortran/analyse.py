atoms =   2125764 # nb = 81
pairs = 161836128
bytes_per_atom  = 8*9 + 4*pairs/atoms # bytes per atom, 8B*(6 reads + 3 write) + 4B*(#interactions per atom)
print("atoms :",atoms)
print("pairs :",pairs)
print("B/atom:",bytes_per_atom)

bandwidth = 95 # GB/s as measured by intel mlc for 2 reads : 1 write
pairs_per_second = 88.6e6
atoms_per_second = pairs_per_second*atoms/pairs
# the minimum data transfer is
G = 1024**3

#how many atoms could we transfer per s
max_atoms_per_second = bandwidth*G/bytes_per_atom
print("maximum atoms per second:", max_atoms_per_second)
print(" actual atoms per second:", atoms_per_second)
print("                   ratio:", atoms_per_second/max_atoms_per_second)

flops_per_pair = 27 
flops_per_atom = flops_per_pair*pairs/atoms
flops_per_second = flops_per_pair*pairs_per_second
peak = 1*1*4*2.8
print("flops_per_pair   :",flops_per_pair)
print("flops_per_atom   :",flops_per_atom) 
print("gflops_per_second:",flops_per_second/1e9)
print("peak performance :",peak)
print("           ratio :",flops_per_second/1e9/peak)

