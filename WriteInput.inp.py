import os
for i in os.listdir():
    if '.bin' in i:
        print(i)
dict = {'y':True, 'Y':True, 'YES': True, 'NO':False, 'n':False}
s = dict[input('Do you want to write the input file?\n')]
if s:
    name = input('Write the name of the Input file\n')
    f = open(name, 'w')
    s = int(input('Dimension of the problem:\n'))
    f.write(f'{s}\n')
    s = input('Density matrix:\n')
    f.write(s+'\n')
    s = input('Hamiltonian:\n')
    f.write(s + '\n')
    s = input('Coupling 1:\n')
    f.write(s+'\n')
    s = input('Coupling 2:\n')
    f.write(s+'\n')
    n = int(input('How many properties?'))
    f.write(f'{n}\n')
    for i in range(n):
        s = input(f'Property{i+1}:\n')
        f.write(s + '\n')
    print('End!')
s = dict[input('Do you want to write the input for the dynamics?\n')]
if s:
    name = input('Write the name of the Input file\n')
    f = open(name, 'w')
    s = input('t0 t1:\n')
    f.write(f'{s}\n')
    s = input('How may steps:\n')
    f.write(s+'\n')
    s = input('When do you want to save ACHTUNG 1 FOR NOW:\n')
    f.write(s+'\n')
    s = input('What algorithm?\n'
              'A:Arnoldi\tD:Diagonalization\t RKR:Runge Kutta\n RKU:Runge Kutta Unitary\n')
    f.write(s+'\n')
    s = input('Temperature:\n')
    f.write(s+'\n')
    s = input('What distribution? BE or Dark for now\n')
    f.write(s + '\n')
    s = input('What distribution? \n'
              'null:null\tconst:constant\tDb:debye')
    f.write(s + '\n')
    s = input('What parameters?\n')
    f.write(s + '\n')
    print('End!')
