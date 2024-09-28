# Original circuit
#   qubits: 9
#   qubits used: 9
#   H: 14
#   cnot: 84
#   X: 0
#   T: 84
#   P: 0
#   Z: 0
#   tdepth (by partitions): 36
#   depth  (by critical paths): 121
#   tdepth (by critical paths): 36
# Optimized circuit
#   qubits: 9
#   qubits used: 9
#   H: 14
#   cnot: 134
#   X: 0
#   T: 40
#   P: 18
#   Z: 2
#   tdepth (by partitions): 18
#   depth  (by critical paths): 123
#   tdepth (by critical paths): 18
#   Time: 0.001 s
.v 1 2 3 4 5 6 7 8 9
.i 1 2 3 4 5 6 7 8 9
.o 1 2 3 4 5 6 7 8 9

BEGIN
H 9
tof 6 3
tof 7 4
tof 8 5
tof 8 9
T 8
T* 5
T* 9
T 7
T* 4
T 6
T* 3
tof 5 8
tof 9 8
tof 9 5
tof 4 7
tof 6 3
tof 7 4
tof 5 9
tof 8 9
tof 8 5
T 8
tof 5 8
tof 2 1
P 9
P* 8
P* 1
tof 8 9
tof 2 1
tof 9 8
H 5
tof 4 7
tof 5 7
tof 5 4
T* 4
T 7
tof 7 4
tof 5 7
tof 4 7
P 9
Z 4
P 5
Z 3
P 1
P 2
tof 5 4
P* 4
tof 5 4
H 7
tof 6 3
tof 7 3
tof 7 6
T* 6
T 3
tof 3 6
tof 7 3
tof 6 3
tof 7 6
P 7
P* 6
tof 6 7
tof 7 6
H 3
tof 3 1
tof 3 2
T 3
T* 1
T* 2
tof 1 3
tof 2 3
tof 1 2
tof 3 2
tof 3 1
T 3
tof 1 3
tof 2 3
H 2
tof 2 7
tof 6 7
tof 6 2
tof 2 6
T* 2
T 7
tof 7 2
tof 6 7
tof 2 7
tof 6 2
P 6
P* 2
tof 2 6
tof 6 2
H 7
tof 7 4
tof 5 4
tof 5 7
tof 7 5
T* 7
T 4
tof 4 7
tof 5 4
tof 7 4
tof 5 7
P 5
P* 7
tof 7 5
tof 5 7
H 4
tof 4 9
tof 8 4
tof 4 8
T 8
T* 9
T* 4
tof 9 8
tof 4 8
tof 9 4
tof 8 4
tof 8 9
T 8
tof 9 8
tof 4 8
H 9
H 4
tof 7 5
tof 4 5
tof 4 7
T* 7
T 5
tof 5 7
tof 4 5
tof 7 5
tof 4 7
P 4
P* 7
tof 7 4
tof 4 7
H 5
tof 2 6
tof 5 6
tof 5 2
T* 2
T 6
tof 6 2
tof 5 6
tof 2 6
tof 5 2
P 5
P* 2
tof 2 5
tof 5 2
H 6
tof 6 3
tof 6 1
T 6
T* 3
T* 1
tof 3 6
tof 1 6
tof 3 1
tof 6 1
tof 6 3
T 6
tof 3 6
tof 1 6
H 1
tof 1 5
tof 1 2
T 1
T* 5
T* 2
tof 5 1
tof 2 1
tof 5 2
tof 1 2
tof 1 5
T 1
tof 5 1
tof 2 1
H 5
tof 5 4
tof 5 7
T 5
T* 4
T* 7
tof 4 5
tof 7 5
tof 4 7
tof 5 7
tof 5 4
T 5
tof 4 5
tof 7 5
H 4
tof 6 2
tof 2 6
tof 6 2
tof 2 3
tof 3 2
tof 2 3
tof 3 1
tof 1 3
tof 3 1
tof 5 8
tof 8 5
tof 5 8
tof 8 4
tof 4 8
tof 8 4
END
