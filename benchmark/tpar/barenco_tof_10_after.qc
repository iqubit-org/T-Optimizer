# Original circuit
#   qubits: 19
#   qubits used: 19
#   H: 34
#   cnot: 224
#   X: 0
#   T: 224
#   P: 0
#   Z: 0
#   tdepth (by partitions): 96
#   depth  (by critical paths): 321
#   tdepth (by critical paths): 96
# Optimized circuit
#   qubits: 19
#   qubits used: 19
#   H: 34
#   cnot: 332
#   X: 0
#   T: 100
#   P: 48
#   Z: 7
#   tdepth (by partitions): 43
#   depth  (by critical paths): 272
#   tdepth (by critical paths): 43
#   Time: 0.006 s
.v 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19
.i 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19
.o 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19

BEGIN
H 19
tof 11 3
tof 12 4
tof 13 5
tof 14 6
tof 15 7
tof 16 8
tof 17 9
tof 18 10
tof 19 18
tof 18 19
T 19
T* 10
T* 18
T 17
T* 9
T 16
T* 8
T 15
T* 7
T 14
T* 6
T 13
T* 5
T 12
T* 4
T 11
T* 3
tof 18 10
tof 10 19
tof 9 17
tof 8 16
tof 7 15
tof 3 11
tof 11 3
tof 12 4
tof 13 5
tof 14 6
tof 15 7
tof 16 8
tof 17 9
tof 19 10
tof 10 18
T 19
tof 10 19
tof 2 1
P 18
P* 19
P* 1
tof 19 18
tof 2 1
tof 18 19
H 10
tof 9 17
tof 10 17
tof 10 9
T* 9
T 17
tof 17 9
tof 10 17
tof 9 17
P 18
Z 9
P 10
Z 16
Z 15
Z 6
Z 5
Z 4
Z 11
P 1
P 2
tof 10 9
P* 9
tof 10 9
H 17
tof 8 16
tof 17 16
tof 17 8
T* 8
T 16
tof 16 8
tof 17 16
tof 8 16
tof 17 8
P 17
P* 8
tof 8 17
tof 17 8
H 16
tof 7 15
tof 16 15
tof 16 7
T* 7
T 15
tof 15 7
tof 16 15
tof 7 15
tof 16 7
P 16
P* 7
tof 7 16
tof 16 7
H 15
tof 14 6
tof 15 6
tof 15 14
T* 14
T 6
tof 6 14
tof 15 6
tof 14 6
tof 15 14
P 15
P* 14
tof 14 15
tof 15 14
H 6
tof 13 5
tof 6 5
tof 6 13
T* 13
T 5
tof 5 13
tof 6 5
tof 13 5
tof 6 13
P 6
P* 13
tof 13 6
tof 6 13
H 5
tof 12 4
tof 5 4
tof 5 12
T* 12
T 4
tof 4 12
tof 5 4
tof 12 4
tof 5 12
P 5
P* 12
tof 12 5
tof 5 12
H 4
tof 3 11
tof 4 11
tof 4 3
T* 3
T 11
tof 11 3
tof 4 11
tof 3 11
tof 4 3
P 4
P* 3
tof 3 4
tof 4 3
H 11
tof 11 1
tof 11 2
T 11
T* 1
T* 2
tof 2 11
tof 1 11
tof 11 1
tof 1 2
T 11
tof 2 11
tof 1 11
H 2
tof 2 4
tof 3 4
tof 3 2
tof 2 3
T* 2
T 4
tof 4 2
tof 3 4
tof 2 4
tof 3 2
P 3
P* 2
tof 2 3
tof 3 2
H 4
tof 4 5
tof 12 5
tof 12 4
tof 4 12
T* 4
T 5
tof 5 4
tof 12 5
tof 4 5
tof 12 4
P 12
P* 4
tof 4 12
tof 12 4
H 5
tof 5 6
tof 13 6
tof 13 5
tof 5 13
T* 5
T 6
tof 6 5
tof 13 6
tof 5 6
tof 13 5
P 13
P* 5
tof 5 13
tof 13 5
H 6
tof 6 15
tof 14 15
tof 14 6
tof 6 14
T* 6
T 15
tof 15 6
tof 14 15
tof 6 15
tof 14 6
P 14
P* 6
tof 6 14
tof 14 6
H 15
tof 15 16
tof 7 16
tof 7 15
tof 15 7
T* 15
T 16
tof 16 15
tof 7 16
tof 15 16
tof 7 15
P 7
P* 15
tof 15 7
tof 7 15
H 16
tof 16 17
tof 8 17
tof 8 16
tof 16 8
T* 16
T 17
tof 17 16
tof 8 17
tof 16 17
tof 8 16
P 8
P* 16
tof 16 8
tof 8 16
H 17
tof 17 9
tof 10 9
tof 10 17
tof 17 10
T* 17
T 9
tof 9 17
tof 10 9
tof 17 9
tof 10 17
P 10
P* 17
tof 17 10
tof 10 17
H 9
tof 9 18
tof 9 19
T 9
T* 18
T* 19
tof 19 9
tof 18 9
tof 9 18
tof 18 19
T 9
tof 18 9
tof 19 9
H 18
H 19
tof 17 10
tof 19 10
tof 19 17
T* 17
T 10
tof 10 17
tof 19 10
tof 17 10
tof 19 17
P 19
P* 17
tof 17 19
tof 19 17
H 10
tof 16 8
tof 10 8
tof 10 16
T* 16
T 8
tof 8 16
tof 10 8
tof 16 8
tof 10 16
P 10
P* 16
tof 16 10
tof 10 16
H 8
tof 15 7
tof 8 7
tof 8 15
T* 15
T 7
tof 7 15
tof 8 7
tof 15 7
tof 8 15
P 8
P* 15
tof 15 8
tof 8 15
H 7
tof 6 14
tof 7 14
tof 7 6
T* 6
T 14
tof 14 6
tof 7 14
tof 6 14
tof 7 6
P 7
P* 6
tof 6 7
tof 7 6
H 14
tof 5 13
tof 14 13
tof 14 5
T* 5
T 13
tof 13 5
tof 14 13
tof 5 13
tof 14 5
P 14
P* 5
tof 5 14
tof 14 5
H 13
tof 4 12
tof 13 12
tof 13 4
T* 4
T 12
tof 12 4
tof 13 12
tof 4 12
tof 13 4
P 13
P* 4
tof 4 13
tof 13 4
H 12
tof 2 3
tof 12 3
tof 12 2
T* 2
T 3
tof 3 2
tof 12 3
tof 2 3
tof 12 2
P 12
P* 2
tof 2 12
tof 12 2
H 3
tof 3 11
tof 3 1
T 3
T* 11
T* 1
tof 1 3
tof 11 3
tof 3 11
tof 11 1
T 3
tof 1 3
tof 11 3
H 1
tof 1 12
tof 2 1
tof 1 2
T 2
T* 12
T* 1
tof 1 2
tof 12 2
tof 2 12
tof 12 1
T 2
tof 12 2
tof 1 2
H 12
tof 12 13
tof 12 4
T 12
T* 13
T* 4
tof 4 12
tof 13 12
tof 12 13
tof 13 4
T 12
tof 13 12
H 13
tof 13 14
tof 13 5
T 13
T* 14
T* 5
tof 5 13
tof 14 13
tof 13 14
tof 14 5
T 13
tof 14 13
H 14
tof 14 7
tof 14 6
T 14
T* 7
T* 6
tof 6 14
tof 7 14
tof 14 7
tof 7 6
T 14
tof 7 14
H 7
tof 7 8
tof 15 7
tof 7 15
T 15
T* 8
T* 7
tof 7 15
tof 8 15
tof 15 8
tof 8 7
T 15
tof 8 15
H 8
tof 8 10
tof 16 8
tof 8 16
T 16
T* 10
T* 8
tof 8 16
tof 10 16
tof 16 10
tof 10 8
T 16
tof 10 16
H 10
tof 10 19
tof 17 10
tof 10 17
T 17
T* 19
T* 10
tof 10 17
tof 19 17
tof 17 19
tof 19 10
T 17
tof 19 17
tof 10 17
H 19
tof 3 2
tof 2 3
tof 3 2
tof 2 11
tof 11 2
tof 2 11
tof 11 1
tof 1 11
tof 11 1
tof 9 17
tof 17 9
tof 9 17
tof 17 10
tof 10 17
tof 17 10
tof 12 4
tof 4 12
tof 13 5
tof 5 13
tof 14 6
tof 6 14
tof 15 7
tof 7 15
tof 16 8
tof 8 16
tof 19 18
tof 18 19
tof 19 18
END
