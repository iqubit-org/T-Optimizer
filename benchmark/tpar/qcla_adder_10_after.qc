# Original circuit
#   qubits: 36
#   qubits used: 36
#   H: 50
#   cnot: 267
#   X: 0
#   T: 238
#   P: 0
#   Z: 0
#   tdepth (by partitions): 102
#   depth  (by critical paths): 79
#   tdepth (by critical paths): 24
# Optimized circuit
#   qubits: 36
#   qubits used: 36
#   H: 50
#   cnot: 648
#   X: 0
#   T: 162
#   P: 22
#   Z: 1
#   tdepth (by partitions): 13
#   depth  (by critical paths): 195
#   tdepth (by critical paths): 13
#   Time: 0.023 s
.v a0 b0 c0 a1 b1 c1 a2 b2 c2 d2 a3 b3 c3 a4 b4 c4 d4 a5 b5 c5 d5 a6 b6 c6 d6 a7 b7 c7 a8 b8 c8 d8 a9 b9 c9 e
.i a0 b0 a1 b1 a2 b2 a3 b3 a4 b4 a5 b5 a6 b6 a7 b7 a8 b8 a9 b9
.o a0 b0 c0 a1 b1 c1 a2 b2 c2 d2 a3 b3 c3 a4 b4 c4 d4 a5 b5 c5 d5 a6 b6 c6 d6 a7 b7 c7 a8 b8 c8 d8 a9 b9 c9 e

BEGIN
H c0
H c1
H c2
H d2
H c3
H c4
H d4
H c5
H d5
H c6
H d6
H c7
H c8
H d8
H c9
tof e c0
tof d2 a2
tof b3 a3
tof d5 a5
tof b7 a7
tof c7 b7
tof b9 a9
tof b9 c9
tof a9 c9
tof c8 a9
tof b8 a9
tof a8 a9
tof b8 d8
tof b7 c7
tof b6 d6
tof c6 a7
tof b5 d5
tof a5 d5
tof c4 d5
tof d2 a3
tof c2 a3
tof b2 a3
tof a2 a3
tof b2 d2
tof b1 c1
tof c0 e
tof b0 e
T a0
T b0
T c0
T* e
T a1
T b1
T* c1
T b2
T* a2
T* d2
T b3
P c3
T c3
T a4
T b4
T d4
T b5
T* a5
T a6
T b6
T* d6
T c7
T* b7
T a8
T b8
T* d8
T b9
T* c9
T c2
T a3
T c4
T* d5
T c6
T* a7
T c8
T a9
tof b0 a0
tof e b0
tof b2 c0
tof d2 c0
tof c1 a1
tof b1 a1
tof a2 b2
tof b3 a2
tof c3 a2
tof c2 a2
tof a3 a2
tof d2 a2
tof a3 b3
tof c3 b3
tof d4 a4
tof b4 a4
tof d5 b5
tof c4 d4
tof c4 b5
tof b6 a6
tof d6 a6
tof c4 a5
tof d5 a5
tof d6 b6
tof a7 c6
tof c6 c7
tof c6 b7
tof b9 b8
tof c8 b8
tof a9 b8
tof d8 a8
tof d8 b8
tof c8 b9
tof a9 b9
tof b9 d8
tof a3 c2
tof d5 c4
tof c6 a7
tof c4 d5
tof b9 a9
tof b9 c8
tof b8 b9
tof a8 b8
tof b8 a9
tof b8 c8
tof b8 d8
tof b8 c9
tof d6 a7
tof b6 d6
tof a7 c6
tof a6 c6
tof a6 a7
tof a6 d6
tof d4 d5
tof d4 c4
tof b4 d5
tof a4 b4
tof a4 c4
tof a4 d4
tof b3 a3
tof b3 c2
tof a2 c3
tof a2 c2
tof a2 b3
tof b2 d2
tof b1 c1
tof b0 c0
tof a0 b0
tof a0 e
T* a0
T* e
T b0
T* a1
T d2
T c0
T b2
T c3
T* a2
T* a3
T* b4
T* d5
T a4
T b5
T* a5
T b6
T* a6
T d6
T c7
T* b7
T b9
T* d8
T a8
T b8
P c9
T c9
T* b3
T* c2
T* d4
T c4
T* a7
T c6
T* c8
T* a9
T c5
tof b0 e
tof e a0
tof c1 a1
tof b2 c0
tof a2 c0
tof a3 c0
tof b3 c0
tof c2 c0
tof a2 b2
tof a3 b2
tof b3 b2
tof c2 b2
tof a3 a2
tof a4 b4
tof a4 d5
tof c4 b5
tof d4 a4
tof d4 b5
tof a6 b6
tof d6 a6
tof b8 b9
tof c9 b9
tof c8 b9
tof a9 b9
tof a8 b9
tof a8 d8
tof c9 b8
tof c8 b8
tof a9 b8
tof b8 a8
tof c2 b3
tof c8 c9
tof a9 c9
tof c6 a7
tof c8 a9
tof a7 c6
tof d4 c4
tof c9 a9
tof b8 c8
tof b8 c9
tof d8 b8
tof b9 b8
tof b9 a8
tof b9 c8
tof c7 b7
tof c7 a7
tof a6 d6
tof a6 c6
tof b6 d6
tof b6 c6
tof b5 a5
tof d5 d4
tof b4 d4
tof d4 a4
tof a2 a3
tof a2 b3
tof a2 c2
tof c3 a3
tof c3 a2
tof d2 b3
tof e b0
tof a0 b0
T* b9
tof a9 b9
tof b3 d2
tof a2 c3
tof b4 d5
tof a6 b6
tof a7 c7
tof d8 b8
tof b9 c9
tof c9 c0
tof b8 c0
P c1
P d2
P c3
P d5
P a5
P c4
P b6
Z b7
P c7
P b8
P c9
P* c0
tof b3 d2
tof a2 c3
tof b4 d5
tof b5 c4
tof a6 b6
tof a7 c7
tof d8 b8
tof d8 c0
tof b9 c9
tof b9 c0
tof b1 a1
tof c9 c0
tof c9 b9
tof b8 c0
tof b8 d8
tof c7 a7
tof c3 a2
H b0
H c2
H b2
H d4
H a4
H c6
H d6
H c8
H a8
tof c1 b0
tof c1 a1
tof c0 c1
tof a2 b2
tof a2 c3
tof c5 d4
tof c6 d4
tof c6 c5
tof a7 c7
tof b7 d6
tof d8 b8
tof a9 a8
tof b9 c9
tof a8 a9
tof a8 b9
tof b8 d8
tof d6 b7
tof d6 a7
tof b2 a2
tof c2 a3
tof c1 c0
tof a1 c0
tof b0 a1
tof b0 c1
tof b0 c0
T c1
T* a1
T* b0
T c0
T a2
T* b2
T* d6
T a7
T* a8
T b9
T* c5
T d4
T c2
T* a3
tof b0 c1
tof c0 a1
tof a4 c1
tof a1 c1
tof b1 c0
tof c0 b0
tof a3 c2
tof c2 a2
tof c2 b2
tof b7 d6
tof a7 d6
tof c6 a7
tof a9 a8
tof b9 a8
tof c8 b9
tof d4 c5
tof c6 d4
tof a3 c2
tof a5 c4
tof a4 c4
tof b5 c4
tof a4 a5
tof a5 a4
tof c4 a4
tof c2 a3
tof a8 b9
tof d6 a7
tof a2 a3
tof b0 c0
tof a1 c0
T* a2
T b2
T c1
T* a4
T* a5
T c4
T b7
T* d6
T a9
T* a8
P c5
T c5
P c6
T c6
T* d4
T* a7
T c8
T* b9
tof c3 b2
tof b2 a2
tof a3 b2
tof a5 c1
tof c4 a4
tof a4 c1
tof b5 c4
tof c4 a5
tof c7 d6
tof d6 b7
tof c9 a8
tof a8 a9
tof d4 c5
tof a7 c6
tof b9 c8
tof c8 b9
tof c6 a7
tof c5 d4
tof a9 a8
tof b7 d6
tof a5 c4
tof a4 c4
tof a2 b2
tof a2 d2
tof c3 d2
tof b3 d2
tof a5 d5
tof b5 d5
tof b4 d5
tof b7 b6
tof c7 b6
tof a6 b6
tof a7 d4
P* d2
P* d5
P* b6
P* d4
tof b3 d2
tof a2 d2
tof c3 d2
tof b4 d5
tof a5 d5
tof b5 d5
tof a6 b6
tof b7 b6
tof c7 b6
tof a7 d4
tof b1 b0
tof d8 b8
tof a5 b5
tof b0 b1
tof d2 b3
H a1
H c5
H a4
tof b3 a1
tof b3 d2
tof c2 b3
tof b2 c2
tof c6 c5
tof a7 a4
tof c6 a7
tof c1 c6
tof b8 d8
tof d8 b8
tof c6 c1
tof a7 c6
tof a4 c6
tof a4 a7
tof c5 c1
tof b3 b2
tof b3 c2
tof a1 b2
tof a1 c2
tof a1 b3
T* b3
T* c2
T b2
T c6
T* a4
T* a7
T c1
T* c5
T* a1
tof b2 b3
tof a1 c2
tof d2 c2
tof c2 b3
tof a7 c6
tof c1 c6
tof c5 c6
tof a1 b2
tof d2 b2
tof c1 a4
tof c5 a4
tof d2 a1
tof c1 c5
tof a7 c5
tof a4 a7
tof c2 b2
tof b3 c2
tof b3 a1
T a4
tof a7 a4
tof c5 a4
P a1
tof b8 d8
tof b7 c7
tof b6 a6
tof d5 b4
tof b3 d2
H c2
H c5
tof b4 d5
tof c6 c2
tof c6 b4
tof c6 d4
tof c6 c1
tof a6 b6
tof c7 b7
tof c6 a4
tof a4 c6
tof b7 c7
tof b6 a6
tof c1 c6
tof c2 d4
tof c2 c6
tof c2 c1
tof c2 a4
tof c2 b4
P c2
T c2
T* c1
T* a4
T c6
T c5
T* d4
T* b4
tof a4 c1
tof c6 c1
tof c1 c2
tof c5 c1
tof c5 a4
tof c5 c6
tof d5 b4
tof b4 c5
tof b4 d4
tof d4 b4
tof c5 d4
tof a4 c6
tof c1 a4
tof c1 c5
tof c1 b4
T* c1
T* b4
T d4
tof d4 c1
tof d4 b4
tof b4 d4
tof c1 d4
tof a6 b6
tof c5 d5
H a4
H d4
H a3
H c4
H d6
H a8
H c6
tof a3 a1
tof d2 b3
tof c2 a3
tof d5 c5
tof c4 b4
tof c4 d5
tof c6 c1
tof a7 c1
tof b6 a6
tof b9 a4
tof a8 b9
tof d8 b8
tof b9 d8
tof c8 b9
tof a8 b9
tof d8 a8
tof a4 a8
tof a4 b9
tof a4 d8
tof a4 c8
tof d4 b6
tof a3 c2
tof d2 c2
tof a1 d2
T a3
T* a1
T* c2
T d2
T c4
T* b4
T* d5
T d4
T* b6
T* a8
T d8
T* a4
T* b9
T c8
T c1
tof c2 a3
tof d2 a1
tof b6 a3
tof d6 a3
tof a1 a3
tof b3 a1
tof d2 c2
tof b3 d2
tof d5 b4
tof b4 c4
tof d6 d4
tof d8 a8
tof a4 a8
tof b9 a8
tof c8 a8
tof b9 d8
tof c8 a4
tof a4 d8
tof c8 b9
tof a7 c1
tof c6 a7
tof b9 c8
tof a4 c8
tof a8 c8
tof d4 b6
tof c4 b4
tof c4 d5
tof c2 d2
tof a1 d2
T c4
T d6
T* d4
T* b6
T a3
T c8
T* a8
T* d8
T c6
T* c1
T* a7
tof b6 d6
tof c5 c4
tof a6 d6
tof d5 c4
tof b4 c4
tof a6 b6
tof b6 d4
tof d8 c8
tof a6 a3
tof b8 c8
tof b8 d8
tof d8 a8
tof c1 c6
tof a7 c1
tof c6 c1
tof c1 a7
tof c8 d8
tof c8 a8
tof b6 a3
tof d4 a3
tof d6 a3
tof d6 b6
tof d6 d4
P a8
tof a9 c9
tof c8 b8
tof d6 a6
tof c4 c5
tof a2 c3
tof a1 b3
H d2
H b4
H b6
H d8
H a4
H a7
H b2
H c6
H c1
H b9
tof b3 a1
tof b2 b3
tof c3 a2
tof a3 c3
tof c5 c4
tof c6 c5
tof b5 a5
tof a6 d6
tof c9 b9
tof b8 c8
tof b9 b8
tof c9 b9
tof c9 a9
tof a9 c9
tof d6 a6
tof c6 b5
tof c3 a3
tof b2 a3
tof b3 c3
tof b1 b0
T b2
T* b3
T* a3
T c3
T c6
T* c5
T* b5
T b8
tof a3 b2
tof c3 b3
tof a6 b3
tof d6 b3
tof c7 b3
tof b7 b3
tof c1 b3
tof b3 b2
tof a2 b3
tof a1 c3
tof a6 c3
tof d6 c3
tof c7 c3
tof b7 c3
tof c1 c3
tof c3 a3
tof a2 c3
tof b5 c6
tof c5 c6
tof a5 c5
tof a5 b5
tof c9 b8
tof a9 b8
tof c1 a6
tof d6 a6
tof c1 c7
tof b7 c7
tof a9 c9
tof b9 c9
tof c6 c5
tof c5 b5
tof a3 c3
tof b2 b3
tof b3 c3
tof b2 a3
T c6
T c1
T* a6
T* c7
T b2
T b9
T* b8
T* c9
tof a6 c1
tof c4 c6
tof c5 c6
tof a5 c6
tof d6 c1
tof b5 c6
tof d6 a6
tof b7 c7
tof c7 a6
tof b8 b9
tof d6 b2
tof b7 b2
tof c8 b9
tof c8 b8
tof a9 c9
tof c9 b8
tof b8 c9
tof b9 c9
tof b9 b8
tof c7 b2
tof a6 c7
tof c1 b2
tof c1 c7
tof c1 a6
tof c0 b0
tof b0 b1
tof c2 a1
tof a1 a3
tof d2 a2
tof a2 b3
tof d5 c4
tof c4 c6
tof b4 a5
tof a5 c5
tof d4 d6
tof d6 c1
tof b6 b7
tof b7 a6
tof a8 c8
tof c8 b9
tof c9 a9
tof a9 b8
tof a9 c9
tof b8 a9
tof c8 a8
tof b9 c8
tof b7 b6
tof a6 b7
tof d6 d4
tof c1 d6
tof a5 b4
tof c5 a5
tof c4 d5
tof c6 c4
tof a2 d2
tof b3 a2
tof a1 c2
tof a3 a1
tof b0 c0
tof b1 b0
tof e b2
tof a0 b2
H c3
H b5
H c7
H b8
tof b0 e
tof e b0
tof b0 e
tof e b2
tof b2 e
tof e b2
tof b2 c2
tof c2 b2
tof b2 c2
tof c2 c3
tof c3 c2
tof c2 c3
tof c3 c6
tof c6 c3
tof c3 c6
tof c6 c7
tof c7 c6
tof c6 c7
tof c7 b9
tof b9 c7
tof c7 b9
tof b9 c9
tof c9 b9
tof b9 c9
tof c9 a4
tof a4 c9
tof c9 a4
tof a4 c4
tof c4 a4
tof a4 c4
tof c4 b5
tof b5 c4
tof c4 b5
tof b5 b4
tof b4 b5
tof b5 b4
tof b4 d5
tof d5 b4
tof b4 d5
tof d5 c1
tof c1 d5
tof d5 c1
tof c1 a3
tof a3 c1
tof c1 a3
tof a3 a2
tof a2 a3
tof a3 a2
tof a2 a1
tof a1 a2
tof a2 a1
tof d6 a6
tof a6 d6
tof d6 a6
tof b7 b6
tof b6 b7
tof b7 b6
tof b6 d4
tof d4 b6
tof b6 d4
tof d4 c5
tof c5 d4
tof d4 c5
tof c5 a7
tof a7 c5
tof c5 a7
tof c8 b8
tof b8 c8
tof c8 b8
tof b8 a8
tof a8 b8
tof b8 a8
tof c0 b1
tof b1 c0
tof c0 b1
tof d2 b3
tof b3 d2
tof d2 b3
END
