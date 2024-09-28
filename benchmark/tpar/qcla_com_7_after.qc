# Original circuit
#   qubits: 24
#   qubits used: 24
#   H: 39
#   cnot: 215
#   X: 15
#   T: 203
#   P: 0
#   Z: 0
#   tdepth (by partitions): 87
#   depth  (by critical paths): 87
#   tdepth (by critical paths): 27
# Optimized circuit
#   qubits: 24
#   qubits used: 24
#   H: 40
#   cnot: 371
#   X: 48
#   T: 94
#   P: 36
#   Z: 0
#   tdepth (by partitions): 12
#   depth  (by critical paths): 154
#   tdepth (by critical paths): 12
#   Time: 0.005 s
.v a0 b0 c0 a1 b1 c1 a2 b2 c2 d2 a3 b3 c3 a4 b4 c4 d4 a5 b5 c5 d5 a6 b6 c6
.i a0 b0 a1 b1 a2 b2 a3 b3 a4 b4 a5 b5 a6 b6
.o a0 b0 c0 a1 b1 c1 a2 b2 c2 d2 a3 b3 c3 a4 b4 c4 d4 a5 b5 c5 d5 a6 b6 c6

BEGIN
tof a0
tof a1
tof a2
tof a4
tof a5
tof a6
H c0
H c1
H c2
H d2
H c3
H c4
H d4
H d5
H c6
tof c0 a0
tof c0 b0
tof c1 a1
tof c1 b1
tof c2 a2
tof b3 a2
tof b2 a2
tof d2 b2
tof c5 d2
tof b3 a3
tof c3 b3
tof b4 a4
tof c4 a4
tof d4 b4
tof d5 a5
tof c6 a6
tof a6 c6
tof a3 c5
tof d2 c5
tof c2 c5
tof c2 a3
tof b2 a3
tof a2 a3
tof a2 c5
tof a0
tof a1
tof a3
tof c5
tof a5
tof c6
tof a6
tof a4
T c0
T* b0
T* b1
T d2
T* b2
P c3
T c3
T* b3
T* b4
T b6
T c2
T a2
T a0
T a1
T a3
T* c5
T a5
T* c6
T a6
T a4
tof a0 b0
tof b0 c0
tof b2 d2
tof a3 d2
tof c5 d2
tof a5 d2
tof b5 d2
tof c1 b1
tof c4 d2
tof d5 d2
tof b3 b2
tof c2 b2
tof a2 b2
tof a3 b2
tof a2 b3
tof a4 b4
tof c4 b4
tof a2 c2
tof c6 b6
tof c5 a3
tof c1 a1
tof a4 a5
tof b5 a5
tof d5 a5
tof a6 c6
tof c4 a4
tof d4 a4
tof d5 b5
tof a2 c5
tof c2 c5
tof c2 a2
tof b6 a6
tof b3 c5
tof b3 c2
tof b2 c5
tof b2 c2
tof b2 b3
tof c0 b0
tof c0 a0
T d4
T* b5
P c6
T c6
T* a6
T c4
T a5
T* c0
T b2
T b4
T* a4
T* b6
T b3
T c2
T d2
tof c4 b5
tof a5 b5
tof a4 b5
tof d5 b5
tof b4 d4
tof b4 b5
tof b6 a6
tof a6 c6
tof b4 a5
tof a4 a5
tof a5 c4
tof a0 c0
tof b0 c0
tof c3 b2
tof a4 b4
tof a3 b3
tof a2 b3
tof a2 c2
tof b4 a4
tof b2 c2
tof a5 d2
tof c4 a5
tof c6 a6
tof b5 c4
tof d4 b5
tof d4 d2
tof d4 a4
tof b2 b3
tof c2 b3
tof a3 b3
tof c2 b2
tof b5 d4
tof c4 d4
tof b4 d4
tof c4 b5
tof c4 d2
tof b5 d2
tof b4 b5
tof d4 b5
tof a3 b2
tof b3 b2
P a0
P b1
P c1
P a3
P c2
P b4
P c4
P* b3
P* d4
P* c0
P* a1
P* b2
P* b5
P* d2
tof c0 a0
tof c2 a3
tof b3 a3
tof b2 a3
tof b3 c2
tof d4 b4
tof b2 c2
tof b5 b4
tof b2 b3
tof b5 d4
tof c4 d2
tof b4 d2
tof c2 b3
tof b4 d4
tof a3 c2
tof a3 b3
tof a0 c0
P d5
P a0
P* a3
P* b3
P* b4
tof c0 a0
tof c4 b4
tof b5 d4
tof b2 c2
tof a1 b1
tof b4 c4
tof a1
tof a3
tof b4
H b0
H a2
H c5
H a5
H a4
tof c1 b0
tof b1 a1
tof c3 a2
tof c3 c5
tof b3 a3
tof d5 a4
tof c4 b4
tof c4 d5
tof a4 d5
tof a4 c4
tof a5 b6
tof b3 c3
tof c5 c3
tof c5 b3
tof b1 c1
tof b0 c1
tof b0 b1
tof b6
tof b0
tof a4
tof c5
T a5
T* b1
T* c4
T* b3
T* a2
T b6
T* b0
T* a4
T* c5
tof b6 a5
tof b0 b1
tof a1 b1
tof a4 c4
tof c5 b3
tof a3 b3
tof b4 c4
tof c5 a2
tof a3 a2
tof c6 b6
tof c3 a2
tof a1 b0
tof b4 a4
tof c1 b0
tof d5 a4
tof a3 c5
tof c3 c5
tof b3 a2
tof c4 a4
tof b3 c5
tof a5 b6
tof b1 b0
tof c1 b1
tof a1 b1
tof c3 b3
tof a3 b3
tof d5 c4
tof b4 c4
P c1
P d5
P c3
P a2
P b1
P c4
P b3
tof b4 c4
tof a3 b3
tof a5 c6
tof b3 a3
tof c4 b4
tof b1 a1
tof c3 a3
tof c3 b3
tof d5 b4
tof d5 c4
tof c1 b1
tof c1 a1
H d2
H b0
H a4
tof c5 b0
tof a6 d2
tof d2 a6
tof b0 c5
tof b0 a2
T* b0
T a2
T a4
T a6
T* d2
tof c5 b0
tof c5 a2
tof d2 a6
tof a5 a4
tof a6 a4
tof a5 c6
tof a6 d2
tof a4 a6
tof a4 c6
tof b0 a2
T* c6
T a6
T* a4
tof a5 c6
tof a4 c6
tof a4 a6
tof a5 a4
tof a6 a4
tof c6 a4
tof c5 a2
P c5
P* a2
tof c6 a5
tof c5 a2
tof c4
H b0
tof a6 b0
tof a6 d2
tof b0 a6
tof b0 d2
T a6
T* d2
T* b0
tof b0 d2
tof d2 a6
tof a6 d2
tof a6 b0
T a6
tof b0 a6
tof d2 a6
tof c6
H a6
H d2
H b0
H a4
tof d2 c5
tof d2 a2
tof d2 b3
tof b0 b6
tof a4 c4
tof a5 c6
tof b0 a5
tof b6 b0
tof b6
tof b3
tof c4
P d2
T d2
T* c5
T* a2
T* b6
T b3
T c4
tof a2 c5
tof c5 d2
tof c3 c5
tof b0 b6
tof b4 c4
tof c3 b3
tof d5 c4
tof c3 a3
tof a4 b4
tof a5 b0
tof a4 d5
tof d2 c5
tof c5 c3
tof d2 a2
tof c5 a3
T d2
T* d5
T a5
T* b0
T* c5
T* a3
T* b4
T* c4
T b6
T* b3
tof b4 d5
tof c4 d5
tof a2 d2
tof c3 d2
tof a4 d5
tof c6 b6
tof b6 a5
tof b6 b0
tof b3 c5
tof c3 a3
tof a3 c5
tof a4 b4
tof a4 c4
tof c3 b3
tof b4 c4
tof a3 b3
tof c5 b3
tof a5 b0
tof a5 b6
tof d5 c4
P a4
tof a5 c6
tof d5 b4
tof c5 a3
tof b1
tof a1
H a2
tof a2 a1
tof a1 c1
tof a1 b1
tof a2 a1
tof c1 a1
tof c1
T* a1
T* b1
T* c1
tof c1 a1
tof a2 b1
tof b1 a1
tof a2 c1
tof b1 c1
tof a1 c1
P a2
tof a1 b1
tof a0
tof a1
tof b1
tof c2
tof b2
tof a3
tof d4
tof b4
H c1
H d2
H b3
H b0
H c4
tof c1 a0
tof c1 c0
tof a2 a1
tof a2 b1
tof b3 c2
tof b3 b2
tof c5 a3
tof c4 b5
tof c4 d4
tof d5 b4
tof b0 d5
tof d2 c5
tof a1 b1
tof d2
tof b0
tof a0
tof a1
tof c2
tof b2
tof b5
tof d4
T c1
T* c0
T d2
T b0
T a0
T a1
T c2
T* b2
T b5
T* d4
tof a0 c0
tof c0 c1
tof b2 d2
tof d4 b0
tof b3 d2
tof c4 b0
tof a2 a1
tof b3 b2
tof b2 c2
tof c5 b2
tof c4 d4
tof d4 b5
tof d5 d4
tof c1 c0
tof c1 a0
T c5
T d2
T d5
T b0
T b3
T* c2
T c4
T* b5
T b2
T d4
T* c1
tof b3 c5
tof c2 c5
tof a3 d2
tof b2 c5
tof b2 d2
tof c4 d5
tof b5 d5
tof b4 b0
tof d4 d5
tof d4 b0
tof c2 b3
tof b5 c4
tof c0 c1
tof a0 c1
tof c4 b5
tof c4 d4
tof b3 c2
tof b3 b2
tof d5 d4
tof c5 b2
tof d2
H c0
H a2
H b2
H c2
H c3
H d4
H b5
H b6
H a4
tof a5
tof a6
H c1
P c1
H c1
P c1
H c1
P c1
tof c1 a2
tof a2 c1
tof c1 a2
tof a2 c5
tof c5 a2
tof a2 c5
tof c5 b6
tof b6 c5
tof c5 b6
tof b6 c6
tof c6 b6
tof b6 c6
tof c6 a6
tof a6 c6
tof c6 a6
tof a6 a5
tof a5 a6
tof a6 a5
tof a5 b0
tof b0 a5
tof a5 b0
tof b0 a0
tof a0 b0
tof b0 a0
tof c2 d2
tof d2 c2
tof c2 b2
tof b2 c2
tof c2 b2
tof b2 b3
tof b3 b2
tof b2 b3
tof b3 a3
tof a3 b3
tof b3 a3
tof d5 a4
tof a4 d5
tof d5 a4
tof c4 d4
tof d4 c4
tof c4 d4
tof d4 b5
tof b5 d4
tof d4 b5
tof b5 b4
tof b4 b5
tof b5 b4
END
