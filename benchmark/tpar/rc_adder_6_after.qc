# Original circuit
#   qubits: 14
#   qubits used: 14
#   H: 22
#   cnot: 104
#   X: 8
#   T: 77
#   P: 0
#   Z: 0
#   tdepth (by partitions): 33
#   depth  (by critical paths): 116
#   tdepth (by critical paths): 33
# Optimized circuit
#   qubits: 14
#   qubits used: 14
#   H: 22
#   cnot: 165
#   X: 14
#   T: 47
#   P: 15
#   Z: 0
#   tdepth (by partitions): 22
#   depth  (by critical paths): 157
#   tdepth (by critical paths): 22
#   Time: 0.002 s
.v b0 a0 e b1 a1 b2 a2 b3 a3 b4 a4 b5 a5 z
.i b0 a0 b1 a1 b2 a2 b3 a3 b4 a4 b5 a5 z
.o b0 a0 e b1 a1 b2 a2 b3 a3 b4 a4 b5 a5 z

BEGIN
tof a1 e
tof a1 b1
tof a2 b2
tof a3 b3
tof a4 b4
tof a5 b5
H a1
tof a1 b0
tof a1 a0
T a1
T* b0
T* a0
T b5
tof b0 a1
tof a0 a1
tof b0 a0
tof a1 a0
tof a1 b0
T a1
tof b0 a1
tof a0 a1
H a0
tof a2 e
H e
tof e a0
tof e b1
T e
T* a0
T* b1
tof a0 e
tof b1 e
tof b1 a0
tof e a0
tof e b1
T e
tof b1 e
tof a0 e
tof b0 a1
tof b1 e
P* a1
P b1
P* e
tof b0 a1
tof b1 e
H a0
tof a3 a2
H a2
tof a2 a0
tof a2 b2
T a2
T* a0
T* b2
tof a0 a2
tof b2 a2
tof b2 a0
tof a2 a0
tof a2 b2
T a2
tof b2 a2
tof a0 a2
P a1
P b0
P e
P b2
P a2
P b3
P b4
tof a2 b2
tof b2 a2
P* b2
tof a2 b2
H a0
tof a4 a3
H a3
tof a3 a0
tof a3 b3
T a3
T* a0
T* b3
tof a0 a3
tof b3 a3
tof b3 a0
tof a3 a0
tof a3 b3
T a3
tof b3 a3
tof a0 a3
tof b3 a3
P b3
P* a3
tof b3 a3
H a0
tof a5 a4
H a4
tof a4 a0
tof a4 b4
T a4
T* a0
T* b4
tof a0 a4
tof b4 a4
tof b4 a0
tof a4 a0
tof a4 b4
T a4
tof b4 a4
tof a0 a4
tof b4 a4
P b4
P* a4
tof b4 a4
H a0
tof a5 z
tof b2
H z
tof b5 a0
tof a0 b5
T b5
T z
T* a0
tof z b5
tof a0 z
tof b5 a0
T* b5
T* a0
T z
tof z b5
tof z a0
tof a0 z
tof b5 z
tof e
tof a3
tof a4
H z
tof a0 b5
tof b4 a4
tof b3 a3
tof a2 b2
tof b1 e
tof a4
H a0
tof a0 b4
tof a0 a4
tof a4
T a0
T* b4
T a4
tof b4 a0
tof a4 a0
tof a4 b4
tof a0 b4
tof a0 a4
T* a0
tof a4 a0
tof b4 a0
tof a3
H b4
H a4
tof a4 b3
tof a4 a3
tof a3
T a4
T* b3
T a3
tof b3 a4
tof a3 a4
tof a3 b3
tof a4 b3
tof a4 a3
T* a4
tof a3 a4
tof b3 a4
H b3
tof a5 b4
H a3
tof a3 a2
tof a3 b2
tof b2
T a3
T* a2
T b2
tof a2 a3
tof b2 a3
tof b2 a2
tof a3 a2
tof a3 b2
T* a3
tof b2 a3
tof a2 a3
H a2
tof b4 b3
tof e
tof a4
H b2
tof b2 b1
tof b2 e
tof e
T b2
T* b1
T e
tof b1 b2
tof e b2
tof e b1
tof b2 b1
tof b2 e
T* b2
tof e b2
tof b1 b2
H b1
tof a3
H e
tof e a1
tof e b0
tof b3 b4
T e
T* a1
T* b0
tof a1 e
tof b0 e
tof a1 b0
tof e b0
tof e a1
T e
tof a1 e
tof b0 e
tof b3 a2
tof b3 b4
H b0
tof a1 e
tof b1 b0
tof a2 b0
tof b1 b2
tof a2 b2
tof a2 b1
tof a2 a3
tof b3 a4
tof b4 a0
tof a5 b5
tof b2
tof a1 b1
tof b1 a1
tof a1 b1
tof b1 b2
tof b2 b1
tof b1 b2
tof b2 a3
tof a3 b2
tof b2 a3
tof a3 b3
tof b3 a3
tof a3 b3
tof b3 a4
tof a4 b3
tof b3 a4
tof a4 b4
tof b4 a4
tof a4 b4
tof b4 a0
tof a0 b4
tof b4 a0
tof e b0
tof b0 e
tof e b0
END
