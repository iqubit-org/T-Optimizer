# Original circuit
#   qubits: 15
#   qubits used: 15
#   H: 20
#   cnot: 90
#   X: 0
#   T: 70
#   P: 0
#   Z: 0
#   tdepth (by partitions): 30
#   depth  (by critical paths): 73
#   tdepth (by critical paths): 21
# Optimized circuit
#   qubits: 15
#   qubits used: 15
#   H: 20
#   cnot: 379
#   X: 0
#   T: 62
#   P: 1
#   Z: 0
#   tdepth (by partitions): 8
#   depth  (by critical paths): 210
#   tdepth (by critical paths): 8
#   Time: 0.001 s
.v cin a0 b0 c0 d0 a1 b1 c1 d1 e1 a2 b2 c2 d2 e2
.i cin a0 b0 a1 b1 a2 b2
.o cin a0 b0 c0 d0 a1 b1 c1 d1 e1 a2 b2 c2 d2 e2

BEGIN
H c0
H c1
H c2
tof c0 a0
tof a2 d0
tof b2 d0
tof b1 e2
tof a1 e2
tof c0 d2
tof c0 e1
tof c0 d1
tof b0 d2
tof b0 d1
tof a0 d2
tof a0 e1
tof cin d2
T e1
T b0
T c0
T* a0
T* d1
T a1
T b1
T c1
T e2
T a2
T b2
T d0
P cin
T cin
T* d2
tof c0 e1
tof a0 e1
tof c0 b0
tof d1 b0
tof a2 b0
tof b2 b0
tof c2 b0
tof d1 c0
tof a2 c0
tof b2 c0
tof c2 c0
tof a2 d1
tof b2 d1
tof c2 d1
tof c1 a1
tof c1 b1
tof b2 a2
tof d0 a2
tof c2 b2
tof d2 cin
tof cin d2
tof b2 d0
tof a2 d0
tof b1 e2
tof b1 c1
tof a1 e2
tof a1 c1
tof d1 cin
tof a0 cin
tof c0 a0
tof c0 d1
tof b0 cin
tof b0 a0
tof b0 c0
T a0
T* a1
T* b1
T c1
T c2
T* d0
T* b2
T b0
tof c0 a0
tof d1 a0
tof c1 a1
tof c1 b1
tof d0 c2
tof b2 d0
tof b2 b0
tof d0 b2
tof c2 b0
tof c2 b2
tof b1 c1
tof a1 c1
H d1
H c1
H b2
tof b0 c1
tof b2 e2
tof c1 b0
tof b1 a1
tof d1 a2
tof a0 c0
H c1
H b0
tof c0 a0
tof d1 c0
tof a1 c0
tof a2 d1
tof a1 d1
tof c1 d1
tof a1 a2
tof c1 a2
tof a1 b1
tof c1 a1
tof b2 c1
tof e2 c1
tof b2 e1
tof e2 e1
tof cin b2
tof e2 b2
tof cin e2
tof b2 e2
tof b2 cin
tof c2 d0
tof c1 e1
tof c1 b2
tof c1 cin
tof a1 cin
tof a1 c1
tof a2 c1
tof d1 c1
tof d1 a2
tof c0 cin
tof c0 e1
tof c0 a1
T cin
T b2
T* c0
T* a1
T* d1
T e1
T a2
T b0
T* c1
tof c0 cin
tof b0 cin
tof c0 b2
tof a1 b2
tof b0 b2
tof a2 c0
tof b0 c0
tof a0 c0
tof a2 a1
tof a0 a1
tof b0 d1
tof a2 e1
tof b0 e1
tof a0 e1
tof b0 c1
tof a2 c1
tof a2 b0
tof a1 e1
tof a1 d1
tof c0 d1
tof c0 a1
tof cin c1
tof cin b2
tof cin e1
tof cin a1
tof cin c0
tof cin a2
T* b0
T* cin
T a2
tof a2 b0
tof b1 b0
tof a2 cin
tof b1 a2
tof cin a2
tof b0 a2
tof e2 c1
tof cin d1
tof d0 b2
tof c2 b2
tof c2 d0
tof a0 cin
tof c0 cin
tof c0 a0
H a1
H a2
H e2
H c1
tof a1 d2
tof a2 d2
tof cin a0
tof d1 a0
tof a2 a0
tof d0 a0
tof d1 cin
tof a2 cin
tof d0 cin
tof a2 d1
tof d0 d1
tof d0 a1
tof b2 a1
tof a2 e1
tof d0 e1
tof b2 e1
tof d0 c2
tof e2 d0
tof b2 d0
tof a1 e1
tof a1 b2
tof a1 d0
tof a0 d1
tof a0 cin
tof c0 cin
T a1
T e2
T* b2
T* d0
T a2
T* a0
T e1
T* d2
tof e2 a1
tof d0 a1
tof b2 e2
tof e1 b2
tof c1 b2
tof c1 d0
tof c1 a2
tof d2 e1
tof c1 e1
tof c1 d2
tof e1 d2
tof a2 a0
tof a2 e1
tof d0 a0
tof d0 d2
tof d0 e1
tof b2 d2
tof b2 e1
tof e2 a0
tof e2 e1
tof e2 d0
T* d2
T e2
T c1
T* a2
T* d0
T b2
tof e2 d2
tof c1 e2
tof d0 e2
tof d0 c1
tof c2 c1
tof d0 a2
tof c2 a2
tof c2 d0
tof c2 b2
tof d0 b2
tof a2 b2
tof c1 b2
tof c1 d0
tof c1 a2
tof d2 e2
tof d1 b1
tof d2 c2
tof a2 a1
tof d2 a1
tof d1 a0
tof c2 a1
tof c1 a1
tof c1 c2
tof b1 b2
tof b0 b2
tof b0 b1
tof cin b1
tof cin d1
tof c0 b1
tof c0 d1
tof c0 cin
H e2
H d0
tof b1 d1
tof b2 d1
tof b1 a0
tof d2 a0
tof b2 a0
tof a2 a0
tof a1 c2
tof d0 e2
tof c2 a1
tof a2 a1
tof a2 c2
tof d2 a1
tof d2 c2
tof b1 b2
tof d1 a0
tof d1 b2
tof cin a0
tof cin b1
tof cin c2
tof cin b2
tof cin d1
H d0
H a1
H b1
tof e2 e1
tof d1 cin
tof b2 cin
tof e2 cin
tof d0 cin
tof b2 d1
tof c2 d1
tof e2 d1
tof e2 b2
tof d0 b2
tof d0 c2
tof a2 d2
tof a0 d2
tof e2 d2
tof d0 d2
tof a0 a2
tof e2 a2
tof d0 a2
tof e2 a0
tof d0 a0
tof a2 a0
tof d2 a2
tof d1 c2
tof cin b2
tof cin c2
tof cin d1
tof c0 c2
tof e1 d1
tof e1 d2
T e2
T d0
T* e1
T* d2
T* cin
T d1
T a1
tof e1 e2
tof a1 e2
tof d2 d0
tof a2 d0
tof a0 d0
tof d2 e1
tof cin e1
tof c0 e1
tof c2 e1
tof a2 e1
tof a0 e1
tof b1 e1
tof cin d2
tof a1 d2
tof c0 d2
tof c2 d2
tof a2 d2
tof a0 d2
tof b1 d2
tof a1 cin
tof c0 cin
tof c2 cin
tof b1 cin
tof a1 d1
tof b1 a1
tof c2 c0
tof b1 c0
tof a0 a2
tof cin d1
tof d2 cin
tof e1 d1
tof e1 d2
tof d0 d1
tof d0 cin
tof d0 e1
tof e2 d2
tof e2 a1
tof e2 a2
tof e2 d0
T* e2
T* d0
T a2
T b1
T* a1
T* c0
T e1
tof a2 e2
tof a0 e2
tof a2 d0
tof b1 d0
tof a1 d0
tof b1 a2
tof a1 a2
tof a0 a2
tof a1 b1
tof c2 c0
tof c2 e1
tof c0 e1
tof a1 c0
tof b1 e1
tof b1 c0
tof b1 a1
tof e2 a2
tof a0 d0
tof e2 d0
tof b2 e1
tof c2 d1
tof c2 b2
tof c0 d1
tof c0 b2
tof c0 c2
H cin
H a2
H a1
tof a2 e1
tof a1 d1
tof cin d2
tof a0 d0
tof a0 a2
tof e2 d0
tof e2 a2
tof b2 e1
tof b2 d0
tof c2 e1
tof c2 d0
tof b1 c2
tof c0 b2
tof b2 c0
tof c0 b2
tof b2 a2
tof a2 b2
tof b2 a2
tof a2 c1
tof c1 a2
tof a2 c1
tof c1 e2
tof e2 c1
tof c1 e2
tof e2 cin
tof cin e2
tof e2 cin
tof cin b1
tof b1 cin
tof cin b1
tof b1 d1
tof d1 b1
tof b1 d1
tof d1 a1
tof a1 d1
tof d1 a1
tof a1 b0
tof b0 a1
tof a1 b0
tof b0 c2
tof c2 b0
tof b0 c2
tof c2 d2
tof d2 c2
tof c2 d2
tof d2 e1
tof e1 d2
tof d2 e1
tof e1 a0
tof a0 e1
tof e1 a0
END
