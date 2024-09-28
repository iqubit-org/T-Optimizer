# Original circuit
#   qubits: 30
#   qubits used: 30
#   H: 28
#   cnot: 196
#   X: 28
#   T: 196
#   P: 0
#   Z: 0
#   tdepth (by partitions): 84
#   depth  (by critical paths): 62
#   tdepth (by critical paths): 18
# Optimized circuit
#   qubits: 30
#   qubits used: 30
#   H: 28
#   cnot: 366
#   X: 27
#   T: 84
#   P: 34
#   Z: 0
#   tdepth (by partitions): 6
#   depth  (by critical paths): 153
#   tdepth (by critical paths): 6
#   Time: 0.004 s
.v c00 c01 c10 c11 e10 e11 c20 c21 c30 c31 a0 a1 e30 e31 c40 c41 c50 c51 e50 e51 c60 c61 c70 c71 b0 b1 b2 b3 e70 e71
.i c00 c01 c10 c11 e10 e11 c20 c21 c30 c31 e30 e31 c40 c41 c50 c51 e50 e51 c60 c61 c70 c71 e70 e71
.o c00 c01 c10 c11 e10 e11 c20 c21 c30 c31 a0 a1 e30 e31 c40 c41 c50 c51 e50 e51 c60 c61 c70 c71 b0 b1 b2 b3 e70 e71

BEGIN
H e10
H e11
H a0
H a1
H e50
H e51
tof c21
tof c61
H b0
H b1
tof c11 c00
tof c11 c01
tof e10 c11
tof e11 c10
tof e11 e10
tof c31 c20
tof c31 c21
tof c31 c30
tof a1 c31
tof c50 c41
tof c51 c40
tof e50 c50
tof c51 c50
tof e51 c51
tof b0 c60
tof b1 c61
tof c70 c60
tof c70 c61
tof b1 c71
tof c71 c70
tof b1 b0
tof b3 b2
tof c70 b0
tof c61 b1
tof c60 b0
tof c61 c70
tof c50 e50
tof c41 e50
tof c40 c50
tof c41 e51
tof c21 c30
tof e10 e11
tof e11 b2
tof e10 b3
tof c11 b3
tof c10 b3
tof b3 b2
tof c10 e10
tof c00 b2
tof c00 b3
tof c00 e11
tof c01
tof c21
tof e50
tof c70
tof c30
tof c41
tof e51
tof b1
tof c61
T* c00
T* c11
T e11
T* c10
T* c20
T* c40
T* b0
T* e10
T* c01
T* c21
T* c31
T* e50
T* c51
T* c70
T* c71
T c30
T c41
T* e51
T b1
T* c61
T b3
T* b2
T* c50
T* c60
tof e11 c00
tof e10 c10
tof c01 c00
tof c10 c00
tof c10 c11
tof c01 c10
tof c30 e11
tof b3 e11
tof a1 e11
tof c10 e11
tof e10 c10
tof c51 c40
tof c41 c40
tof e51 c40
tof c50 c40
tof a0 c20
tof c01 e10
tof c21 e10
tof c70 b0
tof c71 b0
tof c61 b0
tof b3 e10
tof b2 e10
tof c60 b0
tof a1 e10
tof c31 c21
tof c30 c21
tof a0 c21
tof a1 c21
tof c41 e50
tof c50 e50
tof a0 c31
tof a1 c31
tof b1 c70
tof c50 c51
tof a1 c30
tof c50 e51
tof e51 c41
tof c60 b1
tof c61 b1
tof a0 a1
tof b3 b2
tof c30 b2
tof c30 b3
tof c71 b1
tof c51 c41
tof c70 c60
tof c70 c71
tof c70 c61
tof c31 b2
tof e50 c41
tof c31 a0
tof e50 c50
tof e50 e51
tof e50 c51
tof c31 c30
tof c21 b2
tof c01 b2
tof c01 b3
tof c21 a0
tof c21 a1
tof c21 c30
tof e10 b2
tof b0 b1
tof b0 c70
tof c40 c41
tof c20 e10
tof c20 a0
tof c40 c50
tof c10 b2
tof e11 b3
tof e11 e10
tof c10 c01
tof c11 c01
tof c00 c01
tof c00 b2
tof c11 b3
tof c11 e10
tof c01 e11
tof c11 c10
T* c31
T c20
T* a1
T* c50
T e50
T* c51
T* b0
T c60
T* c71
T* c11
T c00
T* c21
T c30
T* c40
T e51
T* c70
T c61
T c10
T* c01
T* e11
T a0
T* e10
T c41
T b1
tof c21 a1
tof c30 a1
tof c20 c31
tof a1 c20
tof e50 c50
tof c21 a1
tof a0 a1
tof c40 e50
tof e51 c51
tof c51 e50
tof c41 c51
tof c70 c71
tof c61 c71
tof c60 b0
tof c71 c60
tof c70 c71
tof c10 c11
tof c01 c11
tof b1 c71
tof c10 c00
tof c01 c00
tof a0 c21
tof c41 c40
tof b1 c70
tof c01 c10
tof a0 e11
tof a0 e10
tof c10 c01
tof c30 e11
tof c21 e10
tof c00 c01
tof c71 c61
tof c11 c01
tof c11 c00
tof c60 c71
tof b0 c60
tof b0 c61
tof b0 c70
tof b0 b1
tof c51 e51
tof e50 c51
tof a1 e11
tof c50 e51
tof c50 c40
tof a1 c30
tof c50 c41
tof c50 e50
tof c20 a1
tof c31 c20
tof c20 e11
tof c31 e10
tof c31 c30
tof c31 c21
tof c31 a0
P c11
P b2
P c01
P c10
P a1
P c21
P a0
P c30
P c51
P c40
P c41
P e51
P c71
P c70
P b1
P c61
H b2
H c10
H c21
H c30
tof c50
H c40
H e51
tof b2
tof c40
H c70
H c61
H e30
H e31
H e10
H e11
tof c21 b2
tof c30 c10
tof c30 c21
tof e31 c30
tof e31 e30
tof c70 e51
tof e10 e51
tof e11 c40
tof e10 c61
tof c61 c70
tof e11 c61
tof e11 e10
tof c61 e10
tof c70 c61
tof e10 e11
tof c40 e11
tof c40 e10
tof c40 c61
tof e51 c70
tof e30 e31
tof c30 e30
tof c21 e30
tof b2 c21
tof c10 c30
tof b2 e31
tof b2 e30
tof c21
tof e11
tof c40
tof e10
tof b2
tof e30
tof e31
tof c61
T* c10
T* c30
T* c70
T* c21
T* e11
T* c40
T e10
T b2
T e30
T* e31
T c61
T* e51
tof c21 c30
tof e30 b2
tof e31 b2
tof c30 c10
tof b2 c30
tof e11 c70
tof e10 c70
tof b2 c21
tof c61 e11
tof e51 e11
tof c40 e11
tof e10 c40
tof e30 b2
tof e31 b2
tof e51 e10
tof e51 c61
tof e10 e51
tof e10 c61
tof b2 e31
tof b2 e30
tof e11 e51
tof c70 c61
tof c21 b2
tof c70 e10
tof c70 e11
tof c70 c40
tof c30 b2
tof c10 b2
tof b2 e31
tof b2 e30
tof c10 c21
T c30
T* c40
T c70
T* c10
T c21
T e11
T* e10
T b2
T e30
T* e31
T c61
T e51
tof c70 c40
tof c21 c30
tof b2 c30
tof e30 c30
tof e31 c30
tof c21 c10
tof e11 c70
tof e10 c70
tof c61 c70
tof e10 e11
tof c61 e11
tof e30 b2
tof e31 b2
tof c61 e10
tof e51 e10
tof e31 e30
tof e30 e31
tof b2 e31
tof e11 e51
tof e11 c61
tof c10 c21
tof c40 c70
tof c30 c21
tof e30 e31
tof c10 c21
P e30
P* e31
P c10
P* c21
P c61
P c40
tof e31 e30
tof c21 c10
tof e51 c61
tof c70 c40
tof c10 c21
tof e30 e31
P e30
P c10
P e51
P* c61
P c70
P* c40
tof c61 e51
tof c40 c70
tof c70 c40
tof e51 c61
H e30
H c10
H c61
H c40
tof e30
H e70
H e71
tof c61 c10
tof c40 e30
tof e71 e30
tof e71 c61
tof c40 c61
tof c61 e71
tof e30 c40
tof c10 e71
tof c10 e70
tof e30 c61
tof c40
tof e30
tof c61
T* e71
T* c40
T e30
T c61
T c10
T* e70
tof c10 c61
tof c61 e71
tof c61 c40
tof c61 e30
tof e70 c61
tof e70 c10
tof c61 e70
tof e30 c10
tof c40 c10
tof e71 c10
tof c10 e70
T* c40
T e30
T* e71
T c10
T* c61
T e70
tof e71 e30
tof e70 c10
tof e30 c40
tof c10 e30
tof e70 c61
tof c10 e70
tof c10 c61
tof c40 e70
tof c40 e71
tof c40 c10
P c61
P e70
P e71
tof e70 c61
tof e30 e71
P* c61
P e30
P* e71
tof e71 e30
tof e70 c61
tof e30 e71
H e70
H e71
tof c60
tof e31 c10
tof c10 e31
tof e31 c10
tof c10 c01
tof c01 c10
tof c10 c01
tof c01 c00
tof c00 c01
tof c01 c00
tof c00 b3
tof b3 c00
tof c00 b3
tof b3 e30
tof e30 b3
tof b3 e30
tof e30 c40
tof c40 e30
tof e30 c40
tof c40 c50
tof c50 c40
tof c40 c50
tof c50 c41
tof c41 c50
tof c50 c41
tof c41 e50
tof e50 c41
tof c41 e50
tof e50 e10
tof e10 e50
tof e50 e10
tof e10 c30
tof c30 e10
tof e10 c30
tof c30 a0
tof a0 c30
tof c30 a0
tof c21 c20
tof c20 c21
tof c21 c20
tof c20 c31
tof c31 c20
tof c20 c31
tof c31 a1
tof a1 c31
tof c31 a1
tof e51 e11
tof e11 e51
tof e51 e11
tof e11 b2
tof b2 e11
tof e11 b2
tof b2 c61
tof c61 b2
tof b2 c61
tof c61 c60
tof c60 c61
tof c61 c60
tof c60 b0
tof b0 c60
tof c60 b0
tof c70 b1
tof b1 c70
tof c70 b1
END
