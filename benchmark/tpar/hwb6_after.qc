# Original circuit
#   qubits: 7
#   qubits used: 7
#   H: 30
#   cnot: 131
#   X: 8
#   T: 105
#   P: 0
#   Z: 0
#   tdepth (by partitions): 45
#   depth  (by critical paths): 166
#   tdepth (by critical paths): 45
# Optimized circuit
#   qubits: 7
#   qubits used: 7
#   H: 26
#   cnot: 270
#   X: 33
#   T: 75
#   P: 13
#   Z: 1
#   tdepth (by partitions): 24
#   depth  (by critical paths): 248
#   tdepth (by critical paths): 24
#   Time: 0.002 s
.v a b c d e f 0
.i a b c d e f
.o a b c d e f 0

BEGIN
tof c a
tof a
H c
tof b a
tof c b
tof 0 c
tof f d
tof c 0
tof b 0
tof a
T 0
T c
T* b
T a
T d
tof b 0
tof a 0
tof a c
tof a b
tof 0 c
tof 0 a
T b
T* 0
tof a b
tof a 0
tof b 0
H a
H 0
tof a b
tof a c
tof 0 a
tof b
P c
T c
T 0
T* a
P b
T b
tof 0 c
tof a c
tof b 0
tof b a
tof 0 b
T a
T* 0
tof b a
tof b 0
tof a 0
tof e f
tof a
H b
H f
tof f a
tof d b
tof f 0
tof f e
tof c f
tof f c
tof d e
tof 0 c
tof b e
tof b d
tof b f
tof a e
tof a c
tof a
tof c
tof e
T* 0
T b
T* f
T a
T* c
T* d
T e
tof a 0
tof c 0
tof d 0
tof d f
tof e a
tof e c
tof d e
tof a e
tof a c
tof b f
tof 0 f
tof 0 d
tof 0 a
tof 0 e
tof 0 b
T 0
T* b
tof d 0
tof d b
P d
Z e
tof 0 b
H d
tof d a
H b
tof b e
tof b c
tof f c
tof f b
tof b f
tof c b
tof e b
tof e
tof b
T f
T* c
T e
T* b
P* a
T* a
tof e f
tof b c
tof b e
tof e b
tof f c
tof f b
H b
tof f
H e
tof b 0
tof a 0
tof d 0
tof d f
tof c e
tof a b
tof d a
tof c a
tof d c
tof a d
tof e c
tof e a
tof e d
tof d
tof a
tof b
T e
T* c
P 0
T 0
T d
T* a
T* b
tof c e
tof d e
tof a e
tof a c
tof d 0
tof a 0
tof b 0
tof d a
tof c d
tof c b
tof f a
P a
tof f a
tof a f
tof 0
H d
H e
tof d 0
tof b 0
tof e d
tof e c
tof c
tof a
T e
T* d
T* 0
T c
tof d e
tof c e
tof 0 d
tof c d
tof b d
tof b 0
tof d 0
tof e 0
tof e c
T* e
tof c e
tof 0 e
tof c a
tof c e
P c
P a
P e
tof a c
tof f c
tof f a
tof a e
tof c f
tof c e
tof c a
tof c
H 0
H c
tof 0 f
T 0
T c
T* f
tof c 0
tof f c
tof 0 f
T* 0
T* f
T c
tof c 0
tof c f
tof f c
tof 0 c
tof e
H c
H f
tof f a
tof f e
tof f 0
tof a 0
tof c 0
tof e
T f
T* a
T c
T* 0
T e
tof a f
tof e f
tof e a
tof e 0
tof c 0
tof f a
tof f e
tof f 0
T* f
tof e f
tof a f
tof d
H a
tof 0 d
tof 0 e
tof d b
tof d
tof b
H 0
tof e d
tof e c
tof 0 c
tof f b
tof 0 b
tof 0 f
tof b f
tof e f
tof c e
tof d f
T* e
T c
T* b
tof c e
tof b e
tof b c
tof 0 c
tof e c
tof e b
T e
tof d e
tof f e
tof 0 e
tof d b
tof 0 b
tof c e
tof f d
tof c d
tof 0 d
tof b d
P 0
P* b
P* d
P e
tof b 0
tof c d
tof f e
tof c e
tof e f
tof d c
tof d f
tof 0 b
tof c
H b
tof b 0
tof f b
H d
tof f e
tof d c
tof d 0
tof f d
tof f b
tof d f
tof c 0
tof c
tof 0
T f
T* d
T c
tof d f
tof c f
tof d c
tof f c
tof f d
T* f
tof e f
tof c f
tof e b
tof f b
tof f e
tof d b
tof d e
tof d 0
tof f
tof d
H c
H b
tof 0 f
tof c f
tof 0 d
tof c e
tof d e
tof f e
tof f 0
tof f
tof 0
T c
T b
T f
T 0
tof b c
tof 0 b
tof 0 f
tof c 0
T* c
T 0
T* b
tof b c
tof b 0
tof 0 b
tof c b
tof c e
tof f e
tof f d
H b
tof c
H c
tof c f
tof c e
tof f
T c
T* e
T f
tof e c
tof f c
tof e f
tof c f
tof c e
T* c
tof f c
tof c d
tof e d
tof e c
H f
tof 0 c
H e
tof c d
tof f c
T f
T e
T* c
tof e f
tof c e
tof f c
T* f
T* c
T e
tof e f
tof e c
tof c e
tof f e
tof f d
H e
tof e c
tof b 0
tof b d
tof 0 d
tof f d
tof e f
tof c e
tof f
tof b
tof a c
tof c a
tof a c
tof c f
tof f c
tof c f
tof f b
tof b f
tof f b
tof b e
tof e b
tof b e
tof e d
tof d e
tof e d
tof d 0
tof 0 d
tof d 0
END
