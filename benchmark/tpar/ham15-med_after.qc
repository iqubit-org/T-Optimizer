# Original circuit
#   qubits: 17
#   qubits used: 17
#   H: 164
#   cnot: 616
#   X: 0
#   T: 574
#   P: 0
#   Z: 0
#   tdepth (by partitions): 246
#   depth  (by critical paths): 817
#   tdepth (by critical paths): 240
# Optimized circuit
#   qubits: 17
#   qubits used: 17
#   H: 92
#   cnot: 759
#   X: 0
#   T: 230
#   P: 91
#   Z: 6
#   tdepth (by partitions): 86
#   depth  (by critical paths): 682
#   tdepth (by critical paths): 84
#   Time: 0.015 s
.v a b c d e f g h i j k l m n o 0 1
.i a b c d e f g h i j k l m n o
.o a b c d e f g h i j k l m n o 0 1

BEGIN
tof o j
tof n o
H 1
tof 1 l
tof 0 n
tof n m
tof 1 n
tof 1 0
tof n 0
tof m n
T 1
T* l
T* n
P 0
T 0
T o
T* m
tof n 1
tof l 1
tof 0 m
tof 1 l
tof l m
tof l n
T 1
tof l 1
tof n 1
H n
H m
tof m 0
tof m o
tof n m
tof m n
tof 0 o
T n
T* m
T* 0
tof 0 m
tof m n
tof n m
tof n 0
T n
tof 0 n
tof m n
tof 0 n
P 1
P l
P 0
P* n
tof n 0
tof 0 n
H m
tof l 1
P* 1
tof l 1
tof 0 o
H e
tof o 0
tof m o
T m
T e
T* o
tof e m
tof o e
tof m o
T* m
T* o
T e
tof 0 m
tof e m
tof e o
tof 0 e
tof o e
tof m e
tof m 0
H e
H o
tof o m
tof n m
tof n o
tof o 0
tof m n
tof m 0
T* o
T* n
T m
tof m o
tof m n
tof n m
tof o m
P m
tof o 0
H n
tof n 1
tof 0 o
tof 1 l
tof 1 0
T* 1
T l
T* 0
tof n 1
tof n l
tof o 0
tof n 0
tof 1 l
tof 1 0
tof n l
P n
P* l
tof l n
tof o 1
tof 0 1
tof 0 o
tof n l
H o
tof n 1
tof l 1
tof o n
tof o 0
T o
T* n
T* 0
T 1
tof 0 o
tof n o
tof l 1
tof 1 0
tof 0 1
tof o n
tof n 1
tof o 0
T o
tof n o
tof 1 o
H 1
H n
tof n 0
tof n 1
tof 0 1
T n
T* 1
T* 0
tof 1 n
tof 0 1
tof 1 0
tof n 0
T* 1
T n
tof n 1
tof 0 n
tof 1 n
P n
H 0
H l
tof m 0
tof l 0
tof l m
tof 0 l
T* m
T* l
T 0
tof 0 m
tof l m
tof l 0
tof m l
P 0
Z l
tof l 0
tof 0 l
P* 0
H m
H b
tof m n
T m
T b
T* n
tof b m
tof n b
tof m n
T* m
T* n
T b
tof b m
tof b n
tof n b
tof m b
H b
H n
tof n 1
tof n o
tof n 0
tof l n
tof 0 l
T* n
T* l
T 0
T* 1
T* o
tof 0 n
tof 1 l
tof l n
tof 1 0
tof o 1
tof 1 o
tof l 0
tof n l
tof n o
tof n 1
P o
H 0
tof 1 n
tof o n
T n
tof 1 n
tof o n
tof 0 n
tof 0 1
P 0
P* n
P* 1
tof 1 0
tof n 0
tof 0 n
tof n 1
P 0
tof n 0
tof 1 0
tof m k
H o
H m
tof o l
T o
T m
T* l
tof m o
tof l m
tof o l
T* o
T* l
T m
tof m o
tof m l
tof l m
tof o m
H m
H l
tof l 0
tof m n
tof m o
tof l m
tof 0 m
T* 0
T* m
P* n
T* n
T* o
tof n m
tof m 0
tof l m
tof o n
tof l o
tof n o
tof 0 o
tof 0 n
T* 0
T m
T* n
tof m 0
tof l m
tof l n
tof 0 m
tof m 0
tof o 0
P l
Z o
P 0
tof 0 o
tof m o
tof m 0
tof o 0
H j
tof j o
tof j 0
P j
T j
T* o
P* 0
T* 0
tof 0 o
tof o j
tof m o
tof m 0
tof j o
tof o m
T j
T* 0
T o
tof o 0
tof m j
tof 0 j
tof m o
H 1
tof l n
tof 1 n
tof 1 l
T* l
T n
tof n l
tof 1 n
tof l n
tof 1 l
P 1
P* l
tof l 1
tof 1 l
H n
H f
tof n 0
T n
T f
T* 0
tof f n
tof 0 f
tof n 0
T* n
T* 0
T f
tof f n
tof f 0
tof 0 f
tof n f
H f
H 0
tof 0 1
tof l 1
tof l 0
tof 0 l
T* 0
T 1
tof 1 0
tof l 1
tof 0 1
tof l 0
P l
P* 0
tof 0 l
tof l 0
H 1
tof 1 j
tof 1 o
P 1
P* j
P* o
tof o 1
tof j 1
tof 1 j
tof j o
P 1
tof j 1
tof o 1
H a
tof j a
tof l j
tof j l
tof a l
tof a j
tof a n
T* a
T n
T* j
tof j n
tof n a
tof l n
tof a n
tof a j
T a
tof j a
tof l a
tof l j
Z l
Z j
tof l j
H o
tof 0 a
tof o a
tof o 0
T* 0
T a
tof a 0
tof o a
tof 0 a
tof o 0
P o
P* 0
tof 0 o
tof o 0
H a
tof a l
tof a j
T a
T* j
T* l
tof l a
tof j a
tof a j
tof j l
T a
tof j a
tof l a
H l
tof l o
tof 0 o
tof 0 l
tof l 0
T* l
T o
tof o l
tof 0 o
tof l o
tof 0 l
P 0
P* l
tof l 0
tof 0 l
H o
tof n 1
tof o 1
tof a n
tof o n
tof n a
T* a
T 1
T* n
tof n 1
tof 1 a
tof o 1
tof o n
tof a 1
tof 1 n
T a
tof 1 a
tof n a
P o
P* a
tof a o
tof n 1
tof o a
H a
tof l 0
tof a 0
tof a l
T* l
T 0
tof 0 l
tof a 0
tof l 0
tof a l
P a
P* l
tof l a
tof a l
H 0
H d
tof 0 n
T 0
T d
T* n
tof d 0
tof n d
tof 0 n
T* 0
T* n
T d
tof d 0
tof d n
tof n d
tof 0 d
H d
H n
tof n a
tof l a
tof l n
tof n l
T* n
T a
tof a n
tof l a
tof n a
tof l n
P l
P* n
tof n l
tof l n
H a
tof a o
tof a 1
P a
P* o
P* 1
tof 1 a
tof o a
tof a o
tof o 1
P a
tof o a
tof 1 a
H c
tof 0 c
tof 0 o
tof c 0
tof c o
P 0
T 0
T* o
P* c
T* c
tof c 0
tof o 0
tof 0 o
tof o c
T 0
tof o 0
tof c 0
H 1
tof n l
tof 1 l
tof 1 n
T* n
T l
tof l n
tof 1 l
tof n l
tof 1 n
P 1
P* n
tof n 1
tof 1 n
H l
tof l c
tof l o
T l
T* o
T* c
tof c l
tof o l
tof l o
tof o c
T l
tof o l
tof c l
H c
tof c a
tof c 0
tof c 1
tof n c
tof 1 n
T* c
T* 1
T n
T* a
T* 0
tof n c
tof a 1
tof 1 c
tof a n
tof 0 a
tof a 0
tof 1 n
tof c 1
tof c 0
tof c a
P 0
H n
tof a c
tof n c
tof n 0
tof c 0
tof c a
T* a
T c
T 0
tof n a
tof n c
tof c 0
tof a c
tof n c
P n
P* c
tof c n
tof n c
H 0
H 1
tof 0 l
T 0
T 1
T* l
tof 1 0
tof l 1
tof 0 l
T* 0
T* l
T 1
tof 1 0
tof 1 l
tof l 1
tof 0 1
H 1
H l
tof n a
tof l a
tof l 1
tof a n
T* n
T a
T* 1
tof 1 a
tof a n
tof l a
tof l 1
tof n a
T n
tof 1 n
tof l n
tof 1 n
P l
Z 1
P* n
tof n l
tof 1 l
tof 1 n
tof l n
P* l
tof n l
H n
tof 0 o
tof n o
tof o 0
T* 0
T o
tof o 0
tof n o
tof 0 o
tof n 0
P n
P* 0
tof 0 n
tof n 0
H o
H a
tof 0 a
tof 0 n
tof a n
tof a 0
T 0
T* a
T* n
tof n a
tof a 0
tof 0 a
tof 0 n
T 0
tof n 0
tof a 0
H a
H n
tof a j
tof c a
tof n 1
tof 1 l
tof l 1
tof 1 n
tof a c
tof a 1
tof j 0
T* n
T l
T* j
T 0
T* a
T 1
tof 1 a
tof l n
tof a l
tof 0 j
tof c 0
tof c 1
tof a 1
tof j 0
tof n l
tof c a
tof c j
tof l n
P l
P* n
Z c
P* j
P* a
tof n l
tof c j
tof c a
tof l n
H 0
H 1
tof n l
tof 1 l
tof 1 n
T* n
T l
tof l n
tof 1 l
tof n l
tof 1 n
P 1
P* n
tof n 1
tof 1 n
H l
H g
tof l j
T l
T g
T* j
tof g l
tof j g
tof l j
T* l
T* j
T g
tof g l
tof g j
tof j g
tof l g
H g
H j
tof j 1
tof n j
tof 1 n
T* j
T* 1
T n
tof n j
tof n 1
tof 1 n
tof j n
P n
H 1
tof 1 c
tof c a
T* c
T a
tof a c
tof 1 a
tof c a
tof 1 c
P 1
P* c
tof c 1
tof 1 c
H j
tof j a
tof l j
tof j l
T l
T* a
T* j
tof j a
tof a l
tof l a
tof l j
T l
tof j l
tof a l
H a
H j
tof a j
tof n a
tof n l
tof a l
tof a n
tof j n
T n
T* l
T* j
T* a
tof l n
tof j l
tof a j
tof l a
tof n a
tof a j
T* l
T n
tof j l
tof j n
tof l n
H j
tof 1 j
tof c j
tof c n
tof n 1
tof 1 n
tof n c
tof j c
tof j n
tof j 1
T* n
T j
P* 1
T* 1
tof j n
tof c j
tof c 1
tof n j
tof c n
P c
P* n
tof c n
H j
tof a 1
tof j 1
tof j a
T* a
T 1
tof j a
tof j 1
tof a 1
tof j 1
P j
P* 1
tof 1 j
tof j 1
H a
H h
tof a l
T a
T h
T* l
tof h a
tof l h
tof a l
T* a
T* l
T h
tof h a
tof h l
tof l h
tof a h
H h
H l
tof l c
tof l j
tof 1 j
tof 1 l
tof l 1
T* l
T j
T* c
tof j l
tof 1 j
tof 1 c
tof l j
tof 1 l
P 1
P* l
tof l 1
tof 1 l
H j
tof 1 c
tof l c
tof l 1
tof j n
tof n 1
tof c 1
T* n
T 1
T c
tof c 1
tof l 1
tof j n
tof j 1
tof l c
tof n 1
tof 1 c
tof j c
P j
P* c
tof c j
tof j c
H l
H n
tof l a
T l
T n
T* a
tof n l
tof a n
tof l a
T* l
T* a
T n
tof n l
tof n a
tof a n
tof l n
H n
H a
tof n j
tof n c
P* n
T* n
P j
T j
T* c
tof a j
tof j n
tof a c
tof n j
tof n c
tof n 1
T* n
T 1
T c
tof 1 n
tof a 1
tof j c
tof n 1
tof 1 c
tof a n
P a
P* n
tof n a
tof a n
H c
tof n a
tof c a
tof c n
T* n
T a
tof a n
tof c a
tof n a
tof c n
P c
P* n
tof n c
tof c n
H a
tof a m
tof m l
T* m
T l
tof l m
tof a l
tof m l
tof a m
P a
P* m
tof m a
tof a m
H l
H k
tof m a
tof a k
tof a 1
tof a c
tof m a
tof k 1
tof k m
tof k c
tof k a
T* a
T* m
T k
T* 1
T* c
tof k a
tof c m
tof c k
tof c 1
tof m k
tof a m
tof m c
P c
H k
tof 1 c
tof k m
tof n k
tof m n
tof c 1
T* k
T* m
T n
T c
tof n k
tof n m
tof 1 c
tof m n
tof k n
tof k c
P n
H m
tof m c
tof m j
P m
P* c
P* j
tof j m
tof c m
tof m c
tof c j
P m
tof c m
tof j m
H a
tof a m
tof a c
T a
T* m
T* c
tof c a
tof m a
tof a m
tof m c
T a
tof m a
tof c a
H c
H 1
tof c k
H j
tof n k
tof j n
tof k j
T* n
T* k
T j
T c
tof j n
tof k j
tof n k
P j
P n
tof n j
tof j n
P* j
tof n j
H k
H i
tof k c
T k
T i
T* c
tof i k
tof c i
tof k c
T* k
T* c
T i
tof i k
tof i c
tof c i
tof k i
H i
H c
tof c j
tof n c
tof c n
T n
T* c
T* j
tof j c
tof c n
tof n c
tof n j
T n
tof j n
tof c n
H c
H j
tof j a
tof j m
T j
T* a
T* m
tof m j
tof a j
tof j a
tof a m
T j
tof a j
tof m j
H m
tof o b
tof d b
tof f 0
tof g e
tof h 0
tof i b
tof l e
tof n e
tof k 0
tof b 0
tof e b
tof e d
tof f d
tof g o
tof h o
tof l i
tof 1 d
tof k o
tof d o
tof i d
tof h f
tof i g
tof l g
tof 1 g
tof k f
tof f e
tof g f
tof l g
tof 1 h
tof k h
tof h g
tof i h
tof n i
tof 1 i
tof 1 l
tof k j
tof n a
tof n k
tof 1 j
tof 1 a
tof l n
tof i n
tof l a
tof i k
tof h 1
tof g 1
tof h k
tof g l
tof f i
tof e i
tof f l
tof e h
tof d g
tof o g
tof d h
tof o f
tof b e
tof 0 e
tof b f
tof 0 d
tof c o
tof o c
tof c o
tof o j
tof j o
tof o j
tof j l
tof l j
tof j l
tof l k
tof k l
tof l k
tof k 1
tof 1 k
tof k 1
tof 1 m
tof m 1
tof 1 m
tof m n
tof n m
tof m n
tof n a
tof a n
tof n a
tof a 0
tof 0 a
tof a 0
END
