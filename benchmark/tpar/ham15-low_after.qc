# Original circuit
#   qubits: 17
#   qubits used: 17
#   H: 46
#   cnot: 259
#   X: 0
#   T: 161
#   P: 0
#   Z: 0
#   tdepth (by partitions): 69
#   depth  (by critical paths): 284
#   tdepth (by critical paths): 69
# Optimized circuit
#   qubits: 17
#   qubits used: 17
#   H: 38
#   cnot: 471
#   X: 0
#   T: 97
#   P: 19
#   Z: 1
#   tdepth (by partitions): 33
#   depth  (by critical paths): 360
#   tdepth (by critical paths): 33
#   Time: 0.004 s
.v a b c d e f g h i j k l m n o 0 1
.i a b c d e f g h i j k l m n o
.o a b c d e f g h i j k l m n o 0 1

BEGIN
tof k a
tof i e
tof k f
tof m e
tof i h
tof m h
tof k g
tof k h
tof m i
tof m l
tof n l
tof m n
tof l n
H i
tof 0 i
tof m l
tof 1 l
tof 1 m
tof m 0
tof m 1
tof l m
tof l 0
tof i 0
P m
T m
T 1
T i
T* l
T* 0
T o
T* n
tof 0 1
tof 1 m
tof 0 l
tof n 0
tof 0 n
tof 1 n
tof m n
tof m l
tof n i
tof m 0
T* m
T 1
tof 1 m
tof n 1
tof m 1
tof 1 0
tof m 0
H n
tof o a
tof 1 d
tof n e
tof n h
tof o g
tof o n
tof o 0
H 1
tof o n
tof 0 m
tof 0 1
tof o 0
tof i o
tof o i
tof 1 0
tof 0 i
tof 1 o
T i
T* 1
T* 0
T o
tof m 1
tof n 0
tof 1 i
tof 0 1
tof m o
tof n o
tof 0 o
tof i 1
tof i o
tof 1 0
tof 1 n
tof i m
H 0
tof 0 d
tof 0 m
H n
tof m i
tof n 0
tof l m
tof o n
tof n o
tof o l
tof m l
tof 0 m
tof 0 n
T n
T o
T* 0
T* l
T m
tof 0 n
tof i l
tof o n
tof l o
tof l 0
tof i m
tof l m
tof 0 m
tof o m
tof o l
tof o 0
tof o i
H l
tof 0 k
H 0
tof i o
tof i 0
tof m l
tof m n
tof l n
tof l m
tof i m
tof 0 i
tof i n
tof 0 l
T i
T* m
T* 0
T* n
T l
tof 0 i
tof o i
tof o m
tof o 0
tof o l
tof 0 n
tof 0 l
tof m n
tof i m
tof m l
tof i 0
P* m
tof i o
H 0
tof k b
tof 0 b
tof 0 k
tof k c
tof k d
tof k e
tof 0 e
tof 0 k
tof k 1
tof k j
tof 0 k
H l
tof 0 i
tof n i
tof l i
tof n 0
tof l o
tof l n
tof 0 o
tof 0 n
tof i o
tof i 0
tof 1 k
T l
T* 0
T* n
T i
tof 0 l
tof n 0
tof n i
tof l 0
tof l i
tof 0 n
tof 0 l
tof l o
Z 0
P* l
P* o
tof l 0
tof k 1
tof 1 k
tof 0 l
tof l o
H n
H i
tof o 0
tof i o
tof n i
tof i n
tof k 1
T n
T* i
T* o
tof o i
tof i n
tof n i
tof n o
T n
tof 0 n
tof o n
tof i n
tof o n
tof 0 n
P o
P* n
tof 0 n
tof n o
tof 1 k
tof k 1
tof o n
H i
H k
tof i m
T i
T k
T* m
tof k i
tof m k
tof i m
T* i
T* m
T k
tof k i
tof k m
tof m k
tof i k
tof o 0
H k
H m
tof 0 o
tof n 0
tof m l
tof m 0
tof n m
tof m n
T* m
T 0
T* l
tof o 0
tof 0 m
tof n 0
tof n l
tof m 0
tof n m
tof o m
P n
P* m
tof o m
tof m n
tof n m
tof n o
H 0
tof 0 n
tof m o
tof o l
tof 0 o
tof l o
T 0
T* n
T* o
T l
tof o n
tof m l
tof n 0
tof l n
tof l o
tof o l
tof n o
tof 0 o
tof o l
T 0
tof o 0
tof l 0
tof 0 n
H l
tof k b
tof k c
tof k d
tof k e
tof k g
tof k j
tof 0 k
H m
H 0
tof m i
T m
T 0
T* i
tof 0 m
tof i 0
tof m i
T* m
T* i
T 0
tof 0 m
tof 0 i
tof i 0
tof m 0
H 0
H i
tof o 0
tof n k
tof i 0
tof i o
tof o n
T* o
T n
T* 0
tof 0 n
tof n o
tof i n
tof i 0
tof o n
tof n 0
T o
tof k o
tof n o
tof 0 o
tof k o
tof 0 n
tof i n
tof o n
P i
P* o
P* n
tof f o
tof k o
tof o i
tof f n
tof h n
tof k n
tof 0 n
tof k f
tof 0 h
tof 0 k
tof n h
tof i o
tof i f
tof i h
tof i n
H o
H 0
tof o m
T o
T 0
T* m
tof 0 o
tof m 0
tof o m
T* o
T* m
T 0
tof 0 o
tof 0 m
tof m 0
tof o 0
H 0
H m
tof m f
tof f h
tof 0 m
tof 0 o
T* f
T h
P* 0
T* 0
T* m
T* o
tof m f
tof h f
tof 0 h
tof o 0
tof o m
tof 0 m
tof 0 o
tof f h
T* f
tof o f
tof m f
tof o f
P m
P f
tof j f
tof h c
tof h d
tof h e
tof o i
tof f j
H m
H h
tof m 0
T m
T h
T* 0
tof h m
tof 0 h
tof m 0
T* m
T* 0
T h
tof h m
tof h 0
tof 0 h
tof m h
H h
H 0
tof j o
tof 0 o
tof 0 j
T* j
T o
tof 0 j
tof 0 o
tof j o
tof h m
tof 0 o
P 0
P* o
P m
tof o 0
tof k a
tof j a
tof h m
tof h a
tof j c
tof h c
tof h d
tof h k
tof k 1
tof 0 o
H j
tof j 0
tof m j
tof j m
T m
T* 0
T* j
tof j 0
tof 0 m
tof m 0
tof m j
T m
tof j m
tof 0 m
H 0
tof 0 j
tof o j
tof o 0
tof 0 o
T* 0
T j
tof j 0
tof o j
tof 0 j
tof o 0
P o
P* 0
tof c o
tof c 0
tof n c
tof k c
tof m c
tof m e
tof h e
tof k n
tof m g
tof m n
tof m k
tof h m
tof m h
tof c n
tof o 0
tof o n
tof o c
H j
H h
tof j m
T j
T h
T* m
tof h j
tof m h
tof j m
T* j
T* m
T h
tof h j
tof h m
tof m h
tof j h
H h
H m
tof m 0
tof h j
tof m n
tof m j
tof m h
tof 0 h
tof 0 j
T m
T* 0
T* n
T h
T* j
tof n m
tof 0 m
tof h n
tof j h
tof h j
tof n j
tof m 0
tof 0 j
tof m n
T m
tof 0 m
tof j m
H j
tof 0 a
tof n a
tof 0 e
tof n e
tof n g
tof n k
tof n m
tof n 0
H n
tof n m
tof n h
T n
T* m
T* h
tof h m
tof m n
tof n m
tof n h
T n
tof h n
tof m n
H m
tof o b
tof d b
tof i a
tof g e
tof c a
tof 1 b
tof f e
tof b a
tof e b
tof e d
tof i d
tof g o
tof c o
tof f 1
tof k d
tof d o
tof 1 d
tof c i
tof 1 g
tof f g
tof k g
tof i e
tof g i
tof f g
tof k c
tof c g
tof 1 c
tof k 1
tof k f
tof n k
tof m k
tof m h
tof m n
tof n m
tof k m
tof k n
tof f 0
tof 1 0
tof f n
tof 1 h
tof c k
tof g k
tof c h
tof g f
tof i 1
tof e 1
tof i f
tof e c
tof d g
tof o g
tof d c
tof o i
tof b e
tof a e
tof b i
tof a d
tof j f
tof f j
tof j f
tof f i
tof i f
tof f i
tof i 1
tof 1 i
tof i 1
tof 1 l
tof l 1
tof 1 l
tof l h
tof h l
tof l h
tof h c
tof c h
tof h c
tof c o
tof o c
tof c o
tof o m
tof m o
tof o m
tof m 0
tof 0 m
tof m 0
END
