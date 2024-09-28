# Original circuit
#   qubits: 20
#   qubits used: 20
#   H: 702
#   cnot: 2500
#   X: 0
#   T: 2457
#   P: 0
#   Z: 0
#   tdepth (by partitions): 1053
#   depth  (by critical paths): 3333
#   tdepth (by critical paths): 996
# Optimized circuit
#   qubits: 20
#   qubits used: 20
#   H: 358
#   cnot: 3427
#   X: 0
#   T: 1019
#   P: 441
#   Z: 17
#   tdepth (by partitions): 386
#   depth  (by critical paths): 2956
#   tdepth (by critical paths): 380
#   Time: 0.218 s
.v a b c d e f g h i j k l m n o 0 1 2 3 4
.i a b c d e f g h i j k l m n o
.o a b c d e f g h i j k l m n o 0 1 2 3 4

BEGIN
H 4
tof 4 f
tof 4 j
tof 4 1
tof j 1
tof f 1
T 4
T* f
T* j
T 1
tof f 4
tof j f
tof j 1
tof 4 f
tof 4 1
tof f j
H j
H 3
tof 3 k
tof j 2
tof j 3
tof 3 j
tof k 2
T j
T* 3
T* k
T 2
tof k j
tof k 3
tof k 2
tof j 3
tof 3 2
tof j k
tof 1 4
tof 1 f
tof 0 m
tof 1 2
tof 0 2
tof m 2
tof m 0
tof l 1
tof f 2
tof f 1
tof 4 f
tof c 4
tof c m
P f
P 3
Z j
Z l
P 0
P* o
P* 1
P* 2
P* c
P* m
P* 4
tof j 3
tof l f
tof 1 f
tof 0 j
tof c l
tof m j
tof l j
tof 1 l
tof h 0
tof c 1
tof c 2
tof h m
tof n 4
tof 4 c
tof 4 m
tof n h
tof h n
tof m 4
tof c m
tof 1 2
tof 1 4
tof 1 m
tof 1 c
tof 0 4
tof 0 m
tof 0 h
tof l m
tof j m
tof j 4
tof l c
tof l 1
tof f m
tof f c
tof f 1
P* f
P* 3
P* l
Z j
P 1
P n
P* 0
P* h
Z 2
P c
P* m
P* 4
tof m f
tof 4 l
tof l f
tof l 3
tof 1 j
tof m l
tof m j
tof c 2
tof m 4
tof h 2
tof 0 h
tof h c
tof n c
tof 1 c
tof c 2
tof n h
tof n 0
tof l 2
tof l c
tof l 1
tof l m
tof l 4
tof f 3
tof f 2
tof f c
tof f 1
tof f m
tof f 4
tof f l
H k
H 2
tof 2 1
tof k c
tof k 2
tof 2 k
tof 1 c
T k
T* 2
T* 1
T c
tof 1 k
tof 1 2
tof 1 c
tof k 2
tof 2 c
tof k 1
tof 2 k
P 2
P* k
tof k 2
tof 2 k
H 1
H c
tof c 0
tof 1 j
tof 1 c
tof c 1
tof 0 j
T 1
T* c
T* 0
T j
tof 0 1
tof 0 c
tof 0 j
tof 1 c
tof c j
tof 1 0
tof c 1
P c
P* 1
tof 1 c
tof c 1
H 0
H j
tof j h
tof 0 j
tof j 0
T 0
T* j
T* h
tof h j
tof j 0
tof 0 j
tof 0 h
T 0
tof h 0
tof j 0
tof h 0
P h
P* 0
tof 0 h
tof h 0
H j
H b
tof j b
tof j o
T j
T* o
T* b
tof b o
tof o j
tof j o
tof o b
T j
tof o j
Z b
Z j
tof b j
H o
tof o h
tof 0 h
tof 0 o
tof o 0
T* o
T h
tof h o
tof 0 h
tof o h
tof 0 o
P 0
P* o
tof o 0
tof 0 o
H h
tof h c
tof 1 c
tof 1 h
tof h 1
T* h
T c
tof c h
tof 1 c
tof h c
tof 1 h
P 1
P* h
tof h 1
tof 1 h
H c
tof c 2
tof k 2
tof k c
tof c k
T* c
T 2
tof 2 c
tof k 2
tof c 2
tof k c
P k
P* c
tof c k
tof k c
H 2
tof 2 l
tof 2 4
tof 2 m
tof 3 m
tof 3 2
tof m 3
T* 2
T* 3
T m
T* l
T* 4
tof m 2
tof l 3
tof 3 2
tof l m
tof 4 l
tof l 4
tof 3 m
tof 2 3
tof 2 4
tof 2 l
P 4
H m
tof l 2
tof 4 2
tof m l
tof l 4
tof 2 4
T m
T* 4
T* l
T 2
tof l m
tof 4 m
tof 2 l
tof l 2
tof 4 l
tof m 4
tof 4 2
T m
tof 4 m
tof 2 m
H 2
H l
tof c k
tof l k
tof l c
T* c
T k
tof k c
tof l k
tof c k
tof l c
P l
P* c
tof c l
tof l c
H k
tof h 1
tof k 1
tof k h
T* h
T 1
tof 1 h
tof k 1
tof h 1
tof k h
P k
P* h
tof h k
tof k h
H 1
tof o 0
tof 1 0
tof 1 o
T* o
T 0
tof 0 o
tof 1 0
tof o 0
tof 1 o
P 1
P* o
tof o 1
tof 1 o
H 0
tof 0 b
tof 0 j
T 0
T* j
T* b
tof b 0
tof j 0
tof 0 j
tof j b
T 0
tof j 0
tof b 0
H b
tof b 1
tof o 1
tof o b
tof b o
T* b
T 1
tof 1 b
tof o 1
tof b 1
tof o b
P o
P* b
tof b o
tof o b
H 1
tof 1 k
tof h k
tof h 1
tof 1 h
T* 1
T k
tof k 1
tof h k
tof 1 k
tof h 1
P h
P* 1
tof 1 h
tof h 1
H k
tof k 4
tof k l
tof c l
tof c k
tof k c
T* k
T l
T* 4
tof l k
tof c l
tof c 4
tof k l
tof c k
P c
P* k
tof k c
tof c k
H l
tof c 4
tof k 4
T 4
tof c 4
tof k 4
tof l m
tof l 4
P l
P* m
P* 4
tof 4 l
tof m l
tof l m
tof m 4
P l
tof m l
tof 4 l
H k
tof 1 h
tof k h
tof k 1
T* 1
T h
tof h 1
tof k h
tof 1 h
tof k 1
P k
P* 1
tof 1 k
tof k 1
H h
tof b o
tof h o
tof h b
T* b
T o
tof o b
tof h o
tof b o
tof h b
P h
P* b
tof b h
tof h b
H o
tof o j
tof o 0
T o
T* 0
T* j
tof j o
tof 0 o
tof o 0
tof 0 j
T o
tof 0 o
tof j o
H j
tof j h
tof b h
tof b j
tof j b
T* j
T h
tof h j
tof b h
tof j h
tof b j
P b
P* j
tof j b
tof b j
H h
tof h m
tof h k
tof 1 k
tof 1 h
tof h 1
T* h
T k
T* m
tof k h
tof 1 k
tof 1 m
tof h k
tof 1 h
P 1
P* h
tof h 1
tof 1 h
H k
tof k f
tof 1 m
tof h m
tof k c
tof h 1
tof c 1
tof m 1
T* c
T 1
T m
T* f
tof 1 c
tof k 1
tof h m
tof k f
tof c 1
tof c m
tof k c
P k
P* c
tof c k
tof k c
H h
tof k f
tof c f
tof j b
tof h b
tof h j
T* j
T b
T f
tof b j
tof h b
tof k f
tof c f
tof j b
tof h j
P h
P* j
tof j h
tof h j
H b
tof b 0
tof b o
T b
T* o
T* 0
tof 0 b
tof o b
tof b o
tof o 0
T b
tof o b
tof 0 b
H 0
tof 0 h
tof j h
tof j 0
tof 0 j
T* 0
T h
tof h 0
tof j h
tof 0 h
tof j 0
P j
P* 0
tof 0 j
tof j 0
H h
tof h k
tof k m
T* k
T m
tof h k
tof h m
tof k m
tof h m
P h
P* m
tof m h
tof h m
H c
tof m h
tof c h
tof c m
T* m
T h
tof h m
tof c h
tof m h
tof c m
P c
P* m
tof m c
tof c m
H h
tof 0 j
tof h j
tof h 0
T* 0
T j
tof j 0
tof h j
tof 0 j
tof h 0
P h
P* 0
tof 0 h
tof h 0
H j
tof j o
tof j b
T j
T* b
T* o
tof o j
tof b j
tof j b
tof b o
T j
tof b j
tof o j
H o
tof o h
tof 0 h
tof 0 o
tof o 0
T* o
T h
tof h o
tof 0 h
tof o h
tof 0 o
P 0
P* o
tof o 0
tof 0 o
H h
tof h f
tof h c
tof m c
tof m h
tof h m
T* h
T c
T* f
tof c h
tof m c
tof m f
tof h c
tof m h
P m
P* h
tof h m
tof m h
H c
tof m f
tof h f
tof c l
tof c k
tof h m
tof k m
tof f m
T* k
T m
T f
T* l
tof m k
tof c m
tof h f
tof c l
tof k m
tof k f
tof c k
P c
P* k
tof k c
tof c k
H h
tof l c
tof k c
tof o 0
tof h 0
tof h o
T* o
T 0
T c
tof 0 o
tof h 0
tof l c
tof k c
tof o 0
tof h o
P h
P* o
tof o h
tof h o
H 0
tof 0 b
tof 0 j
T 0
T* j
T* b
tof b 0
tof j 0
tof 0 j
tof j b
T 0
tof j 0
tof b 0
H b
tof b h
tof o h
tof o b
tof b o
T* b
T h
tof h b
tof o h
tof b h
tof o b
P o
P* b
tof b o
tof o b
H h
tof h c
tof c f
T* c
T f
tof h c
tof h f
tof c f
tof h f
P h
P* f
tof f h
tof h f
H k
tof f h
tof k h
tof k f
T* f
T h
tof h f
tof k h
tof f h
tof k f
P k
P* f
tof f k
tof k f
H h
tof b o
tof h o
tof h b
T* b
T o
tof o b
tof h o
tof b o
tof h b
P h
P* b
tof b h
tof h b
H o
tof o j
tof o 0
T o
T* 0
T* j
tof j o
tof 0 o
tof o 0
tof 0 j
T o
tof 0 o
tof j o
H j
tof j h
tof b h
tof b j
tof j b
T* j
T h
tof h j
tof b h
tof j h
tof b j
P b
P* j
tof j b
tof b j
H h
tof h k
tof f k
tof f h
tof h f
T* h
T k
tof k h
tof f k
tof h k
tof f h
P f
P* h
tof h f
tof f h
H k
tof l c
tof k c
tof k 3
tof k 1
tof l k
tof c l
T* l
T* k
T c
T* 3
T* 1
tof c k
tof 3 l
tof 1 l
tof k l
tof 1 c
tof 3 c
tof 3 1
tof k c
tof l c
tof l 1
tof l 3
T l
tof 3 l
tof 1 l
P 1
H 1
tof h f
tof 1 f
tof 1 h
T* h
T f
tof f h
tof 1 f
tof h f
tof 1 h
P 1
P* h
tof h 1
tof 1 h
H f
tof j b
tof f b
tof f j
T* j
T b
tof b j
tof f b
tof j b
tof f j
P f
P* j
tof j f
tof f j
H b
H k
tof b k
tof b o
T b
T* o
T* k
tof k o
tof o b
tof b o
tof o k
T b
tof o b
P k
P* b
tof k b
H o
tof o f
tof j f
tof j o
tof o j
T* o
T f
tof f o
tof j f
tof o f
tof j o
P j
P* o
tof o j
tof j o
H f
tof f c
tof f 3
tof f 1
tof h f
tof 1 h
T* f
T* 1
T h
T* c
T* 3
tof h f
tof c 1
tof 1 f
tof c h
tof 3 c
tof c 3
tof 1 h
tof f 1
tof f 3
tof f c
P 3
H h
tof c f
tof 3 f
tof h m
tof h l
tof h c
tof c h
tof l h
T* l
T* c
T h
T f
T* m
tof h l
tof m c
tof 3 f
tof f l
tof f c
tof m h
tof f h
tof c h
tof l c
tof c m
tof l f
P m
H 3
tof o j
tof 3 j
tof 3 o
T* o
T j
tof j o
tof 3 j
tof o j
tof 3 o
P 3
P* o
tof o 3
tof 3 o
H j
tof j k
tof j b
T j
T* b
T* k
tof k j
tof b j
tof j b
tof b k
T j
tof b j
tof k j
H b
H k
tof m b
tof k 3
tof o 3
tof o k
tof k o
tof b m
tof b c
T* k
T 3
T* b
T c
tof 3 k
tof o 3
tof c b
tof m c
tof b c
tof k 3
tof m l
tof m b
tof m f
tof o k
P o
P* k
P* m
P* b
P* f
P* l
tof k o
tof m b
tof m f
tof m l
tof o k
H 3
tof f l
tof 3 l
tof l f
T* f
T l
tof 3 f
tof 3 l
tof f l
tof 3 l
P 3
P* l
tof l 3
tof 3 l
H c
tof l 3
tof c 3
tof c l
T* l
T 3
tof 3 l
tof c 3
tof l 3
tof c l
P c
P* l
tof l c
tof c l
H 3
tof k o
tof 3 o
tof 3 k
T* k
T o
tof o k
tof 3 o
tof k o
tof 3 k
P 3
P* k
tof k 3
tof 3 k
H o
tof o 0
tof o j
T o
T* j
T* 0
tof 0 o
tof j o
tof o j
tof j 0
T o
tof j o
tof 0 o
H 0
tof 0 3
tof k 3
tof k 0
tof 0 k
T* 0
T 3
tof 3 0
tof k 3
tof 0 3
tof k 0
P k
P* 0
tof 0 k
tof k 0
H 3
tof 3 m
tof 3 c
tof l c
tof l 3
tof 3 l
T* 3
T c
T* m
tof c 3
tof l c
tof l m
tof 3 c
tof l 3
P l
P* 3
tof 3 l
tof l 3
H c
tof l m
tof 3 m
tof c f
tof c b
tof 3 l
tof b l
tof m l
T* b
T l
T m
T* f
tof m l
tof 3 l
tof c b
tof c l
tof 3 m
tof c f
tof b l
tof l m
tof c m
P c
P* m
tof m c
tof c m
H 3
tof f c
tof m c
tof m f
tof 0 k
tof 3 k
tof 3 0
tof c f
T* 0
T k
T c
tof k 0
tof 3 k
tof f c
tof m c
tof 0 k
tof 3 0
P 3
P* 0
tof 0 3
tof 3 0
H k
tof k j
tof k o
T k
T* o
T* j
tof j k
tof o k
tof k o
tof o j
T k
tof o k
tof j k
H j
tof j 3
tof 0 3
tof 0 j
tof j 0
T* j
T 3
tof 3 j
tof 0 3
tof j 3
tof 0 j
P 0
P* j
tof j 0
tof 0 j
H 3
tof 3 f
tof f l
T* f
T l
tof l f
tof 3 l
tof f l
tof 3 f
P 3
P* f
tof f 3
tof 3 f
H m
tof f 3
tof m 3
tof m f
T* f
T 3
tof 3 f
tof m 3
tof f 3
tof m f
P m
P* f
tof f m
tof m f
H 3
tof j 0
tof 3 0
tof 3 j
T* j
T 0
tof 0 j
tof 3 0
tof j 0
tof 3 j
P 3
P* j
tof j 3
tof 3 j
H 0
tof 0 o
tof 0 k
T 0
T* k
T* o
tof o 0
tof k 0
tof 0 k
tof k o
T 0
tof k 0
tof o 0
H o
tof o 3
tof j 3
tof j o
tof o j
T* o
T 3
tof 3 o
tof j 3
tof o 3
tof j o
P j
P* o
tof o j
tof j o
H 3
tof 3 n
tof 3 m
tof 3 1
tof f 3
tof m f
T* 3
T* m
T f
T* n
T* 1
tof f 3
tof n m
tof m 3
tof n f
tof 1 n
tof n 1
tof m f
tof 3 m
tof 3 1
tof 3 n
P 1
H f
tof c l
tof f l
tof n 3
tof 1 3
tof f h
tof c f
tof l c
T* c
T* f
T l
T 3
T* h
tof l c
tof h f
tof h l
tof n 3
tof 1 3
tof f l
tof c f
tof f h
P h
H 1
tof o j
tof 1 j
tof 1 o
T* o
T j
tof j o
tof 1 j
tof o j
tof 1 o
P 1
P* o
tof o 1
tof 1 o
H j
tof j k
tof j 0
T j
T* 0
T* k
tof k j
tof 0 j
tof j 0
tof 0 k
T j
tof 0 j
tof k j
H k
tof o 1
tof k n
tof k 1
tof o k
tof k o
T* k
T 1
T* n
tof 1 k
tof o 1
tof o n
tof k 1
tof o k
P o
P* k
tof k o
tof o k
H 1
tof 1 3
tof 1 m
tof k o
tof n o
tof 3 n
T* 3
T n
T o
T* m
tof 1 3
tof 1 n
tof k o
tof 1 m
tof 3 n
tof n o
tof 1 n
P 1
P* n
tof n 1
tof 1 n
H k
H 3
tof k 3
tof 1 m
tof n m
tof k j
T k
T* j
T* 3
T m
tof 3 j
tof j k
tof 1 m
tof n m
tof k j
tof j 3
T k
tof j k
P 3
P* k
tof 3 k
H j
tof j 1
tof o 1
tof 1 o
T* o
T 1
tof j o
tof j 1
tof o 1
tof j 1
P j
P* 1
tof 1 j
tof j 1
H n
tof 1 j
tof n j
tof n 1
T* 1
T j
tof j 1
tof n j
tof 1 j
tof n 1
P n
P* 1
tof 1 n
tof n 1
H j
tof j 3
tof j k
T j
T* k
T* 3
tof 3 j
tof k j
tof j k
tof k 3
T j
tof k j
tof 3 j
H k
H 3
tof 3 n
tof 1 n
tof 1 3
tof 3 1
T* 3
T n
tof n 3
tof 1 n
tof 3 n
tof k o
tof 1 3
P 1
P* 3
P k
P* o
tof 3 1
tof k o
tof 1 3
H n
tof n k
tof o m
tof n m
tof m o
T* o
T m
T* k
tof k m
tof m o
tof n m
tof n k
tof o m
T o
tof k o
P n
P* o
tof o n
tof n o
H o
tof 3 1
tof o 1
tof o 3
T* 3
T 1
tof 1 3
tof o 1
tof 3 1
tof o 3
P o
P* 3
tof 3 o
tof o 3
H 1
tof 1 0
tof 1 j
T 1
T* j
T* 0
tof 0 1
tof j 1
tof 1 j
tof j 0
T 1
tof j 1
tof 0 1
H 0
tof 0 o
tof 3 o
tof 3 0
tof 0 3
T* 0
T o
tof o 0
tof 3 o
tof 0 o
tof 3 0
P 3
P* 0
tof 0 3
tof 3 0
H o
tof o k
tof o m
tof k n
T* k
T n
T* m
tof m n
tof n k
tof o n
tof k n
tof k m
T k
tof m k
tof o k
tof o m
P o
P* m
tof m o
tof o m
H m
tof 0 3
tof m 3
tof m 0
T* 0
T 3
tof 3 0
tof m 3
tof 0 3
tof m 0
P m
P* 0
tof 0 m
tof m 0
H 3
H d
tof 3 d
tof 3 1
T 3
T* 1
T* d
tof d 1
tof 1 3
tof 3 1
tof 1 d
T 3
tof 1 3
P* d
P 3
tof d 3
H 1
tof 1 k
tof 1 m
tof 0 m
tof 0 1
tof 1 0
T* 1
T m
T* k
tof m 1
tof 0 m
tof 0 k
tof 1 m
tof 0 1
P 0
P* 1
tof 1 0
tof 0 1
H m
tof m c
tof o k
tof m k
tof 1 0
tof m o
tof 0 o
tof o 0
tof k o
T* 0
T k
T o
T* c
tof k 0
tof m k
tof 1 o
tof m c
tof 0 k
tof 0 o
tof m 0
P m
P* 0
tof 0 m
tof m 0
H 1
tof 1 d
tof m c
tof 0 c
tof 1 3
T 1
T* 3
T* d
T c
tof d 1
tof 3 1
tof m c
tof 0 c
tof 1 3
tof 3 d
T 1
tof 3 1
tof d 1
H d
tof o m
tof d m
tof m o
T* o
T m
tof d o
tof d m
tof o m
tof d m
P d
P* m
tof m d
tof d m
H 0
tof m d
tof 0 d
tof 0 m
T* m
T d
tof d m
tof 0 d
tof m d
tof 0 m
P 0
P* m
tof m 0
tof 0 m
H d
tof d 3
tof d 1
T d
T* 1
T* 3
tof 3 d
tof 1 d
tof d 1
tof 1 3
T d
tof 1 d
tof 3 d
H 3
tof 3 o
tof 3 0
tof m 0
tof m 3
tof 3 m
T* 3
T 0
T* o
tof 0 3
tof m 0
tof m o
tof 3 0
tof m 3
P m
P* 3
tof 3 m
tof m 3
H 0
tof o c
tof 0 c
tof m o
tof 3 o
tof 3 m
tof o m
tof c m
T* m
T c
T o
tof 0 m
tof 0 c
tof 3 o
tof m c
tof c o
tof 0 c
P 0
P* c
tof c 0
tof 0 c
H 3
H k
tof 3 d
T 3
T k
T* d
tof k 3
tof d k
tof 3 d
T* 3
T* d
T k
tof k 3
tof k d
tof d k
tof 3 k
H k
H d
tof o 0
tof d 0
tof d k
tof 0 o
T* o
T 0
T* k
tof k 0
tof 0 o
tof d 0
tof d k
tof o 0
T o
tof k o
tof d o
tof k 0
tof k o
P d
Z k
Z o
P* 0
tof o d
tof k d
tof 0 k
tof 0 o
tof k o
tof d o
tof k 0
P* d
tof o d
H o
tof o j
tof o 3
T o
T* 3
T* j
tof j o
tof 3 o
tof o 3
tof 3 j
T o
tof 3 o
tof j o
H j
tof j 0
tof 0 d
T* 0
T d
tof d 0
tof j d
tof 0 d
tof j 0
P j
P* 0
tof 0 j
tof j 0
H 4
tof h f
tof 4 f
tof 4 h
T* h
T f
tof f h
tof 4 f
tof h f
tof 4 h
P 4
P* h
tof h 4
tof 4 h
H f
tof c k
tof f k
tof f c
T* c
T k
tof k c
tof f k
tof c k
tof f c
P f
P* c
tof c f
tof f c
H k
tof 0 j
tof k j
tof k 0
T* 0
T j
tof j 0
tof k j
tof 0 j
tof k 0
P k
P* 0
tof 0 k
tof k 0
H j
tof j 3
tof j o
T j
T* o
T* 3
tof 3 j
tof o j
tof j o
tof o 3
T j
tof o j
tof 3 j
H 3
tof 3 k
tof 0 k
tof 0 3
tof 3 0
T* 3
T k
tof k 3
tof 0 k
tof 3 k
tof 0 3
P 0
P* 3
tof 3 0
tof 0 3
H k
tof k f
tof c f
tof c k
tof k c
T* k
T f
tof f k
tof c f
tof k f
tof c k
P c
P* k
tof k c
tof c k
H f
tof f m
tof f 4
tof h 4
tof h f
tof f h
T* f
T 4
T* m
tof 4 f
tof h 4
tof h m
tof f 4
tof h f
P h
P* f
tof f h
tof h f
H 4
tof 4 m
tof 4 b
tof f h
tof 4 h
tof m h
T 4
T* m
T* b
T h
tof b 4
tof m 4
tof f h
tof h b
tof b h
tof 4 m
tof m h
tof 4 b
T 4
tof m 4
tof h 4
H h
H f
tof k c
tof f c
tof f k
T* k
T c
tof c k
tof f c
tof k c
tof f k
P f
P* k
tof k f
tof f k
H c
tof 3 0
tof c 0
tof c 3
T* 3
T 0
tof 0 3
tof c 0
tof 3 0
tof c 3
P c
P* 3
tof 3 c
tof c 3
H 0
tof 0 o
tof 0 j
T 0
T* j
T* o
tof o 0
tof j 0
tof 0 j
tof j o
T 0
tof j 0
tof o 0
H o
tof o c
tof 3 c
tof 3 o
tof o 3
T* o
T c
tof c o
tof 3 c
tof o c
tof 3 o
P 3
P* o
tof o 3
tof 3 o
H c
tof k f
tof c b
tof c f
tof k c
tof c k
T* c
T f
T* b
tof f c
tof k f
tof k b
tof c f
tof k c
P k
P* c
tof c k
tof k c
H f
tof f 4
tof f m
tof c k
tof b k
tof 4 b
T* 4
T b
T k
T* m
tof f 4
tof f b
tof c k
tof f m
tof 4 b
tof b k
tof f b
P f
P* b
tof b f
tof f b
H c
tof f m
tof b m
tof o 3
tof c 3
tof c o
T* o
T 3
T m
tof 3 o
tof c 3
tof f m
tof b m
tof o 3
tof c o
P c
P* o
tof o c
tof c o
H 3
tof 3 j
tof 3 0
T 3
T* 0
T* j
tof j 3
tof 0 3
tof 3 0
tof 0 j
T 3
tof 0 3
tof j 3
H j
tof j c
tof o c
tof o j
tof j o
T* j
T c
tof c j
tof o c
tof j c
tof o j
P o
P* j
tof j o
tof o j
H c
tof c f
tof k f
tof f k
T* k
T f
tof c k
tof c f
tof k f
tof c f
P c
P* f
tof f c
tof c f
H b
tof f c
tof b c
tof b f
T* f
T c
tof c f
tof b c
tof f c
tof b f
P b
P* f
tof f b
tof b f
H c
tof j o
tof c o
tof c j
T* j
T o
tof o j
tof c o
tof j o
tof c j
P c
P* j
tof j c
tof c j
H o
tof o 0
tof o 3
T o
T* 3
T* 0
tof 0 o
tof 3 o
tof o 3
tof 3 0
T o
tof 3 o
tof 0 o
H 0
tof 0 c
tof j c
tof j 0
tof 0 j
T* 0
T c
tof c 0
tof j c
tof 0 c
tof j 0
P j
P* 0
tof 0 j
tof j 0
H c
tof c b
tof f b
tof c d
tof f c
tof c f
T* c
T b
T* d
tof b c
tof f b
tof f d
tof c b
tof f c
P f
P* c
tof c f
tof f c
H b
tof b 4
tof k m
tof b m
tof d f
tof c f
tof b k
tof k b
tof m b
T* b
T* k
T m
T f
T* 4
tof m k
tof 4 b
tof k b
tof 4 m
tof d f
tof c f
tof k m
tof b m
tof b 4
P 4
H c
tof 0 j
tof c j
tof c 0
T* 0
T j
tof j 0
tof c j
tof 0 j
tof c 0
P c
P* 0
tof 0 c
tof c 0
H j
H k
tof j o
T j
T k
T* o
tof k j
tof o k
tof j o
T* j
T* o
T k
tof k j
tof k o
tof o k
tof j k
H k
H o
tof 0 c
tof o d
tof o c
tof 0 o
tof o 0
T* o
T c
T* d
tof c o
tof 0 c
tof 0 d
tof o c
tof k f
tof 0 o
tof k d
P 0
P* o
Z k
P* f
P* d
tof o 0
tof k f
tof k d
tof 0 o
H c
tof d f
tof c f
tof o 0
tof 0 d
tof c 0
tof d 0
T* 0
T f
T d
tof o d
tof c d
tof f 0
tof d f
tof c d
tof 0 f
tof f d
tof c 0
P c
P* 0
tof 0 c
tof c 0
H o
H m
tof o j
T o
T m
T* j
tof m o
tof j m
tof o j
T* o
T* j
T m
tof m o
tof m j
tof j m
tof o m
H m
H j
tof m b
tof 4 b
tof j c
tof 4 m
tof j d
tof d f
tof m 4
T* d
T f
T* c
T* m
T b
tof c f
tof f d
tof j f
tof b m
tof 4 b
tof m b
tof d f
tof d c
T d
tof c d
tof j d
tof 4 m
tof 4 d
tof j c
P j
P* c
Z 4
P* m
P* d
tof c j
tof m 4
tof 4 m
tof m d
tof j c
H c
H l
tof c o
T c
P l
T l
T* o
tof l c
tof o l
tof c o
T* c
P* o
T* o
T l
tof l c
tof l o
tof o l
tof c l
H o
tof j d
tof o d
tof o f
tof o j
T* j
T d
T* f
tof f d
tof d j
tof o d
tof o f
tof j d
T j
tof f j
P o
P* j
tof j o
tof o j
H j
H d
tof j c
T j
T d
T* c
tof d j
tof c d
tof j c
T* j
T* c
T d
tof d j
tof d c
tof c d
tof j d
H d
H c
tof f o
tof c o
tof o f
T* f
T o
tof c f
tof c o
tof f o
tof c o
P c
P* o
tof o c
tof c o
H b
tof 0 f
tof b f
tof b 0
T* 0
T f
tof f 0
tof b f
tof 0 f
tof b 0
P b
P* 0
tof 0 b
tof b 0
H f
tof o c
tof f c
tof f o
T* o
T c
tof c o
tof f c
tof o c
tof f o
P f
P* o
tof o f
tof f o
H c
tof c l
tof c j
T c
T* j
T* l
tof l c
tof j c
tof c j
tof j l
T c
tof j c
tof l c
H l
tof l f
tof o f
tof o l
tof l o
T* l
T f
tof f l
tof o f
tof l f
tof o l
P o
P* l
tof l o
tof o l
H f
tof 0 b
tof f m
tof f b
tof 0 f
tof f 0
T* f
T b
T* m
tof b f
tof 0 b
tof 0 m
tof f b
tof 0 f
P 0
P* f
tof f 0
tof 0 f
H b
tof b 4
tof b k
tof f 0
tof b m
tof b 0
tof m 0
tof m b
tof 4 b
T* 4
T* m
T b
T 0
T* k
tof 0 4
tof f 4
tof b 4
tof b m
tof k b
tof f 0
tof 0 b
tof m b
tof 4 b
tof b k
tof 4 0
P k
H f
tof l o
tof f o
tof f l
T* l
T o
tof o l
tof f o
tof l o
tof f l
P f
P* l
tof l f
tof f l
H o
tof o j
tof o c
T o
T* c
T* j
tof j o
tof c o
tof o c
tof c j
T o
tof c o
tof j o
H c
H j
tof k c
tof j 4
tof j f
tof l f
tof l j
tof j l
tof c k
tof c b
T* j
T f
T* 4
T* c
T b
tof f j
tof l f
tof b c
tof l 4
tof k b
tof c b
tof j f
tof k m
tof l j
tof k c
P l
P* j
Z k
P* c
P* m
tof j l
tof k c
tof k m
tof l j
H f
tof f 0
tof j l
tof 4 l
tof 0 4
T* 0
T 4
T l
tof 4 0
tof f 4
tof j l
tof 0 4
tof 0 l
tof f 0
P f
P* 0
tof 0 f
tof f 0
H j
tof j 1
tof j o
T j
T* o
T* 1
tof 1 j
tof o j
tof j o
tof o 1
T j
tof o j
tof 1 j
H 1
tof l f
tof 1 f
tof f l
T* l
T f
tof 1 l
tof 1 f
tof l f
tof 1 f
P 1
P* f
tof f 1
tof 1 f
H b
tof 0 l
tof b l
tof b 0
T* 0
T l
tof l 0
tof b l
tof 0 l
tof b 0
P b
P* 0
tof 0 b
tof b 0
H l
tof f 1
tof l 1
tof l f
T* f
T 1
tof 1 f
tof l 1
tof f 1
tof l f
P l
P* f
tof f l
tof l f
H 1
tof 1 3
tof 1 j
T 1
T* j
T* 3
tof 3 1
tof j 1
tof 1 j
tof j 3
T 1
tof j 1
tof 3 1
H 3
tof 3 l
tof f l
tof f 3
tof 3 f
T* 3
T l
tof l 3
tof f l
tof 3 l
tof f 3
P f
P* 3
tof 3 f
tof f 3
H l
tof l c
tof l b
tof 0 b
tof 0 l
tof l 0
T* l
T b
T* c
tof b l
tof 0 b
tof 0 c
tof l b
tof 0 l
P 0
P* l
tof l 0
tof 0 l
H b
tof c k
tof b k
tof b m
tof 0 c
tof l c
tof l 0
tof b 0
tof c 0
T* 0
T k
T c
T* m
tof l c
tof b c
tof c 0
tof c k
tof b c
tof b m
tof 0 k
tof 0 c
tof b k
P b
P* k
tof k b
tof b k
H l
tof m b
tof k b
tof k m
tof 3 f
tof l f
tof l 3
tof b m
T* 3
T f
T b
tof f 3
tof l f
tof m b
tof k b
tof 3 f
tof l 3
P l
P* 3
tof 3 l
tof l 3
H f
tof f j
tof f 1
T f
T* 1
T* j
tof j f
tof 1 f
tof f 1
tof 1 j
T f
tof 1 f
tof j f
H j
tof j c
tof j l
tof 3 l
tof 3 j
tof j 3
T* j
T l
T* c
tof l j
tof 3 l
tof 3 c
tof j l
tof 3 j
P 3
P* j
tof j 3
tof 3 j
H l
tof l b
tof 0 c
tof l c
tof j 3
tof l 3
tof 0 3
tof c 3
tof c 0
T* 0
T c
T 3
T* b
tof j 3
tof l 3
tof c 0
tof 3 c
tof l 3
tof l b
tof 0 c
tof l 0
P l
P* 0
tof 0 l
tof l 0
H j
tof j 1
tof l b
tof 0 b
tof j f
T j
T* f
T* 1
T b
tof 1 j
tof f j
tof l b
tof 0 b
tof j f
tof f 1
T j
tof f j
tof 1 j
H 1
tof 1 3
tof 3 c
T* c
T 3
tof 3 c
tof 1 3
tof c 3
tof 1 c
P 1
P* c
tof c 1
tof 1 c
H 0
tof c 1
tof 0 1
tof 0 c
T* c
T 1
tof 1 c
tof 0 1
tof c 1
tof 0 c
P 0
P* c
tof c 0
tof 0 c
H 1
H 4
tof 1 4
tof 1 j
T 1
T* j
T* 4
tof 4 j
tof j 1
tof 1 j
tof j 4
T 1
tof j 1
P 4
P* 1
tof 4 1
H j
tof j 0
tof c 0
tof c j
tof j c
T* j
T 0
tof 0 j
tof c 0
tof j 0
tof c j
P c
P* j
tof j c
tof c j
H 0
tof l b
tof 0 b
tof b l
T* l
T b
tof 0 l
tof 0 b
tof l b
tof 0 b
P 0
P* b
tof b 0
tof 0 b
H k
tof b 0
tof k 0
tof k b
T* b
T 0
tof 0 b
tof k 0
tof b 0
tof k b
P k
P* b
tof b k
tof k b
H 0
tof j c
tof 0 c
tof 0 j
T* j
T c
tof c j
tof 0 c
tof j c
tof 0 j
P 0
P* j
tof j 0
tof 0 j
H c
tof c f
tof c 1
T c
T* 1
T* f
tof f c
tof 1 c
tof c 1
tof 1 f
T c
tof 1 c
tof f c
H f
tof f 0
tof f d
tof j f
tof 0 j
T* f
T* 0
T j
T* d
tof j f
tof d 0
tof d j
tof 0 j
tof f 0
tof 0 d
P d
H j
tof j k
tof b k
tof j f
tof b j
tof j b
T* j
T k
T* f
tof k j
tof b k
tof b f
tof j k
tof b j
P b
P* j
tof j b
tof b j
H k
tof l m
tof k m
tof f b
tof j b
tof l k
tof m l
T* l
T* k
T m
T b
tof m l
tof m k
tof f b
tof j b
tof k m
tof l m
P m
H j
tof d 0
tof j 0
tof j d
T* d
T 0
tof 0 d
tof j 0
tof d 0
tof j d
P j
P* d
tof d j
tof j d
H 0
H l
tof 0 c
T 0
T l
T* c
tof l 0
tof c l
tof 0 c
T* 0
T* c
T l
tof l 0
tof l c
tof c l
tof 0 l
H l
H c
tof m l
tof c j
tof d j
tof d c
tof c d
tof l m
T* c
T j
T* l
tof j c
tof d j
tof m l
tof c j
tof m b
tof m 0
tof m f
tof m 3
tof d c
P d
P* c
P* b
P* 3
Z f
P* 0
tof c d
tof m b
tof m 3
tof m f
tof m 0
tof d c
H j
tof j m
tof f b
tof j b
tof j f
T* f
T b
T* m
tof m b
tof b f
tof j b
tof j m
tof f b
T f
tof m f
P j
P* f
tof f j
tof j f
H f
tof c d
tof f d
tof f c
T* c
T d
tof d c
tof f d
tof c d
tof f c
P f
P* c
tof c f
tof f c
H d
H k
tof d 0
T d
T k
T* 0
tof k d
tof 0 k
tof d 0
T* d
T* 0
T k
tof k d
tof k 0
tof 0 k
tof d k
H k
H 0
tof l k
tof 0 b
tof 0 f
tof c 0
tof f c
tof k l
tof k m
T* 0
T* f
T c
T* k
T m
T* b
tof c 0
tof b f
tof f 0
tof b c
tof l k
tof l m
tof k m
tof f c
tof 0 f
tof 0 b
tof l 3
tof l f
tof l m
P b
P* l
P* m
P* f
P* 3
tof l f
tof l 3
tof l m
H c
tof c m
tof c 3
tof c j
tof j c
tof m c
T* m
T* j
T c
T* 3
tof c m
tof 3 j
tof 3 c
tof j c
tof m j
tof j 3
P 3
H k
tof 3 j
tof k j
tof k 3
T* 3
T j
tof j 3
tof k j
tof 3 j
tof k 3
P k
P* 3
tof 3 k
tof k 3
H j
tof b 0
tof j 0
tof j b
T* b
T 0
tof 0 b
tof j 0
tof b 0
tof j b
P j
P* b
tof b j
tof j b
H 0
tof 0 1
tof 0 d
T 0
T* d
T* 1
tof 1 0
tof d 0
tof 0 d
tof d 1
T 0
tof d 0
tof 1 0
H 1
tof 1 j
tof b j
tof b 1
tof 1 b
T* 1
T j
tof j 1
tof b j
tof 1 j
tof b 1
P b
P* 1
tof 1 b
tof b 1
H j
tof j c
tof j k
tof 3 k
tof 3 j
tof j 3
T* j
T k
T* c
tof k j
tof 3 k
tof 3 c
tof j k
tof 3 j
P 3
P* j
tof j 3
tof 3 j
H k
tof c l
tof k l
tof j 3
tof c 3
tof c k
tof l c
T* c
T* k
T l
T 3
tof 3 c
tof j c
tof l c
tof l k
tof j 3
tof 3 l
tof k l
tof c l
tof c 3
P l
H j
tof 1 b
tof j b
tof j 1
T* 1
T b
tof b 1
tof j b
tof 1 b
tof j 1
P j
P* 1
tof 1 j
tof j 1
H b
tof b d
tof b 0
T b
T* 0
T* d
tof d b
tof 0 b
tof b 0
tof 0 d
T b
tof 0 b
tof d b
H d
tof d k
tof 1 j
tof d j
tof d f
tof 1 d
tof j 1
T* d
T* 1
T j
T* k
T* f
tof j d
tof f 1
tof f j
tof f k
tof 1 j
tof d 1
tof 1 f
P f
H j
tof f k
tof j 3
tof j c
tof j m
tof f 1
tof j d
tof c j
tof 3 j
tof k 1
T* 3
T* c
T j
T k
T* m
T* d
tof j 3
tof m c
tof c 3
tof m j
tof 1 k
tof f k
tof d m
tof m d
tof c j
tof 3 c
tof 3 d
tof 3 m
P d
H f
tof f 0
tof m 3
tof d 3
tof f b
T f
T* b
T* 0
T 3
tof 0 f
tof b f
tof m 3
tof d 3
tof f b
tof b 0
T f
tof b f
tof 0 f
H 0
tof 0 1
tof 1 k
T* 1
T k
tof k 1
tof 0 k
tof 1 k
tof 0 1
P 0
P* 1
tof 1 0
tof 0 1
H d
tof 1 0
tof d 0
tof d 1
T* 1
T 0
tof 0 1
tof d 0
tof 1 0
tof d 1
P d
P* 1
tof 1 d
tof d 1
H 0
H a
tof 0 f
T 0
T a
T* f
tof a 0
tof f a
tof 0 f
T* 0
T* f
T a
tof a 0
tof a f
tof f a
tof 0 a
H a
H f
tof f j
tof f m
tof f d
tof 1 f
tof d 1
T* f
T* d
T 1
T* j
T* m
tof 1 f
tof m d
tof d f
tof m 1
tof m j
tof d 1
tof f d
tof f m
tof f j
P m
H 1
tof 1 j
tof f 3
tof 1 3
tof m f
tof 1 f
tof m 1
tof f 1
tof f m
tof 3 1
tof j f
T* 1
T* m
T 3
T f
T* j
tof 3 1
tof 3 m
tof j 3
tof 3 f
tof m 3
tof 1 3
tof 3 j
tof m f
P j
H f
H d
tof f 0
T f
T d
T* 0
tof d f
tof 0 d
tof f 0
T* f
T* 0
T d
tof d f
tof d 0
tof 0 d
tof f d
H d
H 0
tof 1 3
tof 0 3
tof 0 d
tof 0 f
tof d 0
tof d f
tof 3 1
T* 1
T 3
P* 0
T* 0
T* d
T* f
tof d 1
tof 3 1
tof 0 3
tof f 0
tof f d
tof 0 d
tof 0 f
tof 1 3
T 1
tof f 1
tof d 1
tof f 1
tof d f
tof 1 f
P d
P* f
P* 1
tof f d
tof d f
tof d 1
H f
H m
tof f 0
T f
T m
T* 0
tof m f
tof 0 m
tof f 0
T* f
T* 0
T m
tof m f
tof m 0
tof 0 m
tof f m
H m
H 0
tof m 3
tof j 3
tof j m
tof 0 1
tof 1 d
tof m j
T* 1
T d
T* m
T 3
tof d 1
tof 0 d
tof 3 m
tof j 3
tof m 3
tof 1 d
tof j 1
tof j m
tof j d
P 0
Z j
P* m
P* d
P* 1
tof 1 0
tof j 0
tof m j
tof d m
tof m d
tof d 1
tof j 1
tof 0 1
tof j d
tof j m
P* 0
H 3
tof 3 0
tof 3 1
T* 1
T 0
tof 0 1
tof 3 0
tof 1 0
tof 3 1
P 3
P* 1
tof 1 3
tof 3 1
H 0
tof 0 b
tof 0 f
T 0
T* f
T* b
tof b 0
tof f 0
tof 0 f
tof f b
T 0
tof f 0
tof b 0
H b
tof b d
tof b 3
tof b m
tof 1 b
tof 3 1
T* b
T* 3
T 1
T* d
T* m
tof 1 b
tof m 3
tof m 1
tof m d
tof 3 1
tof b 3
tof 3 m
tof 3 d
P m
H 1
tof 1 k
tof d j
tof 1 j
tof 1 c
tof 1 d
tof m 3
tof 1 3
tof d 1
tof d 3
tof j 1
T* 1
T* d
T j
T 3
T* k
T* c
tof j d
tof c 1
tof d 1
tof c j
tof m 3
tof 3 j
tof c k
tof d j
tof 1 j
tof 1 c
P c
H m
tof c k
tof c 1
tof m 4
tof m 0
tof k 1
T m
T* 0
T* 4
T k
tof 4 m
tof 0 m
tof 1 k
tof c k
tof m 0
tof 0 4
T m
tof 0 m
tof 4 m
H 0
H 4
tof 4 3
tof 4 b
tof 4 j
tof j 4
tof 3 4
T* 4
T* j
T 3
T* b
tof 3 4
tof b j
tof j 4
tof b 3
tof j 3
tof 4 j
tof 4 b
tof 0 3
P b
Z 0
P* 3
tof 0 3
H c
tof b 4
tof c 4
tof c b
T* b
T 4
tof 4 b
tof c 4
tof b 4
tof c b
P c
P* b
tof b c
tof c b
H 4
tof 4 f
tof 4 m
T 4
T* m
T* f
tof f 4
tof m 4
tof 4 m
tof m f
T 4
tof m 4
tof f 4
H f
tof f 3
tof f c
tof b c
tof b f
tof f b
T* f
T c
T* 3
tof c f
tof b c
tof b 3
tof f c
tof b f
P b
P* f
tof f b
tof b f
H c
tof k 1
tof c 1
tof c k
tof b 3
tof f 3
tof f b
tof 3 b
tof k c
tof 1 c
T* c
T* k
T 1
T 3
tof 1 c
tof 1 k
tof b 3
tof f 3
tof k 1
tof c 1
P 1
H f
H k
tof f 4
T f
T k
T* 4
tof k f
tof 4 k
tof f 4
T* f
T* 4
T k
tof k f
tof k 4
tof 4 k
tof f k
H k
H 4
tof 4 d
tof 3 b
tof 4 b
tof b 3
T* 3
T b
T* d
tof d b
tof b 3
tof 4 b
tof 4 d
tof 3 b
T 3
tof d 3
P 4
P* 3
P k
tof 3 4
tof 4 3
H 3
H c
tof 3 c
tof 3 f
T 3
T* f
T* c
tof c f
tof f 3
tof 3 f
tof f c
T 3
tof f 3
P c
P* 3
tof c 3
H f
tof f k
tof 4 d
tof f d
tof 4 f
tof d 4
T* 4
T* f
T d
T* k
tof d 4
tof d f
tof k d
tof f d
tof 4 d
tof d k
P k
H 4
tof 4 j
tof 4 3
T 4
T* j
T* 3
tof 3 j
tof j 4
tof 4 j
tof 4 3
T 4
tof 3 4
tof j 4
H j
H 3
tof j d
tof k d
tof k j
tof 1 k
tof k 1
tof 3 1
tof j 1
tof 3 4
tof j k
T 3
T* 1
T* 4
T* j
T d
T* k
tof 4 1
tof 1 3
tof d j
tof k 4
tof k d
tof j d
tof 4 d
tof 4 k
tof 3 1
tof 3 k
tof 3 d
tof 3 4
T 3
tof k 3
tof 1 3
tof k 3
tof k j
P* 3
P* j
tof k 3
H 1
H d
tof d m
tof d 3
T d
T* 3
T* m
P* 1
T* 1
tof m d
tof 3 d
tof d 3
tof 3 m
T d
tof 3 d
tof m d
H m
tof m j
tof m 1
tof j k
T* k
T j
T* 1
tof m k
tof m j
tof m 1
tof k j
tof m j
tof m 1
tof j 1
P m
P* j
P* 1
tof j m
tof m j
tof m 1
H k
tof d k
tof j 1
tof 1 m
tof d 1
tof m 1
tof k d
tof k 1
T d
T* 1
T* k
T m
tof k 1
tof 1 d
tof j m
tof k m
tof d 1
tof d m
tof d k
T d
tof k d
tof 1 d
H 1
H j
tof l 1
tof j c
tof j d
tof 1 l
T j
T* d
T* c
T* 1
tof c j
tof d j
tof l 1
tof j d
tof d c
T j
tof d j
tof c j
H d
H c
tof d l
tof 1 l
tof 1 d
tof c k
tof m k
tof k m
tof d 1
T* m
T k
T* d
T l
tof c m
tof c k
tof 1 d
tof 1 l
tof d l
tof m k
tof l b
tof l 0
P* l
P 1
P* b
P* 0
tof 1 l
tof b 1
tof 0 b
tof c k
tof b 0
tof 1 0
tof l 1
tof 1 b
P c
P* k
P* l
tof k c
tof 0 l
tof c k
H d
tof 4 m
tof d m
tof d 4
T* 4
T m
tof d 4
tof d m
tof 4 m
tof d m
P d
P* m
tof m d
tof d m
H 4
tof k c
tof 4 c
tof 4 k
T* k
T c
tof c k
tof 4 c
tof k c
tof 4 k
P 4
P* k
tof k 4
tof 4 k
H c
tof c 3
tof c j
T c
T* j
T* 3
tof 3 c
tof j c
tof c j
tof j 3
T c
tof j c
tof 3 c
H j
H 3
tof 3 l
tof 3 1
tof 3 4
tof k 4
tof k 3
tof 4 k
T* 3
T* k
T 4
T* l
T* 1
tof 4 3
tof 1 k
tof k 3
tof 1 4
tof 1 l
tof k 4
tof 3 k
tof 3 1
P 1
H 4
tof 1 l
tof 1 3
tof 4 b
tof 4 3
tof 4 d
tof m d
tof m 4
tof d m
tof l 3
T* 4
T* m
T d
T l
T* 3
T* b
tof d 4
tof l m
tof 3 m
tof 1 m
tof m 4
tof 3 l
tof 1 l
tof l d
tof b l
tof b 3
tof l b
tof l 3
tof m d
tof 4 m
tof 4 b
tof 4 3
tof 4 l
P b
H d
tof l 3
tof b 3
tof d 0
tof b l
tof d l
tof 3 l
T d
T* l
T* 0
T 3
tof 0 d
tof l d
tof b 3
tof 3 0
tof 0 3
tof d l
tof l 3
tof d 0
T d
tof l d
tof 3 d
H 3
H 1
H m
tof 1 c
T 1
P m
T m
T* c
tof m 1
tof c m
tof 1 c
T* 1
P* c
T* c
T m
tof m 1
tof m c
tof c m
tof 1 m
H c
tof c d
tof d 4
T* d
T 4
tof 4 d
tof c 4
tof d 4
tof c d
P c
P* d
tof d c
tof c d
H b
tof d c
tof b c
tof b d
T* d
T c
tof c d
tof b c
tof d c
tof b d
P b
P* d
tof d b
tof b d
H c
tof c m
tof c 1
T c
T* 1
T* m
tof m c
tof 1 c
tof c 1
tof 1 m
T c
tof 1 c
tof m c
H m
tof d b
tof m 0
tof m b
tof d m
tof m d
T* m
T b
T* 0
tof b m
tof d b
tof d 0
tof m b
tof d m
P d
P* m
tof m d
tof d m
H b
tof 0 4
tof b 4
tof m d
tof d 0
tof 0 d
tof d b
tof 4 d
T* d
T* b
T 4
T 0
tof 0 d
tof m d
tof 4 d
tof 4 b
tof m 0
tof 0 4
tof b 4
tof d 4
tof d 0
P 4
H m
tof m 1
tof m c
T m
T* c
T* 1
tof 1 m
tof c m
tof m c
tof c 1
T m
tof c m
tof 1 m
H c
H 1
tof 1 0
tof d 1
tof 0 d
T* 0
T* 1
T d
tof d 0
tof d 1
tof 1 d
tof 0 d
tof c m
P d
P c
P* m
tof c m
H 0
tof 0 c
tof 0 m
T 0
T* c
T* m
tof m c
tof c 0
tof 0 c
tof 0 m
T 0
tof m 0
tof c 0
H c
H b
tof m b
tof 0 c
tof 4 0
tof 0 4
tof c 4
tof c 0
tof b m
T* b
P 4
T 4
T* c
T* 0
tof c b
tof 4 b
tof 0 4
tof m c
tof 0 c
tof 4 c
tof b c
tof 4 0
T b
tof c b
P m
P* b
tof m b
H c
tof c 0
tof c b
T c
T* 0
T* b
tof b 0
tof 0 c
tof c 0
tof c b
T c
tof b c
tof 0 c
H 0
H 1
tof 0 c
tof 0 d
P* 0
T* 0
T 1
T* c
T* d
tof 1 0
tof c 1
tof c d
tof 1 d
tof 0 d
tof 0 c
T* 0
T* c
T 1
tof 1 0
tof 1 c
tof c 1
tof 0 1
H 1
tof 4 1
tof 1 b
tof 1 4
T* 1
T b
tof 4 1
tof 4 b
tof 1 b
tof 4 b
Z 4
P* b
tof 4 b
H 1
tof d c
tof 1 c
tof 1 d
T* d
T c
tof 1 d
tof 1 c
tof d c
tof 1 c
P 1
P* c
tof c 1
tof 1 c
H d
tof d m
tof d 0
T d
T* 0
T* m
tof m d
tof 0 d
tof d 0
tof 0 m
T d
tof 0 d
tof m d
H 0
H m
tof m 1
tof c 1
tof c m
tof m c
T* m
T 1
tof 1 m
tof c 1
tof m 1
tof c m
tof 0 4
P c
P* m
P 0
P* 4
tof m c
tof 0 4
tof c m
H 1
tof 1 0
tof 1 b
tof b 4
T* b
T 4
T* 0
tof 0 4
tof 4 b
tof 1 4
tof 1 0
tof b 4
T b
tof 0 b
P 1
P* b
tof b 1
tof 1 b
H b
tof m c
tof b c
tof b m
T* m
T c
tof c m
tof b c
tof m c
tof b m
P b
P* m
tof m b
tof b m
H c
H 4
tof c 4
tof c d
T c
T* d
T* 4
tof 4 c
tof d c
tof f d
tof d 4
tof c d
tof d f
T c
T* d
T 4
tof 4 d
tof f c
tof d c
tof f 4
tof d 4
tof f d
P f
P* d
tof f d
H c
tof c b
tof m c
tof b m
T* c
T* b
T m
tof m c
tof m b
tof b m
tof c m
P m
H b
tof 1 0
tof b 0
tof b 4
tof 1 b
tof 0 1
T* 1
T* b
T 0
T* 4
tof 0 1
tof 4 b
tof b 1
tof 4 0
tof b 0
tof 1 b
tof 1 4
P 4
H f
H c
tof f 1
tof 4 1
tof 4 f
tof f 4
T c
T* f
T 1
tof 1 f
tof f c
tof d f
tof d 1
tof 4 1
tof f 1
tof c f
tof c d
T* c
T* f
T d
tof d c
tof d f
tof f d
tof c d
H d
tof 4 f
P 4
P* f
tof f 4
tof 4 f
H 1
tof 1 m
tof m 1
tof d 1
tof d c
T d
T* m
T* 1
T* c
tof m d
tof c m
tof d m
tof d 1
tof m c
T d
tof c d
tof 1 d
tof 1 c
P 1
P* c
tof c 1
tof 1 c
H d
tof d o
tof d m
T d
T* m
T* o
tof o d
tof m d
tof d m
tof m o
T d
tof m d
tof o d
H o
tof d 1
tof o 1
tof c o
tof d c
tof 1 c
tof 1 d
T* o
T* d
T c
T* 1
tof c o
tof 1 d
tof 1 c
tof d 1
tof o 1
tof d c
P 1
H c
tof c 4
tof c f
tof d o
tof c o
tof c d
tof o d
tof f c
tof 4 c
T* 4
T* f
T c
T* o
tof c f
tof o 4
tof d 4
tof f 4
tof d o
tof o c
tof f c
tof 4 c
tof 4 o
P o
tof d 4
tof 4 d
H c
tof c d
tof c 4
T c
T* 4
T* d
tof d 4
tof 4 c
tof c 4
tof c d
T c
tof 4 c
tof d c
H 4
tof 1 d
tof 1 c
tof d c
tof 4 d
tof 4 1
tof 4 o
P 4
T 4
T* 1
T d
T* o
tof 1 4
tof d 1
tof o d
tof 1 d
tof 4 d
tof d o
tof o 1
Z 1
tof o 1
tof 1 c
H 4
tof c m
tof 4 m
tof 1 o
tof d o
tof 4 k
tof d 1
tof c 4
tof o 1
tof m c
T* 4
T* c
T m
T* k
T o
tof m 4
tof k c
tof o m
tof d m
tof 1 m
tof k m
tof d o
tof 1 o
tof c m
tof 4 c
tof c k
P k
tof o m
H 4
H m
tof m c
tof k m
tof m k
T k
T* m
T* c
tof c m
tof m k
tof k m
tof k c
T k
tof c k
tof m k
H m
H c
tof c 1
tof c o
tof c k
tof c m
tof k m
tof o c
tof 1 c
T* 1
T* o
T c
T* m
T* k
tof c 1
tof k o
tof k c
tof o c
tof 1 o
tof o k
tof o m
P k
P m
H f
tof f m
tof f o
T f
T* o
T* m
tof m f
tof o f
tof f o
tof o m
T f
tof o f
tof m f
H m
H d
tof k f
tof d f
tof d k
T* k
T f
tof f k
tof d f
tof k f
tof m 1
tof d k
P d
P* k
P m
P* 1
tof k d
tof m 1
tof d k
H f
H e
tof f o
T f
T e
T* o
tof e f
tof o e
tof f o
T* f
T* o
T e
tof e f
tof e o
tof o e
tof f e
H e
H o
tof o m
tof o 1
tof o d
tof k o
tof d k
T* o
T* d
T k
T* m
T* 1
tof k o
tof 1 d
tof d o
tof 1 k
tof 1 m
tof d k
tof o d
tof o 1
P 1
H k
tof 1 m
tof o c
tof k c
tof k o
T* o
T c
T m
tof k o
tof k c
tof 1 m
tof o c
tof o m
tof k c
P k
P* c
tof c k
tof k c
H 1
H f
tof 1 d
T 1
T f
T* d
tof f 1
tof d f
tof 1 d
T* 1
T* d
T f
tof f 1
tof f d
tof d f
tof 1 f
H f
H d
tof d m
tof f o
tof f 1
tof f d
tof 1 d
T* m
T f
T* 1
P o
T o
T* d
tof o f
tof f m
tof 1 f
tof d o
tof o 1
tof f d
T* 1
T m
P* f
T* f
T* d
T o
tof o f
tof f 1
tof f m
tof d f
tof o d
tof d o
tof f o
tof 1 m
tof 1 f
P o
tof f 1
H m
tof d m
tof f d
tof m f
tof m d
T* m
T* d
T f
tof f m
tof f d
tof d f
tof m f
P f
H m
tof m 1
tof m d
T m
T* 1
T* d
tof d m
tof 1 m
tof m 1
tof 1 d
T m
tof 1 m
tof d m
H d
H o
tof d f
tof d 1
tof o d
tof k 1
T o
T* d
T* f
P* 1
T* 1
tof f o
tof f d
tof k 1
tof d 1
tof o 1
tof o f
T* d
T o
tof f d
tof f o
tof d o
tof o 1
P o
P* 1
tof 1 o
tof o 1
H d
H f
tof 1 d
tof f k
tof f o
tof o m
tof o f
tof m f
tof d 1
T* f
T* o
T m
T 1
T* d
T* k
tof m f
tof m o
tof d 1
tof k m
tof 1 d
tof o m
tof f m
tof m k
P k
H o
tof o d
tof o m
tof 1 o
tof o f
tof d 1
T* d
T* o
T 1
T* m
T* f
tof 1 d
tof m d
tof o d
tof f o
tof f 1
tof f m
tof m f
tof o m
tof d m
tof d f
tof o 1
P f
tof m o
H m
tof m d
tof m o
tof o f
tof d f
T m
T* d
T* o
T f
tof o d
tof d m
tof f d
tof d f
tof m f
tof m o
T m
tof o m
tof f m
H f
H b
tof b o
tof b m
T b
T* o
T* m
P f
T f
tof m o
tof o b
tof b o
tof b m
T b
tof m b
tof o b
H o
H d
tof d m
tof k m
tof c d
tof c k
tof d c
tof d k
T* k
T m
T* d
tof m k
tof c m
tof c d
tof k m
tof c k
tof c f
P c
Z k
P* f
tof k c
tof c k
tof k f
H m
H l
tof m l
tof k c
tof d c
tof m b
T m
T* b
T* l
T c
tof l b
tof b m
tof d c
tof k l
tof k c
tof b l
tof m l
tof l k
T* b
T m
T* l
tof k b
tof k m
tof k l
tof b m
P k
H m
tof l c
tof m c
tof m l
tof m b
tof l m
tof c m
T* m
T* l
T c
T* b
tof c l
tof b m
tof l m
tof b c
tof l c
tof m c
tof m b
P b
H d
tof b m
tof d m
tof d b
T* b
T m
tof m b
tof d m
tof b m
tof d b
P d
P* b
tof b d
tof d b
H m
H g
tof m f
T m
T g
T* f
tof g m
tof f g
tof m f
T* m
T* f
T g
tof g m
tof g f
tof f g
tof m g
H g
H f
tof f d
tof b f
tof d b
T* f
T* d
T b
tof b f
tof b d
tof d b
tof f b
P b
H d
tof c l
tof d l
tof l c
T* c
T l
tof d c
tof d l
tof c l
tof d l
P d
P* l
tof l d
tof d l
H f
tof f d
tof f m
T f
T* d
T* m
tof m d
tof d f
tof f d
tof f m
T f
tof m f
tof d f
H d
H c
tof m c
tof f m
tof b d
tof b f
tof f b
tof m b
tof d f
tof c b
tof c m
T b
T* c
T* m
T* d
T f
tof m b
tof f d
tof c b
tof d c
tof f d
tof d f
tof m f
tof b c
tof c f
tof c d
tof b m
T b
tof m b
tof f b
H f
tof d f
tof b f
tof c f
tof c m
tof c b
tof l c
tof c l
tof d b
tof f b
tof f l
tof f c
tof f m
tof f d
P* l
T* l
T* f
T d
P* m
T* m
T* c
tof d f
tof b f
tof f l
tof m f
tof c m
tof b d
tof b c
tof m c
tof f c
tof l f
tof f m
tof l d
tof b l
P c
P* l
tof l c
tof c l
tof c b
H d
H n
tof f c
tof m c
tof b f
tof m f
tof d b
tof c f
T d
T n
T* b
T c
tof n d
tof b n
tof m c
tof f c
tof n c
tof d c
tof d b
T* d
T* b
T n
tof f d
tof n d
tof n b
tof f n
tof b n
tof d n
tof d f
H n
H b
tof l d
tof b d
tof b f
tof l b
tof l f
tof d f
tof d l
T* l
T* b
T d
tof d l
tof d b
tof b d
tof l d
P d
tof b f
H f
tof f b
tof f c
T f
T* b
T* c
tof c f
tof b f
tof f b
tof b c
T f
tof b f
tof c f
H c
tof c l
H m
tof b k
tof l k
tof l b
tof d l
tof l d
T d
T* l
T c
T b
T* k
tof m l
tof l d
tof k b
tof l b
tof d b
tof b k
T* d
T l
T b
tof b l
tof l d
tof m l
tof d l
tof d b
tof m b
P m
P* b
tof b m
tof m b
H l
H i
tof l c
T l
T i
T* c
tof i l
tof c i
tof l c
T* l
T* c
T i
tof i l
tof i c
tof c i
tof l i
H i
H c
tof c m
tof c b
T c
T* b
T* m
tof m b
tof b c
tof c b
tof c m
T c
tof m c
tof b c
H b
H m
tof m f
tof m k
T m
T* f
T* k
tof k m
tof f m
tof m f
tof f k
T m
tof f m
tof k m
H k
tof l m
tof f c
tof o m
H d
tof 0 j
tof 4 j
tof 1 a
tof g e
tof n a
tof i j
tof d e
tof l a
tof f e
tof j a
tof e j
tof e 4
tof 1 4
tof g 0
tof n 0
tof d i
tof o 4
tof l 0
tof c i
tof 4 0
tof i 4
tof n 1
tof i g
tof d g
tof o g
tof l 1
tof c g
tof 1 e
tof g 1
tof d g
tof o n
tof l n
tof c g
tof n g
tof i n
tof f i
tof c o
tof o i
tof o d
tof f l
tof c o
tof c l
tof c f
tof o c
tof d f
tof i f
tof d c
tof i l
tof n o
tof g o
tof n l
tof g d
tof 1 i
tof e i
tof 1 d
tof e n
tof 4 g
tof 0 g
tof 4 n
tof 0 1
tof j e
tof a e
tof j 1
tof a 4
tof b j
tof j b
tof b j
tof j d
tof d j
tof j d
tof d 4
tof 4 d
tof d 4
tof 4 2
tof 2 4
tof 4 2
tof 2 3
tof 3 2
tof 2 3
tof 3 h
tof h 3
tof 3 h
tof h n
tof n h
tof h n
tof n c
tof c n
tof n c
tof c 0
tof 0 c
tof c 0
tof k o
tof o k
tof k o
tof o m
tof m o
tof o m
tof m f
tof f m
tof m f
tof f 1
tof 1 f
tof f 1
END
