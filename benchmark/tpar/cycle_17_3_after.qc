# Original circuit
#   qubits: 35
#   qubits used: 35
#   H: 1354
#   cnot: 4742
#   X: 0
#   T: 4739
#   P: 0
#   Z: 0
#   tdepth (by partitions): 2031
#   depth  (by critical paths): 6403
#   tdepth (by critical paths): 1941
# Optimized circuit
#   qubits: 35
#   qubits used: 35
#   H: 685
#   cnot: 6608
#   X: 2
#   T: 1944
#   P: 866
#   Z: 14
#   tdepth (by partitions): 585
#   depth  (by critical paths): 5231
#   tdepth (by critical paths): 562
#   Time: 1.207 s
.v a b c d e f g h i j k l m n o p q r s t 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14
.i a b c d e f g h i j k l m n o p q r s t
.o a b c d e f g h i j k l m n o p q r s t 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14

BEGIN
H 12
H 11
H 10
H 9
H 8
tof s a
tof t c
tof b a
tof c b
tof t c
tof d c
tof 9 f
tof 10 e
tof 8 h
tof 8 5
tof 9 10
tof h 8
tof e 10
P 5
Z 8
P* h
P* j
Z k
P l
P* n
Z o
P p
P* t
P* e
P* f
P* 10
P d
P* c
P* 12
P* b
P* a
tof h 5
tof 8 5
tof d t
tof c t
tof b t
tof 10 e
tof 9 e
tof 9 f
tof 9 10
tof b c
tof s a
tof c b
tof d b
tof d c
tof t a
tof t c
tof 8 h
H h
H 7
tof e f
tof 11 10
tof 7 i
tof 5 s
tof 7 3
tof h 3
tof 11 2
tof e 4
tof 9 6
tof s 5
tof i h
tof g 9
tof f 4
tof 10 2
tof b c
tof a r
tof t s
T 2
P 11
T 11
T* 10
T* 4
T e
T* f
P* g
T* g
T* 6
T* 9
T 7
T* 3
T* i
T h
P b
T b
T c
T* 5
P t
T t
T s
P a
T a
T* r
tof 10 11
tof 11 2
tof e 4
tof f 4
tof f e
tof 9 6
tof i 7
tof h 3
tof h i
tof c b
tof r a
tof a r
tof t s
tof 5 s
tof b c
tof i h
tof 7 3
tof 7 h
tof g 6
tof g 9
tof e f
tof 11 10
tof i 7
P i
P* 7
tof 7 i
tof i 7
H h
H 6
tof 6 j
tof h s
tof h 6
tof 6 h
tof j s
T h
T* 6
T* j
T s
tof j h
tof j 6
tof j s
tof h 6
tof 6 s
tof h j
tof 6 h
P 6
P* h
tof h 6
tof 6 h
H j
H s
tof s k
tof j 4
tof j s
tof s j
tof k 4
T j
T* s
T* k
T 4
tof k j
tof k s
tof k 4
tof j s
tof s 4
tof j k
tof s j
P s
P* j
tof j s
tof s j
H k
H 4
tof 4 l
tof k 3
tof k 4
tof 4 k
tof l 3
T k
T* 4
T* l
T 3
tof l k
tof l 4
tof l 3
tof k 4
tof 4 3
tof k l
tof 4 k
P 4
P* k
tof k 4
tof 4 k
H l
H 3
tof 3 m
tof l 2
tof l 3
tof 3 l
tof m 2
T l
T* 3
T* m
T 2
tof m l
tof m 3
tof m 2
tof l 3
tof 3 2
tof l m
tof 3 l
P 3
P* l
tof l 3
tof 3 l
H m
H 2
tof 2 n
tof m 1
tof 2 m
tof n 1
T 2
T* m
T* n
T 1
tof n 2
tof n m
tof n 1
tof 2 m
tof m 1
tof 2 n
tof m 2
P m
P* 2
tof 2 m
tof m 2
H n
H 1
tof 1 o
tof n 0
tof n 1
tof 1 n
tof o 0
T n
T* 1
T* o
T 0
tof o n
tof o 1
tof o 0
tof n 1
tof 1 0
tof n o
tof 1 n
P 1
P* n
tof n 1
tof 1 n
H o
H 0
tof 0 p
tof o 0
tof 13 o
tof o 13
tof 0 13
tof p o
T 13
T* 0
T* p
T o
tof p 13
tof p 0
tof p o
tof 13 0
tof 0 o
tof 13 p
tof 0 13
P 0
P* 13
tof 13 0
tof 0 13
H p
H q
tof o b
tof p 14
tof p o
tof b o
tof q b
tof q 14
T p
P q
T q
T* o
T* 14
T* b
tof o q
tof b 14
tof a o
tof a 14
tof 5 a
tof a 5
tof o 14
tof o 5
tof o a
tof q b
tof p b
tof b 14
tof q 5
tof q o
tof q a
T q
T* a
T* o
tof p q
tof 5 q
tof 5 a
tof 5 o
H p
tof 13 0
tof p 1
tof p 0
tof 13 p
tof 0 13
T* p
T* 13
T 0
T* 1
tof 0 p
tof 1 13
tof 13 p
tof 1 0
tof 13 0
tof p 13
tof p 1
P 1
H 0
tof 0 m
tof 0 p
tof n p
tof n 0
tof p n
T* 0
T* n
T p
T* m
tof p 0
tof m n
tof n 0
tof m p
tof n p
tof 0 n
tof 0 m
P m
H p
tof p 3
tof p 0
tof 2 p
tof 0 2
T* p
T* 0
T 2
T* 3
tof 2 p
tof 3 0
tof 0 p
tof 3 2
tof 0 2
tof p 0
tof p 3
P 3
H 2
tof 2 4
tof 2 p
tof l 2
tof p l
T* 2
T* p
T l
T* 4
tof l 2
tof 4 p
tof p 2
tof 4 l
tof p l
tof 2 p
tof 2 4
P 4
H l
tof k 2
tof l s
tof l 2
tof k l
tof 2 k
T* l
T* k
T 2
T* s
tof 2 l
tof s k
tof k l
tof s 2
tof k 2
tof l k
tof l s
P s
H 2
tof 2 6
tof 2 l
tof j l
tof j 2
tof l j
T* 2
T* j
T l
T* 6
tof l 2
tof 6 j
tof j 2
tof 6 l
tof j l
tof 2 j
tof 2 6
P 6
H l
tof l i
tof l 2
tof h l
tof 2 h
T* l
T* 2
T h
T* i
tof h l
tof i 2
tof 2 l
tof i h
tof 2 h
tof l 2
tof l i
P i
H h
tof h 8
tof h l
tof 7 h
tof l 7
T* h
T* l
T 7
T* 8
tof 7 h
tof 8 l
tof l h
tof 8 7
tof l 7
tof h l
tof h 8
P 8
H 7
tof 7 h
P 7
P* h
tof h 7
tof 7 h
H h
tof 9 e
tof 10 11
tof b e
tof h 11
tof h e
tof h g
tof f 9
tof b 9
tof 10 f
tof b 10
tof 10 b
tof f 10
tof h b
tof 9 10
tof g 10
tof e 9
tof 11 f
P* h
T* h
T* g
T 10
T 9
T* b
T 11
T f
T e
tof g h
tof b g
tof 11 9
tof f 9
tof b 10
tof b 9
tof 11 b
tof f 11
tof 11 f
tof b f
tof f e
tof 9 e
tof 10 e
tof g 10
tof h g
tof h e
tof g f
tof g 11
tof f b
tof f 9
tof f 11
tof f 10
P* 10
P* 9
P* 11
P* b
tof f 10
tof f 9
tof f 11
tof f b
H f
tof f h
P f
P* h
tof h f
tof f h
H h
tof 8 7
tof h 7
tof h 8
T* 8
T 7
tof 7 8
tof h 7
tof 8 7
tof h 8
P h
P* 8
tof 8 h
tof h 8
H 7
tof i l
tof 7 l
tof 7 i
T* i
T l
tof l i
tof 7 l
tof i l
tof 7 i
P 7
P* i
tof i 7
tof 7 i
H l
tof 6 2
tof l 2
tof l 6
T* 6
T 2
tof 2 6
tof l 2
tof 6 2
tof l 6
P l
P* 6
tof 6 l
tof l 6
H 2
tof s j
tof 2 j
tof 2 s
T* s
T j
tof j s
tof 2 j
tof s j
tof 2 s
P 2
P* s
tof s 2
tof 2 s
H j
tof 4 k
tof j k
tof j 4
T* 4
T k
tof k 4
tof j k
tof 4 k
tof j 4
P j
P* 4
tof 4 j
tof j 4
H k
tof 3 p
tof k p
tof k 3
T* 3
T p
tof p 3
tof k p
tof 3 p
tof k 3
P k
P* 3
tof 3 k
tof k 3
H p
tof m 0
tof p 0
tof p m
T* m
T 0
tof 0 m
tof p 0
tof m 0
tof p m
P p
P* m
tof m p
tof p m
H 0
tof 1 n
tof 0 n
tof 0 1
T* 1
T n
tof n 1
tof 0 n
tof 1 n
tof 0 1
P 0
P* 1
tof 1 0
tof 0 1
H n
H 13
tof q 13
tof e 13
tof n q
tof e q
tof 14 q
tof 14 e
tof n 14
tof q 14
tof q e
tof 13 q
tof 13 14
T n
T 14
T* e
T* q
T* 13
tof e 14
tof 13 q
tof e q
tof 14 13
tof n 13
tof 13 q
tof 14 e
T 14
tof n 14
tof e 14
H e
H n
tof n p
tof n 0
tof 1 0
tof 1 n
tof 0 1
T* n
T* 1
T 0
T* p
tof 0 n
tof p 1
tof 1 n
tof p 0
tof 1 0
tof n 1
tof n p
P p
P e
H 0
tof 0 k
tof 0 n
tof m 0
tof n m
T* 0
T* n
T m
T* k
tof m 0
tof k n
tof n 0
tof k m
tof n m
tof 0 n
tof 0 k
P k
H m
tof m j
tof m 0
tof 3 m
tof 0 3
T* m
T* 0
T 3
T* j
tof 3 m
tof j 0
tof 0 m
tof j 3
tof 0 3
tof m 0
tof m j
P j
H 3
tof 4 m
tof 3 2
tof 3 m
tof 4 3
tof m 4
T* 3
T* 4
T m
T* 2
tof m 3
tof 2 4
tof 4 3
tof 2 m
tof 4 m
tof 3 4
tof 3 2
P 2
H m
tof m l
tof m 3
tof s 3
tof s m
tof 3 s
T* m
T* s
T 3
T* l
tof 3 m
tof l s
tof s m
tof l 3
tof s 3
tof m s
tof m l
P l
H 3
tof 3 7
tof 3 m
tof 6 3
tof m 6
T* 3
T* m
T 6
T* 7
tof 6 3
tof 7 m
tof m 3
tof 7 6
tof m 6
tof 3 m
tof 3 7
P 7
H 6
tof 6 h
tof 6 3
tof i 6
tof 3 i
T* 6
T* 3
T i
T* h
tof i 6
tof h 3
tof 3 6
tof h i
tof 3 i
tof 6 3
tof 6 h
P h
H i
tof 8 6
tof i f
tof i 6
tof 8 i
tof 6 8
T* i
T* 8
T 6
T* f
tof 6 i
tof f 8
tof 8 i
tof f 6
tof 8 6
tof i 8
tof i f
P f
H 6
tof 6 i
P 6
P* i
tof i 6
tof 6 i
H i
tof 11 b
tof i b
tof i 9
tof 13 11
tof 11 13
tof i 13
tof 10 11
tof 9 10
tof b 11
P i
T i
T* 9
T 10
T* 13
T b
T 11
tof 9 i
tof 13 9
tof 13 10
tof b 13
tof b 11
tof 10 11
tof 13 b
tof 9 10
tof i 9
tof i 11
tof 9 b
tof b 13
tof b 10
P* 10
P* 13
tof b 10
tof b 13
H b
tof b i
P b
P* i
tof i b
tof b i
H i
tof f 6
tof i 6
tof i f
T* f
T 6
tof 6 f
tof i 6
tof f 6
tof i f
P i
P* f
tof f i
tof i f
H 6
tof h 8
tof 6 8
tof 6 h
T* h
T 8
tof 8 h
tof 6 8
tof h 8
tof 6 h
P 6
P* h
tof h 6
tof 6 h
H 8
tof 7 3
tof 8 3
tof 8 7
T* 7
T 3
tof 3 7
tof 8 3
tof 7 3
tof 8 7
P 8
P* 7
tof 7 8
tof 8 7
H 3
tof l m
tof 3 m
tof 3 l
T* l
T m
tof m l
tof 3 m
tof l m
tof 3 l
P 3
P* l
tof l 3
tof 3 l
H m
tof 2 s
tof m s
tof m 2
T* 2
T s
tof s 2
tof m s
tof 2 s
tof m 2
P m
P* 2
tof 2 m
tof m 2
H s
tof j 4
tof s 4
tof s j
T* j
T 4
tof 4 j
tof s 4
tof j 4
tof s j
P s
P* j
tof j s
tof s j
H 4
tof k 0
tof 4 0
tof 4 k
T* k
T 0
tof 0 k
tof 4 0
tof k 0
tof 4 k
P 4
P* k
tof k 4
tof 4 k
H 0
tof p n
tof 0 n
tof 0 p
T* p
T n
tof n p
tof 0 n
tof p n
tof 0 p
P 0
P* p
tof p 0
tof 0 p
H n
H 1
tof 14 1
tof 11 1
tof n 14
tof 11 14
tof q 14
tof q n
tof q 11
tof n q
tof 14 n
tof 14 11
tof 1 14
tof 1 n
T q
T n
T* 11
T* 14
T* 1
tof 11 n
tof 1 14
tof 11 14
tof n 1
tof q 1
tof 1 14
tof n 11
T n
tof q n
tof 11 n
H 11
H q
tof q 4
tof q 0
tof p q
tof 0 p
T* q
T* 0
T p
T* 4
tof p q
tof 4 0
tof 0 q
tof 4 p
tof 0 p
tof q 0
tof q 4
P 4
Z 11
H p
tof p s
tof p q
tof k p
tof q k
T* p
T* q
T k
T* s
tof k p
tof s q
tof q p
tof s k
tof q k
tof p q
tof p s
P s
H k
tof j p
tof k m
tof k p
tof j k
tof p j
T* k
T* j
T p
T* m
tof p k
tof m j
tof j k
tof m p
tof j p
tof k j
tof k m
P m
H p
tof p 3
tof p k
tof 2 k
tof 2 p
tof k 2
T* p
T* 2
T k
T* 3
tof k p
tof 3 2
tof 2 p
tof 3 k
tof 2 k
tof p 2
tof p 3
P 3
H k
tof k 8
tof k p
tof l k
tof p l
T* k
T* p
T l
T* 8
tof l k
tof 8 p
tof p k
tof 8 l
tof p l
tof k p
tof k 8
P 8
H l
tof l 6
tof l k
tof 7 l
tof k 7
T* l
T* k
T 7
T* 6
tof 7 l
tof 6 k
tof k l
tof 6 7
tof k 7
tof l k
tof l 6
P 6
H 7
tof h l
tof 7 i
tof 7 l
tof h 7
tof l h
T* 7
T* h
T l
T* i
tof l 7
tof i h
tof h 7
tof i l
tof h l
tof 7 h
tof 7 i
P i
H l
tof l b
tof l 7
tof f 7
tof f l
tof 7 f
T* l
T* f
T 7
T* b
tof 7 l
tof b f
tof f l
tof b 7
tof f 7
tof l f
tof l b
P b
H 7
tof 7 l
P 7
P* l
tof l 7
tof 7 l
H l
tof l 13
tof 10 l
tof 1 10
tof 10 1
tof l 1
tof 13 10
T 1
T* 13
T* l
T 10
tof 13 1
tof l 13
tof l 10
tof 1 13
tof 1 10
tof 13 l
H l
tof l 1
P l
P* 1
tof 1 l
tof l 1
H 1
tof b 7
tof 1 7
tof 1 b
T* b
T 7
tof 7 b
tof 1 7
tof b 7
tof 1 b
P 1
P* b
tof b 1
tof 1 b
H 7
tof i f
tof 7 f
tof 7 i
T* i
T f
tof f i
tof 7 f
tof i f
tof 7 i
P 7
P* i
tof i 7
tof 7 i
H f
tof 6 h
tof f h
tof f 6
T* 6
T h
tof h 6
tof f h
tof 6 h
tof f 6
P f
P* 6
tof 6 f
tof f 6
H h
tof 8 k
tof h k
tof h 8
T* 8
T k
tof k 8
tof h k
tof 8 k
tof h 8
P h
P* 8
tof 8 h
tof h 8
H k
tof 3 p
tof k p
tof k 3
T* 3
T p
tof p 3
tof k p
tof 3 p
tof k 3
P k
P* 3
tof 3 k
tof k 3
H p
tof m 2
tof p 2
tof p m
T* m
T 2
tof 2 m
tof p 2
tof m 2
tof p m
P p
P* m
tof m p
tof p m
H 2
tof s j
tof 2 j
tof 2 s
T* s
T j
tof j s
tof 2 j
tof s j
tof 2 s
P 2
P* s
tof s 2
tof 2 s
H j
tof 4 q
tof j q
tof j 4
T* 4
T q
tof q 4
tof j q
tof 4 q
tof j 4
P j
P* 4
tof 4 j
tof j 4
H q
H 0
tof q 0
tof 14 0
tof q n
tof 14 q
tof 14 10
tof q 14
tof n 10
tof 0 q
tof 0 10
T 14
T q
T* n
T* 0
T* 10
tof n q
tof 10 0
tof n 0
tof q 10
tof 14 10
tof 10 0
tof q n
T q
tof 14 q
tof n q
H n
H 14
tof 14 2
tof 14 j
tof 4 14
tof j 4
T* 14
T* j
T 4
T* 2
tof 4 14
tof 2 j
tof j 14
tof 2 4
tof j 4
tof 14 j
tof 14 2
P 2
P* n
H 4
tof s 14
tof 4 p
tof 4 14
tof s 4
tof 14 s
T* 4
T* s
T 14
T* p
tof 14 4
tof p s
tof s 4
tof p 14
tof s 14
tof 4 s
tof 4 p
P p
H 14
tof 14 k
tof 14 4
tof m 4
tof m 14
tof 4 m
T* 14
T* m
T 4
T* k
tof 4 14
tof k m
tof m 14
tof k 4
tof m 4
tof 14 m
tof 14 k
P k
H 4
tof 4 h
tof 4 14
tof 3 4
tof 14 3
T* 4
T* 14
T 3
T* h
tof 3 4
tof h 14
tof 14 4
tof h 3
tof 14 3
tof 4 14
tof 4 h
P h
H 3
tof 3 f
tof 3 4
tof 8 3
tof 4 8
T* 3
T* 4
T 8
T* f
tof 8 3
tof f 4
tof 4 3
tof f 8
tof 4 8
tof 3 4
tof 3 f
P f
H 8
tof 6 3
tof 8 7
tof 8 3
tof 6 8
tof 3 6
T* 8
T* 6
T 3
T* 7
tof 3 8
tof 7 6
tof 6 8
tof 7 3
tof 6 3
tof 8 6
tof 8 7
P 7
H 3
tof 3 1
tof 3 8
tof i 8
tof i 3
tof 8 i
T* 3
T* i
T 8
T* 1
tof 8 3
tof 1 i
tof i 3
tof 1 8
tof i 8
tof 3 i
tof 3 1
P 1
H 8
tof 8 l
tof 8 3
tof b 8
tof 3 b
T* 8
T* 3
T b
T* l
tof b 8
tof l 3
tof 3 8
tof l b
tof 3 b
tof 8 3
tof 8 l
P l
H b
tof b d
tof 10 c
tof d c
tof b 10
tof d 10
tof c 10
P b
T b
T* 10
T* d
T c
tof 10 b
tof d 10
tof d c
tof b 10
tof b c
tof 10 d
tof d 8
P* 8
tof d 8
H d
tof l 8
tof d 8
tof d l
T* l
T 8
tof 8 l
tof d 8
tof l 8
tof d l
P d
P* l
tof l d
tof d l
H 8
tof 1 3
tof 8 3
tof 8 1
T* 1
T 3
tof 3 1
tof 8 3
tof 1 3
tof 8 1
P 8
P* 1
tof 1 8
tof 8 1
H 3
tof 7 i
tof 3 i
tof 3 7
T* 7
T i
tof i 7
tof 3 i
tof 7 i
tof 3 7
P 3
P* 7
tof 7 3
tof 3 7
H i
tof f 6
tof i 6
tof i f
T* f
T 6
tof 6 f
tof i 6
tof f 6
tof i f
P i
P* f
tof f i
tof i f
H 6
tof h 4
tof 6 4
tof 6 h
T* h
T 4
tof 4 h
tof 6 4
tof h 4
tof 6 h
P 6
P* h
tof h 6
tof 6 h
H 4
tof k 14
tof 4 14
tof 4 k
T* k
T 14
tof 14 k
tof 4 14
tof k 14
tof 4 k
P 4
P* k
tof k 4
tof 4 k
H 14
tof p m
tof 14 m
tof 14 p
T* p
T m
tof m p
tof 14 m
tof p m
tof 14 p
P 14
P* p
tof p 14
tof 14 p
H m
tof 2 s
tof m s
tof m 2
T* 2
T s
tof s 2
tof m s
tof 2 s
tof m 2
P m
P* 2
tof 2 m
tof m 2
H s
H j
tof s j
tof 0 j
tof s q
tof c q
tof 0 q
tof 0 s
tof 0 c
tof s 0
tof q c
tof j s
tof j q
T 0
T s
T* c
T* j
T* q
tof c s
tof q j
tof c j
tof s q
tof 0 q
tof q j
tof s c
T s
tof 0 s
tof c s
H c
H 0
tof 2 m
tof 0 14
tof 0 m
tof 2 0
tof m 2
T* 0
T* 2
T m
T* 14
tof m 0
tof 14 2
tof 2 0
tof 14 m
tof 2 m
tof 0 2
tof 0 14
P 14
H m
tof m 4
tof m 0
tof p 0
tof p m
tof 0 p
T* m
T* p
T 0
T* 4
tof 0 m
tof 4 p
tof p m
tof 4 0
tof p 0
tof m p
tof m 4
P 4
H 0
tof 0 6
tof 0 m
tof k 0
tof m k
T* 0
T* m
T k
T* 6
tof k 0
tof 6 m
tof m 0
tof 6 k
tof m k
tof 0 m
tof 0 6
P 6
H k
tof k i
tof k 0
tof h k
tof 0 h
T* k
T* 0
T h
T* i
tof h k
tof i 0
tof 0 k
tof i h
tof 0 h
tof k 0
tof k i
P i
H h
tof f k
tof h 3
tof h k
tof f h
tof k f
T* h
T* f
T k
T* 3
tof k h
tof 3 f
tof f h
tof 3 k
tof f k
tof h f
tof h 3
P 3
H k
tof k 8
tof k h
tof 7 h
tof 7 k
tof h 7
T* k
T* 7
T h
T* 8
tof h k
tof 8 7
tof 7 k
tof 8 h
tof 7 h
tof k 7
tof k 8
P 8
H h
tof h d
tof h k
tof 1 h
tof k 1
T* h
T* k
T 1
T* d
tof 1 h
tof d k
tof k h
tof d 1
tof k 1
tof h k
tof h d
P d
H 1
tof 1 b
tof 1 10
tof 1 h
tof l h
tof l 1
tof q l
tof l q
tof 1 l
tof h q
tof 10 l
tof b l
T* 1
T* q
T h
T* b
T* 10
T l
tof h 1
tof b q
tof q 1
tof b h
tof 10 b
tof 10 l
tof b 10
tof q h
tof 1 q
tof 1 l
tof 1 10
tof 1 b
P 10
H h
tof h b
tof l 1
tof b 1
tof l h
tof h l
tof b h
tof 1 h
T l
T* h
T* b
T 1
tof h l
tof b h
tof b 1
tof l h
tof l 1
tof h b
H b
H 10
tof d q
tof 10 q
tof 10 d
tof 12 b
tof 12 13
tof 12 9
tof 12 g
tof b 13
tof b g
tof b 12
tof b 9
T* d
T q
T* b
T* g
T* 9
T* 13
tof q d
tof 10 q
tof 12 b
tof 12 g
tof 12 9
tof 12 13
tof d q
tof 10 d
P 10
P* d
Z 12
tof d 10
tof 10 d
H q
tof 8 k
tof q k
tof q 8
T* 8
T k
tof k 8
tof q k
tof 8 k
tof q 8
P q
P* 8
tof 8 q
tof q 8
H k
tof 3 7
tof k 7
tof k 3
T* 3
T 7
tof 7 3
tof k 7
tof 3 7
tof k 3
P k
P* 3
tof 3 k
tof k 3
H 7
tof i f
tof 7 f
tof 7 i
T* i
T f
tof f i
tof 7 f
tof i f
tof 7 i
P 7
P* i
tof i 7
tof 7 i
H f
tof 6 0
tof f 0
tof f 6
T* 6
T 0
tof 0 6
tof f 0
tof 6 0
tof f 6
P f
P* 6
tof 6 f
tof f 6
H 0
tof 4 m
tof 0 m
tof 0 4
T* 4
T m
tof m 4
tof 0 m
tof 4 m
tof 0 4
P 0
P* 4
tof 4 0
tof 0 4
H m
tof 14 p
tof m p
tof m 14
T* 14
T p
tof p 14
tof m p
tof 14 p
tof m 14
P m
P* 14
tof 14 m
tof m 14
H p
H 2
tof j s
tof p j
tof p 1
tof s j
tof 2 s
tof 2 1
T p
T 2
T* j
T* 1
T* s
tof j 2
tof j 1
tof 1 s
tof 2 1
tof p 1
tof 1 s
tof 2 j
T 2
tof p 2
tof j 2
H j
H p
tof p 0
tof p m
tof 14 m
tof 14 p
tof m 14
T* p
T* 14
T m
T* 0
tof m p
tof 0 14
tof 14 p
tof 0 m
tof 14 m
tof p 14
tof p 0
P 0
P j
H m
tof m f
tof m p
tof 4 m
tof p 4
T* m
T* p
T 4
T* f
tof 4 m
tof f p
tof p m
tof f 4
tof p 4
tof m p
tof m f
P f
H 4
tof 4 7
tof 4 m
tof 6 4
tof m 6
T* 4
T* m
T 6
T* 7
tof 6 4
tof 7 m
tof m 4
tof 7 6
tof m 6
tof 4 m
tof 4 7
P 7
H 6
tof i 4
tof 6 k
tof 6 4
tof i 6
tof 4 i
T* 6
T* i
T 4
T* k
tof 4 6
tof k i
tof i 6
tof k 4
tof i 4
tof 6 i
tof 6 k
P k
H 4
tof 4 q
tof 4 6
tof 3 6
tof 3 4
tof 6 3
T* 4
T* 3
T 6
T* q
tof 6 4
tof q 3
tof 3 4
tof q 6
tof 3 6
tof 4 3
tof 4 q
P q
H 6
tof 6 10
tof 6 4
tof 8 6
tof 4 8
T* 6
T* 4
T 8
T* 10
tof 8 6
tof 10 4
tof 4 6
tof 10 8
tof 4 8
tof 6 4
tof 6 10
P 10
H 8
tof 8 l
tof 8 h
tof 8 6
tof d 6
tof 1 8
tof d 1
tof 8 1
tof 6 d
tof h 8
tof l 8
T* 1
T* d
T 6
T* l
T* h
T 8
tof 6 1
tof l d
tof d 1
tof l 6
tof h l
tof h 8
tof l h
tof d 6
tof 1 d
tof 1 8
tof 1 h
tof 1 l
P h
H 6
tof 6 1
tof l 1
tof l 6
tof 1 l
T* l
T* 6
T 1
tof 1 l
tof 1 6
tof 6 1
tof l 1
P 1
H h
tof 10 d
tof h d
tof h 10
T* 10
T d
tof d 10
tof h d
tof 10 d
tof h 10
P h
P* 10
tof 10 h
tof h 10
H d
tof q 4
tof d 4
tof d q
T* q
T 4
tof 4 q
tof d 4
tof q 4
tof d q
P d
P* q
tof q d
tof d q
H 4
tof k 3
tof 4 3
tof 4 k
T* k
T 3
tof 3 k
tof 4 3
tof k 3
tof 4 k
P 4
P* k
tof k 4
tof 4 k
H 3
tof 7 i
tof 3 i
tof 3 7
T* 7
T i
tof i 7
tof 3 i
tof 7 i
tof 3 7
P 3
P* 7
tof 7 3
tof 3 7
H i
tof f m
tof i m
tof i f
T* f
T m
tof m f
tof i m
tof f m
tof i f
P i
P* f
tof f i
tof i f
H m
tof 0 p
tof m p
tof m 0
T* 0
T p
tof p 0
tof m p
tof 0 p
tof m 0
P m
P* 0
tof 0 m
tof m 0
H p
H 14
tof s 2
tof 8 p
tof s 8
tof 8 s
tof p s
tof p 8
tof 2 8
tof 14 2
tof 14 p
T s
T 14
T* 8
T* p
T* 2
tof 8 14
tof 8 p
tof p 2
tof 14 p
tof s p
tof p 2
tof 14 8
T 14
tof s 14
tof 8 14
H 8
H s
tof s i
tof s m
tof 0 s
tof m 0
T* s
T* m
T 0
T* i
tof 0 s
tof i m
tof m s
tof i 0
tof m 0
tof s m
tof s i
P i
Z 8
H 0
tof 0 3
tof 0 s
tof f 0
tof s f
T* 0
T* s
T f
T* 3
tof f 0
tof 3 s
tof s 0
tof 3 f
tof s f
tof 0 s
tof 0 3
P 3
H f
tof 7 0
tof f 4
tof f 0
tof 7 f
tof 0 7
T* f
T* 7
T 0
T* 4
tof 0 f
tof 4 7
tof 7 f
tof 4 0
tof 7 0
tof f 7
tof f 4
P 4
H 0
tof 0 d
tof 0 f
tof k f
tof k 0
tof f k
T* 0
T* k
T f
T* d
tof f 0
tof d k
tof k 0
tof d f
tof k f
tof 0 k
tof 0 d
P d
H f
tof f h
tof f 0
tof q f
tof 0 q
T* f
T* 0
T q
T* h
tof q f
tof h 0
tof 0 f
tof h q
tof 0 q
tof f 0
tof f h
P h
H q
tof q 6
tof q l
tof q f
tof 10 f
tof 10 q
tof p 10
tof 10 p
tof q 10
tof f p
tof l 10
tof 6 10
T* q
T* p
T f
T* 6
T* l
T 10
tof f q
tof 6 p
tof p q
tof 6 f
tof l 6
tof l 10
tof 6 l
tof p f
tof q p
tof q 10
tof q l
tof q 6
P l
H f
tof f q
tof 6 q
tof 6 f
tof q 6
T* 6
T* f
T q
tof q 6
tof q f
tof f q
tof 6 q
P q
H l
tof h p
tof l p
tof l h
T* h
T p
tof p h
tof l p
tof h p
tof l h
P l
P* h
tof h l
tof l h
H p
tof d 0
tof p 0
tof p d
T* d
T 0
tof 0 d
tof p 0
tof d 0
tof p d
P p
P* d
tof d p
tof p d
H 0
tof 4 k
tof 0 k
tof 0 4
T* 4
T k
tof k 4
tof 0 k
tof 4 k
tof 0 4
P 0
P* 4
tof 4 0
tof 0 4
H k
tof 3 7
tof k 7
tof k 3
T* 3
T 7
tof 7 3
tof k 7
tof 3 7
tof k 3
P k
P* 3
tof 3 k
tof k 3
H 7
tof i s
tof 7 s
tof 7 i
T* i
T s
tof s i
tof 7 s
tof i s
tof 7 i
P 7
P* i
tof i 7
tof 7 i
H s
H m
tof s m
tof 10 m
tof s 14
tof 10 14
tof 2 14
tof 10 s
tof 2 10
tof 10 2
tof s 2
tof 14 10
tof m s
tof m 14
T 2
T s
T* 10
T* m
T* 14
tof 10 s
tof 10 m
tof m 14
tof s m
tof 2 m
tof m 14
tof s 10
T s
tof 2 s
tof 10 s
H 10
H 2
tof 2 k
tof 2 7
tof i 7
tof q 10
tof i 2
tof 10 q
tof 7 i
T* 2
T* i
T 7
T* k
T* 10
tof 7 2
tof k i
tof i 2
tof k 7
tof q 10
tof i 7
tof 2 i
tof 2 k
P k
P* q
H 7
tof 3 2
tof 7 0
tof 7 2
tof 3 7
tof 2 3
T* 7
T* 3
T 2
T* 0
tof 2 7
tof 0 3
tof 3 7
tof 0 2
tof 3 2
tof 7 3
tof 7 0
P 0
H 2
tof 2 p
tof 2 7
tof 4 7
tof 4 2
tof 7 4
T* 2
T* 4
T 7
T* p
tof 7 2
tof p 4
tof 4 2
tof p 7
tof 4 7
tof 2 4
tof 2 p
P p
H 7
tof 7 l
tof 7 2
tof d 7
tof 2 d
T* 7
T* 2
T d
T* l
tof d 7
tof l 2
tof 2 7
tof l d
tof 2 d
tof 7 2
tof 7 l
P l
H d
tof d f
tof d 6
tof d 7
tof h 7
tof m d
tof h m
tof d m
tof 7 h
tof 6 d
tof f d
T* m
T* h
T 7
T* f
T* 6
T d
tof 7 m
tof f h
tof h m
tof f 7
tof 6 f
tof 6 d
tof f 6
tof h 7
tof m h
tof m d
tof m 6
tof m f
P 6
H 7
tof 7 m
tof f m
tof 7 8
tof f 7
tof m f
T* f
T* 7
T m
T* 8
tof m f
tof m 7
tof 8 m
tof 7 m
tof f m
tof m 8
P 8
H 6
tof l h
tof 6 h
tof 6 l
T* l
T h
tof h l
tof 6 h
tof l h
tof 6 l
P 6
P* l
tof l 6
tof 6 l
H h
tof p 2
tof h 2
tof h p
T* p
T 2
tof 2 p
tof h 2
tof p 2
tof h p
P h
P* p
tof p h
tof h p
H 2
tof 0 4
tof 2 4
tof 2 0
T* 0
T 4
tof 4 0
tof 2 4
tof 0 4
tof 2 0
P 2
P* 0
tof 0 2
tof 2 0
H 4
tof k 3
tof 4 3
tof 4 k
T* k
T 3
tof 3 k
tof 4 3
tof k 3
tof 4 k
P 4
P* k
tof k 4
tof 4 k
H 3
H i
tof 3 i
tof d i
tof 3 s
tof d 3
tof d 14
tof 3 d
tof s 14
tof i 3
tof i 14
T d
T 3
T* s
T* i
T* 14
tof s 3
tof s i
tof i 14
tof 3 i
tof d i
tof i 14
tof 3 s
T 3
tof d 3
tof s 3
H s
H d
tof d 2
tof d 4
tof 1 s
tof k d
tof s 1
tof 4 k
T* d
T* 4
T k
T* 2
T* s
tof k d
tof 2 4
tof 4 d
tof 2 k
tof 1 s
tof 4 k
tof d 4
tof d 2
P 2
H k
tof k h
tof k d
tof 0 d
tof 0 k
tof d 0
T* k
T* 0
T d
T* h
tof d k
tof h 0
tof 0 k
tof h d
tof 0 d
tof k 0
tof k h
P h
H d
tof d 6
tof d k
tof p d
tof k p
T* d
T* k
T p
T* 6
tof p d
tof 6 k
tof k d
tof 6 p
tof k p
tof d k
tof d 6
P 6
H p
tof p 7
tof p f
tof p d
tof l d
tof l p
tof l i
tof p i
tof d l
tof f i
tof 7 i
T* p
T* l
T d
T* 7
T* f
T i
tof d p
tof 7 l
tof l p
tof 7 d
tof f 7
tof f i
tof 7 f
tof l d
tof p l
tof p i
tof p f
tof p 7
P f
H d
tof d p
tof 7 p
tof d j
tof 7 d
tof p 7
T* 7
T* d
T p
T* j
tof p 7
tof p d
tof j p
tof d p
tof 7 p
tof p j
P j
H f
tof 6 l
tof f l
tof f 6
T* 6
T l
tof l 6
tof f l
tof 6 l
tof f 6
P f
P* 6
tof 6 f
tof f 6
H l
tof h k
tof l k
tof l h
T* h
T k
tof k h
tof l k
tof h k
tof l h
P l
P* h
tof h l
tof l h
H k
tof 2 0
tof k 0
tof k 2
T* 2
T 0
tof 0 2
tof k 0
tof 2 0
tof k 2
P k
P* 2
tof 2 k
tof k 2
H 0
H 4
tof 0 3
tof 0 i
tof 0 14
tof 3 14
tof 4 14
tof 4 i
T 0
T 4
T* 3
T* i
T* 14
tof 3 4
tof 3 i
tof i 14
tof 4 i
tof 0 i
tof i 14
tof 4 3
T 4
tof 0 4
tof 3 4
H 3
H 0
tof 0 l
tof 0 k
tof 2 k
tof 2 0
tof k 2
T* 0
T* 2
T k
T* l
tof k 0
tof l 2
tof 2 0
tof l k
tof 2 k
tof 0 2
tof 0 l
P l
Z 3
H k
tof k f
tof k 0
tof h k
tof 0 h
T* k
T* 0
T h
T* f
tof h k
tof f 0
tof 0 k
tof f h
tof 0 h
tof k 0
tof k f
P f
H h
tof h d
tof h 7
tof h k
tof 6 k
tof 6 h
tof i 6
tof 6 i
tof h 6
tof k i
tof 7 6
tof d 6
T* h
T* i
T k
T* d
T* 7
T 6
tof k h
tof d i
tof i h
tof d k
tof 7 d
tof 7 6
tof d 7
tof i k
tof h i
tof h 6
tof h 7
tof h d
P 7
H k
tof k h
tof d h
tof k c
tof d k
tof h d
T* d
T* k
T h
T* c
tof h d
tof h k
tof c h
tof k h
tof d h
tof h c
P c
H 7
tof f i
tof 7 i
tof 7 f
T* f
T i
tof i f
tof 7 i
tof f i
tof 7 f
P 7
P* f
tof f 7
tof 7 f
H i
tof l 0
tof i 0
tof i l
T* l
T 0
tof 0 l
tof i 0
tof l 0
tof i l
P i
P* l
tof l i
tof i l
H 0
H 2
tof 6 4
tof 0 6
tof 0 14
tof 4 6
tof 2 4
tof 2 14
T 0
T 2
T* 6
T* 14
T* 4
tof 6 2
tof 4 14
tof 6 14
tof 2 4
tof 0 4
tof 4 14
tof 2 6
T 2
tof 0 2
tof 6 2
H 6
H 0
tof 0 7
tof 0 i
tof 12 6
tof 13 6
tof l 0
tof 13 12
tof 6 12
tof 6 13
tof i l
T* 0
T* i
T l
T* 7
T* 13
T 6
tof l 0
tof 7 i
tof i 0
tof 7 l
tof 6 13
tof 12 6
tof 13 6
tof i l
tof 0 i
tof 0 7
tof 12 13
P 7
P* 12
P* 13
tof 12 13
H l
tof l k
tof l d
tof l 0
tof f 0
tof 4 l
tof f 4
tof l 4
tof 0 f
tof d l
tof k l
T* 4
T* f
T 0
T* k
T* d
T l
tof 0 4
tof k f
tof f 4
tof k 0
tof d k
tof d l
tof k d
tof f 0
tof 4 f
tof 4 l
tof 4 d
tof 4 k
P d
H 0
tof 0 4
tof k 4
tof 0 n
tof k 0
tof 4 k
T* k
T* 0
T 4
T* n
tof 4 k
tof 4 0
tof n 4
tof 0 4
tof k 4
tof 4 n
P n
H d
tof 7 f
tof d f
tof d 7
T* 7
T f
tof f 7
tof d f
tof 7 f
tof d 7
P d
P* 7
tof 7 d
tof d 7
H f
H i
tof 2 i
tof l i
tof f 2
tof l 2
tof 14 2
tof 14 l
tof f 14
tof 2 14
tof 2 l
tof i 2
tof i 14
T f
T 14
T* l
T* 2
T* i
tof l 14
tof i 2
tof l 2
tof 14 i
tof f i
tof i 2
tof 14 l
T 14
tof f 14
tof l 14
H l
H f
tof f 0
tof f k
tof f d
tof 7 d
tof 6 l
tof i f
tof 7 f
tof 7 i
tof f i
tof l 6
tof d 7
tof l 9
tof k i
tof 0 i
T* f
T* 7
T d
T* 0
T* k
T i
T* l
T 9
tof d f
tof 0 7
tof k 7
tof i 7
tof 7 f
tof k d
tof i 0
tof 0 d
tof i k
tof 9 l
tof 6 9
tof l 9
tof k i
tof 0 i
tof 7 d
tof 6 l
tof 6 13
P i
P* l
P* 13
tof 6 l
tof 6 13
H d
tof d k
tof 0 k
tof d 11
tof 0 d
tof k 0
T* 0
T* d
T k
T* 11
tof k 0
tof k d
tof 11 k
tof d k
tof 0 k
tof k 11
P 11
H i
H 7
tof 14 7
tof f 7
tof i 14
tof f 14
tof 2 14
tof 2 i
tof 2 f
tof i 2
tof 14 i
tof 14 f
tof 7 14
tof 7 i
T 2
T i
T* f
T* 14
T* 7
tof f i
tof 7 14
tof f 14
tof i 7
tof 2 7
tof 7 14
tof i f
T i
tof 2 i
tof f i
H f
H 2
tof 2 d
tof 0 d
tof 9 f
tof g f
tof 2 e
tof g 9
tof f 9
tof f g
tof 0 2
tof d 0
T* 0
T* 2
T d
T* g
T f
T* e
tof d 0
tof d 2
tof f g
tof e d
tof 9 f
tof g f
tof 2 d
tof 0 d
tof d e
tof 9 g
tof 9 l
tof 9 13
P e
P 9
P* g
P* l
P* 13
tof 9 g
tof 9 l
tof 9 13
H 0
tof i 0
tof 7 2
tof 0 2
tof 7 i
tof i 7
tof 2 i
tof 0 7
T 7
T* i
T* 0
T 2
tof i 7
tof 0 i
tof 0 2
tof 7 i
tof 7 2
tof i 0
H 0
tof i 7
H 2
tof f 0
tof 0 f
tof 0 b
T* 0
T b
tof b 0
tof f b
tof 0 b
tof f 0
tof f g
tof f l
tof f 13
Z f
P* 0
P* g
P* l
P* 13
P 2
tof f 0
tof f g
tof f l
tof f 13
H b
tof b 3
P b
P* 3
tof 3 b
tof b 3
H 3
tof s 1
tof 3 1
tof 3 s
T* s
T 1
tof 3 s
tof 3 1
tof s 1
tof 3 1
P 3
P* 1
tof 1 3
tof 3 1
H s
tof 10 q
tof s q
tof s 10
T* 10
T q
tof s 10
tof s q
tof 10 q
tof s q
P s
P* q
tof q s
tof s q
H 10
tof 8 m
tof 10 m
tof 10 8
T* 8
T m
tof m 8
tof 10 m
tof 8 m
tof 10 8
P 10
P* 8
tof 8 10
tof 10 8
H m
tof j p
tof m p
tof m j
T* j
T p
tof p j
tof m p
tof j p
tof m j
P m
P* j
tof j m
tof m j
H p
tof c h
tof p h
tof p c
T* c
T h
tof h c
tof p h
tof c h
tof p c
P p
P* c
tof c p
tof p c
H h
tof n 4
tof h 4
tof h n
T* n
T 4
tof 4 n
tof h 4
tof n 4
tof h n
P h
P* n
tof n h
tof h n
H 4
tof 11 k
tof 4 k
tof 4 11
T* 11
T k
tof k 11
tof 4 k
tof 11 k
tof 4 11
P 4
P* 11
tof 11 4
tof 4 11
H k
tof e d
tof k d
tof k e
T* e
T d
tof d e
tof k d
tof e d
tof k e
P k
P* e
tof e k
tof k e
H d
tof d 5
tof d a
tof d 14
tof a 14
tof 5 14
T d
T* a
T* 5
T 14
tof a d
tof 5 a
tof 5 14
tof d a
tof d 14
tof a 5
H 5
tof e k
tof 5 4
tof 5 k
tof e 5
tof k e
T* 5
T* e
T k
T* 4
tof k 5
tof 4 e
tof 4 k
tof e k
tof 5 e
tof e 4
P 4
H k
tof k h
tof k e
tof 11 e
tof 11 k
tof e 11
T* k
T* 11
T e
T* h
tof e k
tof h 11
tof h e
tof 11 e
tof k 11
tof 11 h
P h
H e
tof n 11
tof e p
tof e 11
tof n e
tof 11 n
T* e
T* n
T 11
T* p
tof 11 e
tof p n
tof p 11
tof n 11
tof e n
tof n p
P p
H 11
tof 11 m
tof 11 n
tof c n
tof c 11
tof n c
T* 11
T* c
T n
T* m
tof n 11
tof m c
tof m n
tof c n
tof 11 c
tof c m
P m
H n
tof j c
tof n 10
tof n c
tof j n
tof c j
T* n
T* j
T c
T* 10
tof c n
tof 10 j
tof 10 c
tof j c
tof n j
tof j 10
P 10
H c
tof c s
tof c j
tof 8 j
tof 8 c
tof j 8
T* c
T* 8
T j
T* s
tof j c
tof s 8
tof s j
tof 8 j
tof c 8
tof 8 s
P s
H j
tof q 8
tof j 3
tof j 8
tof q j
tof 8 q
T* j
T* q
T 8
T* 3
tof 8 j
tof 3 q
tof 3 8
tof q 8
tof j q
tof q 3
P 3
H 8
tof 8 b
tof 8 q
tof 1 q
tof 1 8
tof q 1
T* 8
T* 1
T q
T* b
tof q 8
tof b 1
tof b q
tof 1 q
tof 8 1
tof 1 b
P b
H q
tof q 1
P q
P* 1
tof 1 q
tof q 1
H 1
tof 1 f
tof 0 f
tof l 6
tof g 9
tof 1 9
tof 1 6
tof 1 12
tof l 13
tof 0 1
tof 14 13
tof g l
tof 0 g
tof 0 14
tof g 14
tof l 14
tof 13 14
tof 1 0
tof 12 14
tof 6 13
tof 9 l
tof f g
T* 0
T* 12
T 14
T 6
T 9
T* 1
T f
T g
T l
T 13
tof 14 12
tof 1 0
tof 9 6
tof f 6
tof g 6
tof l 6
tof 1 14
tof 1 6
tof f 9
tof g 9
tof l 9
tof l 1
tof g f
tof l g
tof l 13
tof f g
tof g 13
tof g l
tof 1 g
tof 9 1
tof g 13
tof g l
tof 6 13
tof 14 13
tof 12 14
tof 0 14
tof 14 13
tof 0 l
tof 0 9
tof 0 f
tof l f
tof l 9
tof l 12
tof l 6
tof l 1
tof l g
P* 12
P* 6
P* 9
P* 1
P* f
P* g
tof l 12
tof l 6
tof l 9
tof l 1
tof l f
tof l g
H l
tof l 14
P l
P* 14
tof 14 l
tof l 14
H 14
tof b q
tof 14 q
tof 14 b
T* b
T q
tof q b
tof 14 q
tof b q
tof 14 b
P 14
P* b
tof b 14
tof 14 b
H q
tof 3 8
tof q 8
tof q 3
T* 3
T 8
tof 8 3
tof q 8
tof 3 8
tof q 3
P q
P* 3
tof 3 q
tof q 3
H 8
tof s j
tof 8 j
tof 8 s
T* s
T j
tof j s
tof 8 j
tof s j
tof 8 s
P 8
P* s
tof s 8
tof 8 s
H j
tof 10 c
tof j c
tof j 10
T* 10
T c
tof c 10
tof j c
tof 10 c
tof j 10
P j
P* 10
tof 10 j
tof j 10
H c
tof m n
tof c n
tof c m
T* m
T n
tof n m
tof c n
tof m n
tof c m
P c
P* m
tof m c
tof c m
H n
tof p 11
tof n 11
tof n p
T* p
T 11
tof 11 p
tof n 11
tof p 11
tof n p
P n
P* p
tof p n
tof n p
H 11
tof h e
tof 11 e
tof 11 h
T* h
T e
tof e h
tof 11 e
tof h e
tof 11 h
P 11
P* h
tof h 11
tof 11 h
H e
tof 4 k
tof e k
tof e 4
T* 4
T k
tof k 4
tof e k
tof 4 k
tof e 4
P e
P* 4
tof 4 e
tof e 4
H k
H 5
tof k d
tof k 13
tof 5 13
T k
T 5
T* d
T* 13
tof d 5
tof 5 k
tof d 13
tof k 13
tof 5 d
T* k
T 5
tof d k
tof d 5
tof k 5
H d
H 5
tof 5 11
tof 5 e
tof 4 e
tof 4 5
tof e 4
T* 5
T* 4
T e
T* 11
tof e 5
tof 11 4
tof 11 e
tof 4 e
tof 5 4
tof 4 11
P 11
P d
H e
tof h 4
tof e n
tof e 4
tof h e
tof 4 h
T* e
T* h
T 4
T* n
tof 4 e
tof n h
tof n 4
tof h 4
tof e h
tof h n
P n
H 4
tof 4 c
tof 4 h
tof p h
tof p 4
tof h p
T* 4
T* p
T h
T* c
tof h 4
tof c p
tof c h
tof p h
tof 4 p
tof p c
P c
H h
tof m p
tof h j
tof h p
tof m h
tof p m
T* h
T* m
T p
T* j
tof p h
tof j m
tof j p
tof m p
tof h m
tof m j
P j
H p
tof p 8
tof p m
tof 10 m
tof 10 p
tof m 10
T* p
T* 10
T m
T* 8
tof m p
tof 8 10
tof 8 m
tof 10 m
tof p 10
tof 10 8
P 8
H m
tof s 10
tof m q
tof m 10
tof s m
tof 10 s
T* m
T* s
T 10
T* q
tof 10 m
tof q s
tof q 10
tof s 10
tof m s
tof s q
P q
H 10
tof 10 14
tof 10 s
tof 3 s
tof 3 10
tof s 3
T* 10
T* 3
T s
T* 14
tof s 10
tof 14 3
tof 14 s
tof 3 s
tof 10 3
tof 3 14
P 14
H s
tof b 3
tof s l
tof s 3
tof b s
tof 3 b
T* s
T* b
T 3
T* l
tof 3 s
tof l b
tof l 3
tof b 3
tof s b
tof b l
P l
H 3
tof 3 b
P 3
P* b
tof b 3
tof 3 b
H b
tof b g
tof f g
tof 9 1
tof b 1
tof b 6
tof 9 12
tof f b
tof 13 12
tof f 9
tof f 13
tof 9 13
tof 12 13
tof b f
tof 6 13
tof 1 12
tof g 9
P* f
T* f
T* 6
T 13
T 1
T* b
T g
T 9
T 12
tof 13 6
tof b f
tof g 1
tof 9 1
tof b 13
tof b 1
tof g b
tof 9 g
tof g 9
tof b 9
tof 9 12
tof 1 12
tof 13 12
tof 6 13
tof f 13
tof 13 12
tof f 9
tof f g
tof 9 g
tof 9 6
tof 9 1
tof 9 b
P* 6
P* 1
P* g
P* b
tof 9 6
tof 9 1
tof 9 g
tof 9 b
H 9
tof 9 13
P 9
P* 13
tof 13 9
tof 9 13
H 13
tof l 3
tof 13 3
tof 13 l
T* l
T 3
tof 3 l
tof 13 3
tof l 3
tof 13 l
P 13
P* l
tof l 13
tof 13 l
H 3
tof 14 s
tof 3 s
tof 3 14
T* 14
T s
tof s 14
tof 3 s
tof 14 s
tof 3 14
P 3
P* 14
tof 14 3
tof 3 14
H s
tof q 10
tof s 10
tof s q
T* q
T 10
tof 10 q
tof s 10
tof q 10
tof s q
P s
P* q
tof q s
tof s q
H 10
tof 8 m
tof 10 m
tof 10 8
T* 8
T m
tof m 8
tof 10 m
tof 8 m
tof 10 8
P 10
P* 8
tof 8 10
tof 10 8
H m
tof j p
tof m p
tof m j
T* j
T p
tof p j
tof m p
tof j p
tof m j
P m
P* j
tof j m
tof m j
H p
tof c h
tof p h
tof p c
T* c
T h
tof h c
tof p h
tof c h
tof p c
P p
P* c
tof c p
tof p c
H h
tof n 4
tof h 4
tof h n
T* n
T 4
tof 4 n
tof h 4
tof n 4
tof h n
P h
P* n
tof n h
tof h n
H 4
tof 11 e
tof 4 e
tof 4 11
T* 11
T e
tof e 11
tof 4 e
tof 11 e
tof 4 11
P 4
P* 11
tof 11 4
tof 4 11
H e
H 5
tof e k
tof e 12
tof 5 12
T e
T 5
T* k
T* 12
tof k 5
tof 5 e
tof k 12
tof e 12
tof 5 k
T* e
T 5
tof k e
tof k 5
tof e 5
H k
H 5
tof 11 4
tof 5 h
tof 5 4
tof 11 5
tof 4 11
T* 5
T* 11
T 4
T* h
tof 4 5
tof h 11
tof h 4
tof 11 4
tof 5 11
tof 11 h
P h
Z k
H 4
tof 4 p
tof 4 11
tof n 11
tof n 4
tof 11 n
T* 4
T* n
T 11
T* p
tof 11 4
tof p n
tof p 11
tof n 11
tof 4 n
tof n p
P p
H 11
tof c n
tof 11 m
tof 11 n
tof c 11
tof n c
T* 11
T* c
T n
T* m
tof n 11
tof m c
tof m n
tof c n
tof 11 c
tof c m
P m
H n
tof n 10
tof n c
tof j c
tof j n
tof c j
T* n
T* j
T c
T* 10
tof c n
tof 10 j
tof 10 c
tof j c
tof n j
tof j 10
P 10
H c
tof 8 j
tof c s
tof c j
tof 8 c
tof j 8
T* c
T* 8
T j
T* s
tof j c
tof s 8
tof s j
tof 8 j
tof c 8
tof 8 s
P s
H j
tof j 3
tof j 8
tof q 8
tof q j
tof 8 q
T* j
T* q
T 8
T* 3
tof 8 j
tof 3 q
tof 3 8
tof q 8
tof j q
tof q 3
P 3
H 8
tof 14 q
tof 8 13
tof 8 q
tof 14 8
tof q 14
T* 8
T* 14
T q
T* 13
tof q 8
tof 13 14
tof 13 q
tof 14 q
tof 8 14
tof 14 13
P 13
H q
tof q 9
tof q 14
tof l 14
tof l q
tof 14 l
T* q
T* l
T 14
T* 9
tof 14 q
tof 9 l
tof 9 14
tof l 14
tof q l
tof l 9
P 9
H 14
tof 14 l
P 14
P* l
tof l 14
tof 14 l
H l
tof l b
tof g b
tof l 1
tof 12 6
tof g l
tof g 6
tof g 12
tof 6 12
tof l g
tof 1 12
tof b 6
P g
T g
T* 1
T 12
T* l
T b
T 6
tof 12 1
tof l g
tof l 12
tof b l
tof b 6
tof 12 6
tof l b
tof 1 12
tof g 12
tof 12 6
tof g b
tof b 1
tof b l
P* 1
P* l
tof b 1
tof b l
H b
tof b 12
P b
P* 12
tof 12 b
tof b 12
H 12
tof 9 14
tof 12 14
tof 12 9
T* 9
T 14
tof 14 9
tof 12 14
tof 9 14
tof 12 9
P 12
P* 9
tof 9 12
tof 12 9
H 14
tof 13 q
tof 14 q
tof 14 13
T* 13
T q
tof q 13
tof 14 q
tof 13 q
tof 14 13
P 14
P* 13
tof 13 14
tof 14 13
H q
tof 3 8
tof q 8
tof q 3
T* 3
T 8
tof 8 3
tof q 8
tof 3 8
tof q 3
P q
P* 3
tof 3 q
tof q 3
H 8
tof s j
tof 8 j
tof 8 s
T* s
T j
tof j s
tof 8 j
tof s j
tof 8 s
P 8
P* s
tof s 8
tof 8 s
H j
tof 10 c
tof j c
tof j 10
T* 10
T c
tof c 10
tof j c
tof 10 c
tof j 10
P j
P* 10
tof 10 j
tof j 10
H c
tof m n
tof c n
tof c m
T* m
T n
tof n m
tof c n
tof m n
tof c m
P c
P* m
tof m c
tof c m
H n
tof p 11
tof n 11
tof n p
T* p
T 11
tof 11 p
tof n 11
tof p 11
tof n p
P n
P* p
tof p n
tof n p
H 11
tof h 4
tof 11 4
tof 11 h
T* h
T 4
tof 4 h
tof 11 4
tof h 4
tof 11 h
P 11
P* h
tof h 11
tof 11 h
H 4
H 5
tof 4 e
tof 4 6
tof 5 6
T 4
T 5
T* e
T* 6
tof e 5
tof 5 4
tof e 6
tof 4 6
tof 5 e
T* 4
T 5
tof e 4
tof e 5
tof 4 5
H e
H 5
tof 5 n
tof 5 11
tof h 11
tof h 5
tof 11 h
T* 5
T* h
T 11
T* n
tof 11 5
tof n h
tof n 11
tof h 11
tof 5 h
tof h n
P n
P* e
H 11
tof p h
tof 11 c
tof 11 h
tof p 11
tof h p
T* 11
T* p
T h
T* c
tof h 11
tof c p
tof c h
tof p h
tof 11 p
tof p c
P c
H h
tof h j
tof h p
tof m p
tof m h
tof p m
T* h
T* m
T p
T* j
tof p h
tof j m
tof j p
tof m p
tof h m
tof m j
P j
H p
tof 10 m
tof p 8
tof p m
tof 10 p
tof m 10
T* p
T* 10
T m
T* 8
tof m p
tof 8 10
tof 8 m
tof 10 m
tof p 10
tof 10 8
P 8
H m
tof m q
tof m 10
tof s 10
tof s m
tof 10 s
T* m
T* s
T 10
T* q
tof 10 m
tof q s
tof q 10
tof s 10
tof m s
tof s q
P q
H 10
tof 3 s
tof 10 14
tof 10 s
tof 3 10
tof s 3
T* 10
T* 3
T s
T* 14
tof s 10
tof 14 3
tof 14 s
tof 3 s
tof 10 3
tof 3 14
P 14
H s
tof s 12
tof s 3
tof 13 3
tof 13 s
tof 3 13
T* s
T* 13
T 3
T* 12
tof 3 s
tof 12 13
tof 12 3
tof 13 3
tof s 13
tof 13 12
P 12
H 3
tof 3 b
tof 3 13
tof 9 3
tof 13 9
T* 3
T* 13
T 9
T* b
tof 9 3
tof b 13
tof b 9
tof 13 9
tof 3 13
tof 13 b
P b
H 9
tof 9 13
P 9
P* 13
tof 13 9
tof 9 13
H 13
tof 13 l
tof 1 l
tof 6 l
tof 6 1
tof 1 6
tof 13 6
tof l 1
T 13
T* 1
T* 6
T l
tof 6 13
tof 6 1
tof 6 l
tof 13 1
tof 1 l
tof 13 6
H 6
tof 6 1
P 6
P* 1
tof 1 6
tof 6 1
H 1
tof b 9
tof 1 9
tof 1 b
T* b
T 9
tof 9 b
tof 1 9
tof b 9
tof 1 b
P 1
P* b
tof b 1
tof 1 b
H 9
tof 12 3
tof 9 3
tof 9 12
T* 12
T 3
tof 3 12
tof 9 3
tof 12 3
tof 9 12
P 9
P* 12
tof 12 9
tof 9 12
H 3
tof 14 s
tof 3 s
tof 3 14
T* 14
T s
tof s 14
tof 3 s
tof 14 s
tof 3 14
P 3
P* 14
tof 14 3
tof 3 14
H s
tof q 10
tof s 10
tof s q
T* q
T 10
tof 10 q
tof s 10
tof q 10
tof s q
P s
P* q
tof q s
tof s q
H 10
tof 8 m
tof 10 m
tof 10 8
T* 8
T m
tof m 8
tof 10 m
tof 8 m
tof 10 8
P 10
P* 8
tof 8 10
tof 10 8
H m
tof j p
tof m p
tof m j
T* j
T p
tof p j
tof m p
tof j p
tof m j
P m
P* j
tof j m
tof m j
H p
tof c h
tof p h
tof p c
T* c
T h
tof h c
tof p h
tof c h
tof p c
P p
P* c
tof c p
tof p c
H h
tof n 11
tof h 11
tof h n
T* n
T 11
tof 11 n
tof h 11
tof n 11
tof h n
P h
P* n
tof n h
tof h n
H 11
H 5
tof 11 4
tof 11 l
tof 5 l
T 11
T 5
T* 4
T* l
tof 4 5
tof 5 11
tof 4 l
tof 11 l
tof 5 4
T* 11
T 5
tof 4 11
tof 4 5
tof 11 5
H 4
H 5
tof n h
tof 5 p
tof 5 h
tof n 5
tof h n
T* 5
T* n
T h
T* p
tof h 5
tof p n
tof p h
tof n h
tof 5 n
tof n p
P p
H h
tof h m
tof h n
tof c n
tof c h
tof n c
T* h
T* c
T n
T* m
tof n h
tof m c
tof m n
tof c n
tof h c
tof c m
P m
H n
tof j c
tof n 10
tof n c
tof j n
tof c j
T* n
T* j
T c
T* 10
tof c n
tof 10 j
tof 10 c
tof j c
tof n j
tof j 10
P 10
H c
tof c s
tof c j
tof 8 j
tof 8 c
tof j 8
T* c
T* 8
T j
T* s
tof j c
tof s 8
tof s j
tof 8 j
tof c 8
tof 8 s
P s
H j
tof q 8
tof j 3
tof j 8
tof q j
tof 8 q
T* j
T* q
T 8
T* 3
tof 8 j
tof 3 q
tof 3 8
tof q 8
tof j q
tof q 3
P 3
H 8
tof 8 9
tof 8 q
tof 14 q
tof 14 8
tof q 14
T* 8
T* 14
T q
T* 9
tof q 8
tof 9 14
tof 9 q
tof 14 q
tof 8 14
tof 14 9
P 9
H q
tof q 1
tof q 14
tof 12 q
tof 14 12
T* q
T* 14
T 12
T* 1
tof 12 q
tof 1 14
tof 1 12
tof 14 12
tof q 14
tof 14 1
P 1
H 12
tof 12 6
tof 12 14
tof b 12
tof 14 b
T* 12
T* 14
T b
T* 6
tof b 12
tof 6 14
tof 6 b
tof 14 b
tof 12 14
tof 14 6
P 6
H b
tof b t
tof b 7
tof l b
tof b l
tof 7 b
tof t b
P l
T l
T* t
T* 7
T b
tof t l
tof 7 t
tof 7 b
tof l t
tof l b
tof t 7
tof 7 14
P* 14
tof 7 14
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
tof 1 12
tof 14 12
tof 14 1
T* 1
T 12
tof 12 1
tof 14 12
tof 1 12
tof 14 1
P 14
P* 1
tof 1 14
tof 14 1
H 12
tof 9 q
tof 12 q
tof 12 9
T* 9
T q
tof q 9
tof 12 q
tof 9 q
tof 12 9
P 12
P* 9
tof 9 12
tof 12 9
H q
tof 3 8
tof q 8
tof q 3
T* 3
T 8
tof 8 3
tof q 8
tof 3 8
tof q 3
P q
P* 3
tof 3 q
tof q 3
H 8
tof s j
tof 8 j
tof 8 s
T* s
T j
tof j s
tof 8 j
tof s j
tof 8 s
P 8
P* s
tof s 8
tof 8 s
H j
tof 10 c
tof j c
tof j 10
T* 10
T c
tof c 10
tof j c
tof 10 c
tof j 10
P j
P* 10
tof 10 j
tof j 10
H c
tof m n
tof c n
tof c m
T* m
T n
tof n m
tof c n
tof m n
tof c m
P c
P* m
tof m c
tof c m
H n
tof p h
tof n h
tof n p
T* p
T h
tof h p
tof n h
tof p h
tof n p
P n
P* p
tof p n
tof n p
H h
H 5
tof h 11
tof h b
tof 5 b
T h
T 5
T* 11
T* b
tof 11 5
tof 5 h
tof 11 b
tof h b
tof 5 11
T* h
T 5
tof 11 h
tof 11 5
tof h 5
H 11
H 5
tof 5 c
tof 5 n
tof p n
tof p 5
tof n p
T* 5
T* p
T n
T* c
tof n 5
tof c p
tof c n
tof p n
tof 5 p
tof p c
P c
P 11
H n
tof m p
tof n j
tof n p
tof m n
tof p m
T* n
T* m
T p
T* j
tof p n
tof j m
tof j p
tof m p
tof n m
tof m j
P j
H p
tof p 8
tof p m
tof 10 m
tof 10 p
tof m 10
T* p
T* 10
T m
T* 8
tof m p
tof 8 10
tof 8 m
tof 10 m
tof p 10
tof 10 8
P 8
H m
tof s 10
tof m q
tof m 10
tof s m
tof 10 s
T* m
T* s
T 10
T* q
tof 10 m
tof q s
tof q 10
tof s 10
tof m s
tof s q
P q
H 10
tof 10 12
tof 10 s
tof 3 s
tof 3 10
tof s 3
T* 10
T* 3
T s
T* 12
tof s 10
tof 12 3
tof 12 s
tof 3 s
tof 10 3
tof 3 12
P 12
H s
tof s 14
tof s 3
tof 9 s
tof 3 9
T* s
T* 3
T 9
T* 14
tof 9 s
tof 14 3
tof 14 9
tof 3 9
tof s 3
tof 3 14
P 14
H 9
tof 9 7
tof 9 3
tof 1 9
tof 3 1
T* 9
T* 3
T 1
T* 7
tof 1 9
tof 7 3
tof 7 1
tof 3 1
tof 9 3
tof 3 7
P 7
H 1
tof 1 l
tof 6 3
tof 1 t
tof 1 3
tof b 1
tof 6 b
tof 1 b
tof t 1
tof 3 6
tof l 1
T* b
T* 6
T 3
T* l
T* t
T 1
tof 3 b
tof l 6
tof 6 b
tof l 3
tof t l
tof t 1
tof l t
tof 6 3
tof b 6
tof b 1
tof b t
tof b l
P t
H 3
tof 3 b
tof 3 l
tof 3 1
tof l 1
tof b 1
T 3
T* b
T* l
T 1
tof b 3
tof l b
tof l 1
tof 3 b
tof 3 1
tof b l
H l
H t
tof 7 6
tof t 6
tof t 7
tof 0 l
tof 0 13
tof 0 g
tof 0 f
tof 0 2
tof l 13
tof l f
tof l g
tof l 2
tof l 0
T* 7
T 6
P* 0
T* 0
T* 2
T* l
T* f
T* g
T* 13
tof 6 7
tof t 6
tof l 0
tof l 2
tof f l
tof g f
tof 13 g
tof g 13
tof l f
tof f 13
tof f g
tof 0 2
tof 0 13
tof 0 g
tof 0 f
tof 0 l
tof 7 6
tof t 7
P t
P* 7
tof 7 t
tof t 7
H 6
tof 14 9
tof 6 9
tof 6 14
T* 14
T 9
tof 9 14
tof 6 9
tof 14 9
tof 6 14
P 6
P* 14
tof 14 6
tof 6 14
H 9
tof 12 s
tof 9 s
tof 9 12
T* 12
T s
tof s 12
tof 9 s
tof 12 s
tof 9 12
P 9
P* 12
tof 12 9
tof 9 12
H s
tof q 10
tof s 10
tof s q
T* q
T 10
tof 10 q
tof s 10
tof q 10
tof s q
P s
P* q
tof q s
tof s q
H 10
tof 8 m
tof 10 m
tof 10 8
T* 8
T m
tof m 8
tof 10 m
tof 8 m
tof 10 8
P 10
P* 8
tof 8 10
tof 10 8
H m
tof j p
tof m p
tof m j
T* j
T p
tof p j
tof m p
tof j p
tof m j
P m
P* j
tof j m
tof m j
H p
tof c n
tof p n
tof p c
T* c
T n
tof n c
tof p n
tof c n
tof p c
P p
P* c
tof c p
tof p c
H n
H 5
tof n h
tof 1 n
tof n 1
tof 5 n
T 1
T 5
T* h
T* n
tof h 5
tof 5 1
tof h n
tof 1 n
tof 5 h
T* 1
T 5
tof h 1
tof h 5
tof 1 5
H h
H 5
tof c p
tof 5 m
tof 5 p
tof c 5
tof p c
T* 5
T* c
T p
T* m
tof p 5
tof m c
tof m p
tof c p
tof 5 c
tof c m
P m
Z h
H p
tof p 10
tof p c
tof j c
tof j p
tof c j
T* p
T* j
T c
T* 10
tof c p
tof 10 j
tof 10 c
tof j c
tof p j
tof j 10
P 10
H c
tof 8 j
tof c s
tof c j
tof 8 c
tof j 8
T* c
T* 8
T j
T* s
tof j c
tof s 8
tof s j
tof 8 j
tof c 8
tof 8 s
P s
H j
tof j 9
tof j 8
tof q 8
tof q j
tof 8 q
T* j
T* q
T 8
T* 9
tof 8 j
tof 9 q
tof 9 8
tof q 8
tof j q
tof q 9
P 9
H 8
tof 8 6
tof 8 q
tof 12 8
tof q 12
T* 8
T* q
T 12
T* 6
tof 12 8
tof 6 q
tof 6 12
tof q 12
tof 8 q
tof q 6
P 6
H 12
tof 12 t
tof 12 q
tof 14 12
tof q 14
T* 12
T* q
T 14
T* t
tof 14 12
tof t q
tof t 14
tof q 14
tof 12 q
tof q t
P t
H 14
tof 14 3
tof 7 q
tof 14 b
tof 14 q
tof n 14
tof 7 n
tof 14 n
tof b 14
tof q 7
tof 3 14
T* n
T* 7
T q
T* 3
T* b
T 14
tof q n
tof 3 7
tof 7 n
tof 3 q
tof b 3
tof b 14
tof 3 b
tof 7 q
tof n 7
tof n 14
tof n b
tof n 3
P b
H q
tof q n
tof 3 q
tof n 3
T* n
T* q
T 3
tof 3 n
tof 3 q
tof q 3
tof n 3
P 3
H b
tof t 7
tof b 7
tof b t
T* t
T 7
tof 7 t
tof b 7
tof t 7
tof b t
P b
P* t
tof t b
tof b t
H 7
tof 6 12
tof 7 12
tof 7 6
T* 6
T 12
tof 12 6
tof 7 12
tof 6 12
tof 7 6
P 7
P* 6
tof 6 7
tof 7 6
H 12
tof 9 8
tof 12 8
tof 12 9
T* 9
T 8
tof 8 9
tof 12 8
tof 9 8
tof 12 9
P 12
P* 9
tof 9 12
tof 12 9
H 8
tof s j
tof 8 j
tof 8 s
T* s
T j
tof j s
tof 8 j
tof s j
tof 8 s
P 8
P* s
tof s 8
tof 8 s
H j
tof 10 c
tof j c
tof j 10
T* 10
T c
tof c 10
tof j c
tof 10 c
tof j 10
P j
P* 10
tof 10 j
tof j 10
H c
tof m p
tof c p
tof c m
T* m
T p
tof p m
tof c p
tof m p
tof c m
P c
P* m
tof m c
tof c m
H p
H 5
tof p 1
tof 14 p
tof p 14
tof 5 p
T 14
T 5
T* 1
T* p
tof 1 5
tof 5 14
tof 1 p
tof 14 p
tof 5 1
T* 14
T 5
tof 1 14
tof 1 5
tof 14 5
H 1
H 5
tof 5 j
tof 5 c
tof m c
tof m 5
tof c m
T* 5
T* m
T c
T* j
tof c 5
tof j m
tof j c
tof m c
tof 5 m
tof m j
P j
P* 1
H c
tof 10 m
tof c 8
tof c m
tof 10 c
tof m 10
T* c
T* 10
T m
T* 8
tof m c
tof 8 10
tof 8 m
tof 10 m
tof c 10
tof 10 8
P 8
H m
tof m 12
tof m 10
tof s 10
tof s m
tof 10 s
T* m
T* s
T 10
T* 12
tof 10 m
tof 12 s
tof 12 10
tof s 10
tof m s
tof s 12
P 12
H 10
tof 10 7
tof 10 s
tof 9 10
tof s 9
T* 10
T* s
T 9
T* 7
tof 9 10
tof 7 s
tof 7 9
tof s 9
tof 10 s
tof s 7
P 7
H 9
tof 9 b
tof 9 s
tof 6 9
tof s 6
T* 9
T* s
T 6
T* b
tof 6 9
tof b s
tof b 6
tof s 6
tof 9 s
tof s b
P b
H 6
tof 6 q
tof t s
tof 6 n
tof 6 s
tof p 6
tof t p
tof 6 p
tof n 6
tof s t
tof q 6
T* p
T* t
T s
T* q
T* n
T 6
tof s p
tof q t
tof t p
tof q s
tof n q
tof n 6
tof q n
tof t s
tof p t
tof p 6
tof p n
tof p q
P n
H s
tof s p
tof s 1
tof q s
tof p q
T* p
T* s
T q
T* 1
tof q p
tof q s
tof 1 q
tof s q
tof p q
tof q 1
P 1
H n
tof b t
tof n t
tof n b
T* b
T t
tof t b
tof n t
tof b t
tof n b
P n
P* b
tof b n
tof n b
H t
tof 7 9
tof t 9
tof t 7
T* 7
T 9
tof 9 7
tof t 9
tof 7 9
tof t 7
P t
P* 7
tof 7 t
tof t 7
H 9
tof 12 10
tof 9 10
tof 9 12
T* 12
T 10
tof 10 12
tof 9 10
tof 12 10
tof 9 12
P 9
P* 12
tof 12 9
tof 9 12
H 10
tof 8 m
tof 10 m
tof 10 8
T* 8
T m
tof m 8
tof 10 m
tof 8 m
tof 10 8
P 10
P* 8
tof 8 10
tof 10 8
H m
tof j c
tof m c
tof m j
T* j
T c
tof c j
tof m c
tof j c
tof m j
P m
P* j
tof j m
tof m j
H c
H 5
tof c 14
tof 6 c
tof c 6
tof 5 c
T 6
T 5
T* 14
T* c
tof 14 5
tof 5 6
tof 14 c
tof 6 c
tof 5 14
T* 6
T 5
tof 14 6
tof 14 5
tof 6 5
H 14
H 5
tof 5 10
tof 5 m
tof 3 14
tof j 5
tof 14 3
tof m j
T* 5
T* m
T j
T* 10
T* 14
tof j 5
tof 10 m
tof 10 j
tof 3 14
tof m j
tof 5 m
tof m 10
P 10
H j
tof j 9
tof j m
tof 8 m
tof 8 j
tof m 8
T* j
T* 8
T m
T* 9
tof m j
tof 9 8
tof 9 m
tof 8 m
tof j 8
tof 8 9
P 9
H m
tof m t
tof m 8
tof 12 m
tof 8 12
T* m
T* 8
T 12
T* t
tof 12 m
tof t 8
tof t 12
tof 8 12
tof m 8
tof 8 t
P t
H 12
tof 12 n
tof 12 8
tof 7 12
tof 8 7
T* 12
T* 8
T 7
T* n
tof 7 12
tof n 8
tof n 7
tof 8 7
tof 12 8
tof 8 n
P n
H 7
tof 7 s
tof b 8
tof 7 p
tof 7 8
tof c 7
tof b c
tof 7 c
tof p 7
tof 8 b
tof s 7
T* c
T* b
T 8
T* s
T* p
T 7
tof 8 c
tof s b
tof b c
tof s 8
tof p s
tof p 7
tof s p
tof b 8
tof c b
tof c 7
tof c p
tof c s
P p
H 8
tof 8 c
tof 8 h
tof s 8
tof c s
T* c
T* 8
T s
T* h
tof s c
tof s 8
tof h s
tof 8 s
tof c s
tof s h
P h
H p
tof n b
tof p b
tof p n
T* n
T b
tof b n
tof p b
tof n b
tof p n
P p
P* n
tof n p
tof p n
H b
tof t 12
tof b 12
tof b t
T* t
T 12
tof 12 t
tof b 12
tof t 12
tof b t
P b
P* t
tof t b
tof b t
H 12
tof 9 m
tof 12 m
tof 12 9
T* 9
T m
tof m 9
tof 12 m
tof 9 m
tof 12 9
P 12
P* 9
tof 9 12
tof 12 9
H m
tof 10 j
tof m j
tof m 10
T* 10
T j
tof j 10
tof m j
tof 10 j
tof m 10
P m
P* 10
tof 10 m
tof m 10
H j
H 5
tof j 6
tof j 7
tof 5 7
T j
T 5
T* 6
T* 7
tof 6 5
tof 5 j
tof 6 7
tof j 7
tof 5 6
T* j
T 5
tof 6 j
tof 6 5
tof j 5
H 6
H 5
tof 5 12
tof 5 m
tof 10 m
tof 10 5
tof m 10
T* 5
T* 10
T m
T* 12
tof m 5
tof 12 10
tof 12 m
tof 10 m
tof 5 10
tof 10 12
P 12
Z 6
H m
tof m b
tof m 10
tof 9 m
tof 10 9
T* m
T* 10
T 9
T* b
tof 9 m
tof b 10
tof b 9
tof 10 9
tof m 10
tof 10 b
P b
H 9
tof 9 p
tof 9 10
tof t 9
tof 10 t
T* 9
T* 10
T t
T* p
tof t 9
tof p 10
tof p t
tof 10 t
tof 9 10
tof 10 p
P p
H t
tof t 8
tof n 10
tof t c
tof t 10
tof n t
tof n 7
tof t 7
tof c 7
tof 10 n
tof 8 7
T* t
T* n
T 10
T* 8
T* c
T 7
tof 10 t
tof 8 n
tof n t
tof 8 10
tof c 8
tof c 7
tof 8 c
tof n 10
tof t n
tof t 7
tof t c
tof t 8
P c
H 10
tof 10 t
tof 10 11
tof 8 10
tof t 8
T* t
T* 10
T 8
T* 11
tof 8 t
tof 8 10
tof 11 8
tof 10 8
tof t 8
tof 8 11
P 11
H c
tof p n
tof c n
tof c p
T* p
T n
tof n p
tof c n
tof p n
tof c p
P c
P* p
tof p c
tof c p
H n
tof b 9
tof n 9
tof n b
T* b
T 9
tof 9 b
tof n 9
tof b 9
tof n b
P n
P* b
tof b n
tof n b
H 9
tof 12 m
tof 9 m
tof 9 12
T* 12
T m
tof m 12
tof 9 m
tof 12 m
tof 9 12
P 9
P* 12
tof 12 9
tof 9 12
H m
H 5
tof m j
tof 7 m
tof m 7
tof 5 m
T 7
T 5
T* j
T* m
tof j 5
tof 5 7
tof j m
tof 7 m
tof 5 j
T* 7
T 5
tof j 7
tof j 5
tof 7 5
H j
H 5
tof 5 n
tof 5 9
tof 13 j
tof g j
tof 12 5
tof g 13
tof j 13
tof j g
tof 9 12
T* 5
T* 9
T 12
T* n
T* g
T j
tof 12 5
tof n 9
tof n 12
tof j g
tof 13 j
tof g j
tof 9 12
tof 5 9
tof 9 n
tof 13 g
P n
P* 13
P* g
tof 13 g
H 12
tof 12 c
tof 12 9
tof b 12
tof 9 b
T* 12
T* 9
T b
T* c
tof b 12
tof c 9
tof c b
tof 9 b
tof 12 9
tof 9 c
P c
H b
tof b 10
tof p 9
tof b t
tof b 9
tof m b
tof p m
tof b m
tof t b
tof 9 p
tof 10 b
T* m
T* p
T 9
T* 10
T* t
T b
tof 9 m
tof 10 p
tof p m
tof 10 9
tof t 10
tof t b
tof 10 t
tof p 9
tof m p
tof m b
tof m t
tof m 10
P t
H 9
tof 9 m
tof 9 4
tof 10 9
tof m 10
T* m
T* 9
T 10
T* 4
tof 10 m
tof 10 9
tof 4 10
tof 9 10
tof m 10
tof 10 4
P 4
H t
tof c p
tof t p
tof t c
T* c
T p
tof p c
tof t p
tof c p
tof t c
P t
P* c
tof c t
tof t c
H p
tof n 12
tof p 12
tof p n
T* n
T 12
tof 12 n
tof p 12
tof n 12
tof p n
P p
P* n
tof n p
tof p n
H 12
H 5
tof 12 7
tof b 12
tof 12 b
tof 5 12
T b
T 5
T* 7
T* 12
tof 7 5
tof 5 b
tof 7 12
tof b 12
tof 5 7
T* b
T 5
tof 7 b
tof 7 5
tof b 5
H 7
H 5
tof 5 t
tof 5 p
tof j 7
tof n 5
tof 7 j
tof p n
tof 7 f
T* 5
T* p
T n
T* t
T* 7
T f
tof n 5
tof t p
tof t n
tof f 7
tof j f
tof 7 f
tof p n
tof 5 p
tof p t
tof j 7
tof j g
P t
P* 7
P* g
tof j 7
tof j g
H n
tof n 9
tof c p
tof n m
tof n p
tof c n
tof c 12
tof n 12
tof m 12
tof p c
tof 9 12
T* n
T* c
T p
T* 9
T* m
T 12
tof p n
tof 9 c
tof c n
tof 9 p
tof m 9
tof m 12
tof 9 m
tof c p
tof n c
tof n 12
tof n m
tof n 9
P m
H p
tof p n
tof p e
tof 9 p
tof n 9
T* n
T* p
T 9
T* e
tof 9 n
tof 9 p
tof e 9
tof p 9
tof n 9
tof 9 e
P e
H m
tof t c
tof m c
tof m t
T* t
T c
tof c t
tof m c
tof t c
tof m t
P m
P* t
tof t m
tof m t
H c
H 5
tof c b
tof 12 c
tof c 12
tof 5 c
T 12
T 5
T* b
T* c
tof b 5
tof 5 12
tof b c
tof 12 c
tof 5 b
T* 12
T 5
tof b 12
tof b 5
tof 12 5
H b
H 5
tof 5 p
tof t m
tof c m
tof 5 n
tof 5 m
tof f b
tof t 5
tof c t
tof t c
tof 5 c
tof b f
tof b l
tof n c
tof m t
tof p c
T* 5
T* t
T m
T* p
T* n
T c
T* b
T l
tof n p
tof c p
tof p 5
tof p t
tof p m
tof n p
tof c n
tof l b
tof f l
tof b l
tof n c
tof p c
tof t m
tof 5 m
tof f b
tof f 7
tof f g
P c
P f
P* b
P* 7
P* g
tof f b
tof f 7
tof f g
H 5
tof 5 n
tof 5 k
tof p 5
tof n p
T* n
T* 5
T p
T* k
tof p n
tof p 5
tof k p
tof 5 p
tof n p
tof p k
P k
H c
H t
tof c t
tof m t
tof c 12
tof m c
tof c m
tof t c
T m
T c
T* 12
T* t
tof 12 c
tof c m
tof 12 t
tof m t
tof c 12
T* m
T c
tof 12 m
tof 12 c
tof m c
H 12
H c
tof c 5
tof l 12
tof c d
tof 12 l
tof n c
tof 12 0
tof 5 n
T* 5
T* c
T n
T* 12
T 0
T* d
tof n 5
tof n c
tof 0 12
tof d n
tof l 0
tof 12 0
tof c n
tof 5 n
tof n d
tof l 12
tof l b
tof l 7
tof l g
P d
Z l
P* 12
P* b
P* 7
P* g
tof l 12
tof l b
tof l 7
tof l g
H 5
tof 5 c
tof m 5
tof m t
tof 5 m
tof c t
T m
T* c
T* 5
T t
tof c m
tof 5 c
tof 5 t
tof m c
tof m t
tof c 5
H 5
tof 0 5
tof 5 0
tof 5 2
T* 5
T 2
tof 2 5
tof 0 2
tof 5 2
tof 0 12
tof 0 b
tof 0 7
tof 0 g
tof 0 5
P* 0
P* 5
P* 12
P* b
P* 7
P* g
tof c m
tof 0 5
tof 0 12
tof 0 b
tof 0 7
tof 0 g
H 2
tof 2 6
P 2
P* 6
tof 6 2
tof 2 6
H 6
tof 14 3
tof 6 3
tof 6 14
T* 14
T 3
tof 6 14
tof 6 3
tof 14 3
tof 6 3
P 6
P* 3
tof 3 6
tof 6 3
H 14
tof 1 q
tof 14 q
tof 14 1
T* 1
T q
tof q 1
tof 14 q
tof 1 q
tof 14 1
P 14
P* 1
tof 1 14
tof 14 1
H q
tof h s
tof q s
tof q h
T* h
T s
tof s h
tof q s
tof h s
tof q h
P q
P* h
tof h q
tof q h
H s
tof 11 8
tof s 8
tof s 11
T* 11
T 8
tof 8 11
tof s 8
tof 11 8
tof s 11
P s
P* 11
tof 11 s
tof s 11
H 8
tof 4 10
tof 8 10
tof 8 4
T* 4
T 10
tof 10 4
tof 8 10
tof 4 10
tof 8 4
P 8
P* 4
tof 4 8
tof 8 4
H 10
tof e 9
tof 10 9
tof 10 e
T* e
T 9
tof 9 e
tof 10 9
tof e 9
tof 10 e
P 10
P* e
tof e 10
tof 10 e
H 9
tof k p
tof 9 p
tof 9 k
T* k
T p
tof p k
tof 9 p
tof k p
tof 9 k
P 9
P* k
tof k 9
tof 9 k
H p
tof d n
tof p n
tof p d
T* d
T n
tof n d
tof p n
tof d n
tof p d
P p
P* d
tof d p
tof p d
H n
tof n o
tof t a
tof o a
tof n t
tof o t
tof a t
T n
T* o
T* t
T a
tof o n
tof t o
tof t a
tof n o
tof n a
tof o t
H o
H t
tof t 9
tof t p
tof d t
tof p d
T* t
T* p
T d
T* 9
tof d t
tof 9 p
tof 9 d
tof p d
tof t p
tof p 9
P 9
H d
tof d 10
tof d p
tof k d
tof p k
T* d
T* p
T k
T* 10
tof k d
tof 10 p
tof 10 k
tof p k
tof d p
tof p 10
P 10
H k
tof k 8
tof k p
tof e k
tof p e
T* k
T* p
T e
T* 8
tof e k
tof 8 p
tof 8 e
tof p e
tof k p
tof p 8
P 8
H e
tof e s
tof e p
tof 4 e
tof p 4
T* e
T* p
T 4
T* s
tof 4 e
tof s p
tof s 4
tof p 4
tof e p
tof p s
P s
H 4
tof 4 q
tof 4 p
tof 11 4
tof p 11
T* 4
T* p
T 11
T* q
tof 11 4
tof q p
tof q 11
tof p 11
tof 4 p
tof p q
P q
H 11
tof 11 14
tof 11 p
tof h 11
tof p h
T* 11
T* p
T h
T* 14
tof h 11
tof 14 p
tof 14 h
tof p h
tof 11 p
tof p 14
P 14
H h
tof h 6
tof h p
tof 1 h
tof p 1
T* h
T* p
T 1
T* 6
tof 1 h
tof 6 p
tof 6 1
tof p 1
tof h p
tof p 6
P 6
H 1
tof 1 2
tof 1 p
tof 3 1
tof p 3
T* 1
T* p
T 3
T* 2
tof 3 1
tof 2 p
tof 2 3
tof p 3
tof 1 p
tof p 2
P 2
H 3
tof 3 p
P 3
P* p
tof p 3
tof 3 p
H p
tof 12 l
tof 5 0
tof p 0
tof p l
tof 7 j
tof b f
tof p f
tof p j
tof p 13
tof g 13
tof 7 g
tof 5 p
tof a g
tof b 7
tof 12 b
tof 5 12
tof a 5
tof 5 a
tof 12 5
tof b 5
tof 7 5
tof p a
tof g 5
tof 13 5
tof j g
tof f 7
tof l b
tof 0 12
T a
T* 5
T 13
T j
T f
T l
T* p
T 0
T 12
T b
T 7
T g
tof 13 5
tof p a
tof 0 a
tof 12 a
tof f j
tof l j
tof 0 p
tof 12 p
tof b j
tof 7 j
tof p 13
tof p j
tof b l
tof 7 f
tof l f
tof 12 0
tof 7 p
tof 7 12
tof 7 b
tof 7 g
tof p g
tof p 7
tof p 12
tof p b
tof l b
tof f p
tof j g
tof 13 g
tof 5 13
tof a 13
tof 13 g
tof a 7
tof a 0
tof a f
tof a l
tof 7 a
tof 7 l
tof 7 f
tof 7 5
tof 7 j
tof 7 p
tof 7 b
tof 7 12
P* 5
P* j
P* f
P* p
P* l
P* b
P* a
P* 12
tof 7 5
tof 7 j
tof 7 f
tof 7 p
tof 7 l
tof 7 b
tof 7 a
tof 7 12
H 7
tof 7 13
P 7
P* 13
tof 13 7
tof 7 13
H 13
tof 2 3
tof 13 3
tof 13 2
T* 2
T 3
tof 3 2
tof 13 3
tof 2 3
tof 13 2
P 13
P* 2
tof 2 13
tof 13 2
H 3
tof 6 1
tof 3 1
tof 3 6
T* 6
T 1
tof 1 6
tof 3 1
tof 6 1
tof 3 6
P 3
P* 6
tof 6 3
tof 3 6
H 1
tof 14 h
tof 1 h
tof 1 14
T* 14
T h
tof h 14
tof 1 h
tof 14 h
tof 1 14
P 1
P* 14
tof 14 1
tof 1 14
H h
tof q 11
tof h 11
tof h q
T* q
T 11
tof 11 q
tof h 11
tof q 11
tof h q
P h
P* q
tof q h
tof h q
H 11
tof s 4
tof 11 4
tof 11 s
T* s
T 4
tof 4 s
tof 11 4
tof s 4
tof 11 s
P 11
P* s
tof s 11
tof 11 s
H 4
tof 8 e
tof 4 e
tof 4 8
T* 8
T e
tof e 8
tof 4 e
tof 8 e
tof 4 8
P 4
P* 8
tof 8 4
tof 4 8
H e
tof 10 k
tof e k
tof e 10
T* 10
T k
tof k 10
tof e k
tof 10 k
tof e 10
P e
P* 10
tof 10 e
tof e 10
H k
tof 9 d
tof k d
tof k 9
T* 9
T d
tof d 9
tof k d
tof 9 d
tof k 9
P k
P* 9
tof 9 k
tof k 9
H d
H t
tof d n
tof d g
tof t g
T d
T t
T* n
T* g
tof n t
tof t d
tof n g
tof d g
tof t n
T* d
T t
tof n d
tof n t
tof d t
H n
H t
tof t e
tof t k
tof 9 t
tof k 9
T* t
T* k
T 9
T* e
tof 9 t
tof e k
tof e 9
tof k 9
tof t k
tof k e
P e
H 9
tof 9 4
tof 9 k
tof 10 9
tof k 10
T* 9
T* k
T 10
T* 4
tof 10 9
tof 4 k
tof 4 10
tof k 10
tof 9 k
tof k 4
P 4
H 10
tof 10 11
tof 10 k
tof 8 10
tof k 8
T* 10
T* k
T 8
T* 11
tof 8 10
tof 11 k
tof 11 8
tof k 8
tof 10 k
tof k 11
P 11
H 8
tof 8 h
tof 8 k
tof s 8
tof k s
T* 8
T* k
T s
T* h
tof s 8
tof h k
tof h s
tof k s
tof 8 k
tof k h
P h
H s
tof s 1
tof s k
tof q s
tof k q
T* s
T* k
T q
T* 1
tof q s
tof 1 k
tof 1 q
tof k q
tof s k
tof k 1
P 1
H q
tof q 3
tof q k
tof 14 q
tof k 14
T* q
T* k
T 14
T* 3
tof 14 q
tof 3 k
tof 3 14
tof k 14
tof q k
tof k 3
P 3
H 14
tof 14 13
tof 14 k
tof 6 14
tof k 6
T* 14
T* k
T 6
T* 13
tof 6 14
tof 13 k
tof 13 6
tof k 6
tof 14 k
tof k 13
P 13
H 6
tof 6 7
tof 6 k
tof 2 6
tof k 2
T* 6
T* k
T 2
T* 7
tof 2 6
tof 7 k
tof 7 2
tof k 2
tof 6 k
tof k 7
P 7
H 2
tof 2 k
P 2
P* k
tof k 2
tof 2 k
H k
tof l b
tof a 12
tof k 12
tof k b
tof 5 j
tof f p
tof k p
tof k j
tof f 5
tof a k
tof g 5
tof l f
tof a l
tof g a
tof a g
tof l a
tof f a
tof k g
tof 5 a
tof j a
tof p 5
tof b f
tof 12 l
T* g
T* a
T j
T p
T b
T* k
T 12
T l
T f
T 5
tof j a
tof k g
tof b p
tof 12 p
tof l p
tof f p
tof k j
tof k p
tof 12 b
tof l b
tof f b
tof f k
tof l 12
tof f l
tof f 5
tof 12 l
tof l 5
tof l f
tof k l
tof b k
tof l 5
tof l f
tof p 5
tof j 5
tof a j
tof g j
tof j 5
tof g f
tof g b
tof g 12
tof f 12
tof f b
tof f a
tof f p
tof f k
tof f l
P* a
P* p
P* b
P* k
P* 12
P* l
tof f a
tof f p
tof f b
tof f k
tof f 12
tof f l
H 0
H f
tof f j
P f
P* j
tof j f
tof f j
H j
tof 7 2
tof j 2
tof j 7
T* 7
T 2
tof 2 7
tof j 2
tof 7 2
tof j 7
P j
P* 7
tof 7 j
tof j 7
H 2
tof 13 6
tof 2 6
tof 2 13
T* 13
T 6
tof 6 13
tof 2 6
tof 13 6
tof 2 13
P 2
P* 13
tof 13 2
tof 2 13
H 6
tof 3 14
tof 6 14
tof 6 3
T* 3
T 14
tof 14 3
tof 6 14
tof 3 14
tof 6 3
P 6
P* 3
tof 3 6
tof 6 3
H 14
tof 1 q
tof 14 q
tof 14 1
T* 1
T q
tof q 1
tof 14 q
tof 1 q
tof 14 1
P 14
P* 1
tof 1 14
tof 14 1
H q
tof h s
tof q s
tof q h
T* h
T s
tof s h
tof q s
tof h s
tof q h
P q
P* h
tof h q
tof q h
H s
tof 11 8
tof s 8
tof s 11
T* 11
T 8
tof 8 11
tof s 8
tof 11 8
tof s 11
P s
P* 11
tof 11 s
tof s 11
H 8
tof 4 10
tof 8 10
tof 8 4
T* 4
T 10
tof 10 4
tof 8 10
tof 4 10
tof 8 4
P 8
P* 4
tof 4 8
tof 8 4
H 10
tof e 9
tof 10 9
tof 10 e
T* e
T 9
tof 9 e
tof 10 9
tof e 9
tof 10 e
P 10
P* e
tof e 10
tof 10 e
H 9
H t
tof 9 d
tof 9 5
tof t 5
T 9
T t
T* d
T* 5
tof d t
tof t 9
tof d 5
tof 9 5
tof t d
T* 9
T t
tof d 9
tof d t
tof 9 t
H d
H t
tof t 8
tof t 10
tof e t
tof 10 e
T* t
T* 10
T e
T* 8
tof e t
tof 8 10
tof 8 e
tof 10 e
tof t 10
tof 10 8
P 8
H e
tof e s
tof e 10
tof 4 e
tof 10 4
T* e
T* 10
T 4
T* s
tof 4 e
tof s 10
tof s 4
tof 10 4
tof e 10
tof 10 s
P s
H 4
tof 4 q
tof 4 10
tof 11 4
tof 10 11
T* 4
T* 10
T 11
T* q
tof 11 4
tof q 10
tof q 11
tof 10 11
tof 4 10
tof 10 q
P q
H 11
tof 11 14
tof 11 10
tof h 11
tof 10 h
T* 11
T* 10
T h
T* 14
tof h 11
tof 14 10
tof 14 h
tof 10 h
tof 11 10
tof 10 14
P 14
H h
tof h 6
tof h 10
tof 1 h
tof 10 1
T* h
T* 10
T 1
T* 6
tof 1 h
tof 6 10
tof 6 1
tof 10 1
tof h 10
tof 10 6
P 6
H 1
tof 1 2
tof 1 10
tof 3 1
tof 10 3
T* 1
T* 10
T 3
T* 2
tof 3 1
tof 2 10
tof 2 3
tof 10 3
tof 1 10
tof 10 2
P 2
H 3
tof 3 j
tof 3 10
tof 13 3
tof 10 13
T* 3
T* 10
T 13
T* j
tof 13 3
tof j 10
tof j 13
tof 10 13
tof 3 10
tof 10 j
P j
H 13
tof 13 f
tof 13 10
tof 7 13
tof 10 7
T* 13
T* 10
T 7
T* f
tof 7 13
tof f 10
tof f 7
tof 10 7
tof 13 10
tof 10 f
P f
H 7
tof 7 10
P 7
P* 10
tof 10 7
tof 7 10
H 10
tof b k
tof 12 l
tof 10 l
tof 10 k
tof 10 p
tof a p
tof 12 10
tof 12 b
tof 12 5
tof b 5
tof a 5
tof 10 12
tof p a
tof k 5
tof l b
P* 12
T* 12
T* a
T p
T k
T* 10
T l
T b
T 5
tof p a
tof 10 12
tof l k
tof b k
tof 10 p
tof 10 k
tof l 10
tof b l
tof l b
tof 10 b
tof b 5
tof k 5
tof p 5
tof a p
tof 12 p
tof p 5
tof 12 b
tof 12 l
tof b l
tof b a
tof b k
tof b 10
P* a
P* k
P* l
P* 10
tof b a
tof b k
tof b l
tof b 10
H g
H b
tof b p
P b
P* p
tof p b
tof b p
H p
tof f 7
tof p 7
tof p f
T* f
T 7
tof 7 f
tof p 7
tof f 7
tof p f
P p
P* f
tof f p
tof p f
H 7
tof j 13
tof 7 13
tof 7 j
T* j
T 13
tof 13 j
tof 7 13
tof j 13
tof 7 j
P 7
P* j
tof j 7
tof 7 j
H 13
tof 2 3
tof 13 3
tof 13 2
T* 2
T 3
tof 3 2
tof 13 3
tof 2 3
tof 13 2
P 13
P* 2
tof 2 13
tof 13 2
H 3
tof 6 1
tof 3 1
tof 3 6
T* 6
T 1
tof 1 6
tof 3 1
tof 6 1
tof 3 6
P 3
P* 6
tof 6 3
tof 3 6
H 1
tof 14 h
tof 1 h
tof 1 14
T* 14
T h
tof h 14
tof 1 h
tof 14 h
tof 1 14
P 1
P* 14
tof 14 1
tof 1 14
H h
tof q 11
tof h 11
tof h q
T* q
T 11
tof 11 q
tof h 11
tof q 11
tof h q
P h
P* q
tof q h
tof h q
H 11
tof s 4
tof 11 4
tof 11 s
T* s
T 4
tof 4 s
tof 11 4
tof s 4
tof 11 s
P 11
P* s
tof s 11
tof 11 s
H 4
tof 8 e
tof 4 e
tof 4 8
T* 8
T e
tof e 8
tof 4 e
tof 8 e
tof 4 8
P 4
P* 8
tof 8 4
tof 4 8
H e
H t
tof e 9
tof 5 e
tof e 5
tof t e
T 5
T t
T* 9
T* e
tof 9 t
tof t 5
tof 9 e
tof 5 e
tof t 9
T* 5
T t
tof 9 5
tof 9 t
tof 5 t
H 9
H t
tof t 11
tof t 4
tof 8 t
tof 4 8
T* t
T* 4
T 8
T* 11
tof 8 t
tof 11 4
tof 11 8
tof 4 8
tof t 4
tof 4 11
P 11
H 8
tof 8 h
tof 8 4
tof s 8
tof 4 s
T* 8
T* 4
T s
T* h
tof s 8
tof h 4
tof h s
tof 4 s
tof 8 4
tof 4 h
P h
H s
tof s 1
tof s 4
tof q s
tof 4 q
T* s
T* 4
T q
T* 1
tof q s
tof 1 4
tof 1 q
tof 4 q
tof s 4
tof 4 1
P 1
H q
tof q 3
tof q 4
tof 14 q
tof 4 14
T* q
T* 4
T 14
T* 3
tof 14 q
tof 3 4
tof 3 14
tof 4 14
tof q 4
tof 4 3
P 3
H 14
tof 14 13
tof 14 4
tof 6 14
tof 4 6
T* 14
T* 4
T 6
T* 13
tof 6 14
tof 13 4
tof 13 6
tof 4 6
tof 14 4
tof 4 13
P 13
H 6
tof 6 7
tof 6 4
tof 2 6
tof 4 2
T* 6
T* 4
T 2
T* 7
tof 2 6
tof 7 4
tof 7 2
tof 4 2
tof 6 4
tof 4 7
P 7
H 2
tof 2 p
tof 2 4
tof j 2
tof 4 j
T* 2
T* 4
T j
T* p
tof j 2
tof p 4
tof p j
tof 4 j
tof 2 4
tof 4 p
P p
H j
tof j b
tof j 4
tof f j
tof 4 f
T* j
T* 4
T f
T* b
tof f j
tof b 4
tof b f
tof 4 f
tof j 4
tof 4 b
P b
H f
tof f 4
P f
P* 4
tof 4 f
tof f 4
H 4
tof a k
tof l 10
tof 4 10
tof 4 k
tof e l
tof l e
tof 4 e
tof a l
tof k a
tof 10 l
P 4
T 4
T* a
T k
T* e
T 10
T l
tof k a
tof e 4
tof e k
tof 10 e
tof 10 l
tof k l
tof e 10
tof a k
tof 4 k
tof k l
tof 4 10
tof 10 a
tof 10 e
P* a
P* e
tof 10 a
tof 10 e
H 12
H 10
tof 10 k
P 10
P* k
tof k 10
tof 10 k
H k
tof b f
tof k f
tof k b
T* b
T f
tof f b
tof k f
tof b f
tof k b
P k
P* b
tof b k
tof k b
H f
tof p j
tof f j
tof f p
T* p
T j
tof j p
tof f j
tof p j
tof f p
P f
P* p
tof p f
tof f p
H j
tof 7 2
tof j 2
tof j 7
T* 7
T 2
tof 2 7
tof j 2
tof 7 2
tof j 7
P j
P* 7
tof 7 j
tof j 7
H 2
tof 13 6
tof 2 6
tof 2 13
T* 13
T 6
tof 6 13
tof 2 6
tof 13 6
tof 2 13
P 2
P* 13
tof 13 2
tof 2 13
H 6
tof 3 14
tof 6 14
tof 6 3
T* 3
T 14
tof 14 3
tof 6 14
tof 3 14
tof 6 3
P 6
P* 3
tof 3 6
tof 6 3
H 14
tof 1 q
tof 14 q
tof 14 1
T* 1
T q
tof q 1
tof 14 q
tof 1 q
tof 14 1
P 14
P* 1
tof 1 14
tof 14 1
H q
tof h s
tof q s
tof q h
T* h
T s
tof s h
tof q s
tof h s
tof q h
P q
P* h
tof h q
tof q h
H s
tof 11 8
tof s 8
tof s 11
T* 11
T 8
tof 8 11
tof s 8
tof 11 8
tof s 11
P s
P* 11
tof 11 s
tof s 11
H 8
H t
tof 8 5
tof l 8
tof 8 l
tof t 8
T l
T t
T* 5
T* 8
tof 5 t
tof t l
tof 5 8
tof l 8
tof t 5
T* l
T t
tof 5 l
tof 5 t
tof l t
H 5
H t
tof t q
tof t s
tof 11 t
tof s 11
T* t
T* s
T 11
T* q
tof 11 t
tof q s
tof q 11
tof s 11
tof t s
tof s q
P q
H 11
tof 11 14
tof 11 s
tof h 11
tof s h
T* 11
T* s
T h
T* 14
tof h 11
tof 14 s
tof 14 h
tof s h
tof 11 s
tof s 14
P 14
H h
tof h 6
tof h s
tof 1 h
tof s 1
T* h
T* s
T 1
T* 6
tof 1 h
tof 6 s
tof 6 1
tof s 1
tof h s
tof s 6
P 6
H 1
tof 1 2
tof 1 s
tof 3 1
tof s 3
T* 1
T* s
T 3
T* 2
tof 3 1
tof 2 s
tof 2 3
tof s 3
tof 1 s
tof s 2
P 2
H 3
tof 3 j
tof 3 s
tof 13 3
tof s 13
T* 3
T* s
T 13
T* j
tof 13 3
tof j s
tof j 13
tof s 13
tof 3 s
tof s j
P j
H 13
tof 13 f
tof 13 s
tof 7 13
tof s 7
T* 13
T* s
T 7
T* f
tof 7 13
tof f s
tof f 7
tof s 7
tof 13 s
tof s f
P f
H 7
tof 7 k
tof 7 s
tof p 7
tof s p
T* 7
T* s
T p
T* k
tof p 7
tof k s
tof k p
tof s p
tof 7 s
tof s k
P k
H p
tof p 10
tof p s
tof b p
tof s b
T* p
T* s
T b
T* 10
tof b p
tof 10 s
tof 10 b
tof s b
tof p s
tof s 10
P 10
H b
tof b s
P b
P* s
tof s b
tof b s
H s
tof s e
tof a e
tof 8 e
tof a s
tof 8 a
tof a 8
tof s 8
tof e a
T 8
T* a
T* s
T e
tof s 8
tof s a
tof s e
tof 8 a
tof a e
tof 8 s
H 4
H s
tof s a
P s
P* a
tof a s
tof s a
H a
tof 10 b
tof a b
tof a 10
T* 10
T b
tof b 10
tof a b
tof 10 b
tof a 10
P a
P* 10
tof 10 a
tof a 10
H b
tof k p
tof b p
tof b k
T* k
T p
tof p k
tof b p
tof k p
tof b k
P b
P* k
tof k b
tof b k
H p
tof f 7
tof p 7
tof p f
T* f
T 7
tof 7 f
tof p 7
tof f 7
tof p f
P p
P* f
tof f p
tof p f
H 7
tof j 13
tof 7 13
tof 7 j
T* j
T 13
tof 13 j
tof 7 13
tof j 13
tof 7 j
P 7
P* j
tof j 7
tof 7 j
H 13
tof 2 3
tof 13 3
tof 13 2
T* 2
T 3
tof 3 2
tof 13 3
tof 2 3
tof 13 2
P 13
P* 2
tof 2 13
tof 13 2
H 3
tof 6 1
tof 3 1
tof 3 6
T* 6
T 1
tof 1 6
tof 3 1
tof 6 1
tof 3 6
P 3
P* 6
tof 6 3
tof 3 6
H 1
tof 14 h
tof 1 h
tof 1 14
T* 14
T h
tof h 14
tof 1 h
tof 14 h
tof 1 14
P 1
P* 14
tof 14 1
tof 1 14
H h
tof q 11
tof h 11
tof h q
T* q
T 11
tof 11 q
tof h 11
tof q 11
tof h q
P h
P* q
tof q h
tof h q
H 11
H t
tof 11 l
tof e 11
tof 11 e
tof t 11
T e
T t
T* l
T* 11
tof l t
tof t e
tof l 11
tof e 11
tof t l
T* e
T t
tof l e
tof l t
tof e t
H l
H t
tof t 1
tof t h
tof q t
tof h q
T* t
T* h
T q
T* 1
tof q t
tof 1 h
tof 1 q
tof h q
tof t h
tof h 1
P 1
H q
tof q 3
tof q h
tof 14 q
tof h 14
T* q
T* h
T 14
T* 3
tof 14 q
tof 3 h
tof 3 14
tof h 14
tof q h
tof h 3
P 3
H 14
tof 14 13
tof 14 h
tof 6 14
tof h 6
T* 14
T* h
T 6
T* 13
tof 6 14
tof 13 h
tof 13 6
tof h 6
tof 14 h
tof h 13
P 13
H 6
tof 6 7
tof 6 h
tof 2 6
tof h 2
T* 6
T* h
T 2
T* 7
tof 2 6
tof 7 h
tof 7 2
tof h 2
tof 6 h
tof h 7
P 7
H 2
tof 2 p
tof 2 h
tof j 2
tof h j
T* 2
T* h
T j
T* p
tof j 2
tof p h
tof p j
tof h j
tof 2 h
tof h p
P p
H j
tof j b
tof j h
tof f j
tof h f
T* j
T* h
T f
T* b
tof f j
tof b h
tof b f
tof h f
tof j h
tof h b
P b
H f
tof f a
tof f h
tof k f
tof h k
T* f
T* h
T k
T* a
tof k f
tof a h
tof a k
tof h k
tof f h
tof h a
P a
H k
tof 10 h
tof k s
tof k h
tof 10 k
tof h 10
T* k
T* 10
T h
T* s
tof h k
tof s 10
tof s h
tof 10 h
tof k 10
tof 10 s
P s
H h
tof h r
tof h m
tof 11 h
tof h 11
tof m h
tof r h
P 11
T 11
T* r
T* m
T h
tof r 11
tof m r
tof m h
tof 11 r
tof 11 h
tof r m
tof m 10
P* 10
tof m 10
H 8
H m
tof s 10
tof m 10
tof m s
T* s
T 10
tof 10 s
tof m 10
tof s 10
tof m s
P m
P* s
tof s m
tof m s
H 10
tof a k
tof 10 k
tof 10 a
T* a
T k
tof k a
tof 10 k
tof a k
tof 10 a
P 10
P* a
tof a 10
tof 10 a
H k
tof b f
tof k f
tof k b
T* b
T f
tof f b
tof k f
tof b f
tof k b
P k
P* b
tof b k
tof k b
H f
tof p j
tof f j
tof f p
T* p
T j
tof j p
tof f j
tof p j
tof f p
P f
P* p
tof p f
tof f p
H j
tof 7 2
tof j 2
tof j 7
T* 7
T 2
tof 2 7
tof j 2
tof 7 2
tof j 7
P j
P* 7
tof 7 j
tof j 7
H 2
tof 13 6
tof 2 6
tof 2 13
T* 13
T 6
tof 6 13
tof 2 6
tof 13 6
tof 2 13
P 2
P* 13
tof 13 2
tof 2 13
H 6
tof 3 14
tof 6 14
tof 6 3
T* 3
T 14
tof 14 3
tof 6 14
tof 3 14
tof 6 3
P 6
P* 3
tof 3 6
tof 6 3
H 14
tof 1 q
tof 14 q
tof 14 1
T* 1
T q
tof q 1
tof 14 q
tof 1 q
tof 14 1
P 14
P* 1
tof 1 14
tof 14 1
H q
H t
tof q e
tof q h
tof t h
T q
T t
T* e
T* h
tof e t
tof t q
tof e h
tof q h
tof t e
T* q
T t
tof e q
tof e t
tof q t
H e
H t
tof t 6
tof t 14
tof 1 t
tof 14 1
T* t
T* 14
T 1
T* 6
tof 1 t
tof 6 14
tof 6 1
tof 14 1
tof t 14
tof 14 6
P 6
H 1
tof 1 2
tof 1 14
tof 3 1
tof 14 3
T* 1
T* 14
T 3
T* 2
tof 3 1
tof 2 14
tof 2 3
tof 14 3
tof 1 14
tof 14 2
P 2
H 3
tof 3 j
tof 3 14
tof 13 3
tof 14 13
T* 3
T* 14
T 13
T* j
tof 13 3
tof j 14
tof j 13
tof 14 13
tof 3 14
tof 14 j
P j
H 13
tof 13 f
tof 13 14
tof 7 13
tof 14 7
T* 13
T* 14
T 7
T* f
tof 7 13
tof f 14
tof f 7
tof 14 7
tof 13 14
tof 14 f
P f
H 7
tof 7 k
tof 7 14
tof p 7
tof 14 p
T* 7
T* 14
T p
T* k
tof p 7
tof k 14
tof k p
tof 14 p
tof 7 14
tof 14 k
P k
H p
tof p 10
tof p 14
tof b p
tof 14 b
T* p
T* 14
T b
T* 10
tof b p
tof 10 14
tof 10 b
tof 14 b
tof p 14
tof 14 10
P 10
H b
tof a 14
tof b m
tof b 14
tof a b
tof 14 a
T* b
T* a
T 14
T* m
tof 14 b
tof m a
tof m 14
tof a 14
tof b a
tof a m
P m
H 14
tof 14 11
tof 14 r
tof 14 a
tof s a
tof h 14
tof s h
tof 14 h
tof a s
tof r 14
tof 11 14
T* h
T* s
T a
T* 11
T* r
T 14
tof a h
tof 11 s
tof s h
tof 11 a
tof r 11
tof r 14
tof 11 r
tof s a
tof h s
tof h 14
tof h r
tof h 11
P r
H a
tof a h
tof a 11
tof a 14
tof 11 14
tof h 14
T a
T* h
T* 11
T 14
tof h a
tof 11 h
tof 11 14
tof a h
tof a 14
tof h 11
H 11
H r
tof m s
tof r s
tof r m
T* m
T s
tof s m
tof r s
tof m s
tof r m
P r
P* m
tof m r
tof r m
H s
tof 10 b
tof s b
tof s 10
T* 10
T b
tof b 10
tof s b
tof 10 b
tof s 10
P s
P* 10
tof 10 s
tof s 10
H b
tof k p
tof b p
tof b k
T* k
T p
tof p k
tof b p
tof k p
tof b k
P b
P* k
tof k b
tof b k
H p
tof f 7
tof p 7
tof p f
T* f
T 7
tof 7 f
tof p 7
tof f 7
tof p f
P p
P* f
tof f p
tof p f
H 7
tof j 13
tof 7 13
tof 7 j
T* j
T 13
tof 13 j
tof 7 13
tof j 13
tof 7 j
P 7
P* j
tof j 7
tof 7 j
H 13
tof 2 3
tof 13 3
tof 13 2
T* 2
T 3
tof 3 2
tof 13 3
tof 2 3
tof 13 2
P 13
P* 2
tof 2 13
tof 13 2
H 3
tof 6 1
tof 3 1
tof 3 6
T* 6
T 1
tof 1 6
tof 3 1
tof 6 1
tof 3 6
P 3
P* 6
tof 6 3
tof 3 6
H 1
H t
tof 1 q
tof 14 1
tof 1 14
tof t 1
T 14
T t
T* q
T* 1
tof q t
tof t 14
tof q 1
tof 14 1
tof t q
T* 14
T t
tof q 14
tof q t
tof 14 t
H q
H t
tof t 13
tof t 3
tof 6 t
tof 3 6
T* t
T* 3
T 6
T* 13
tof 6 t
tof 13 3
tof 13 6
tof 3 6
tof t 3
tof 3 13
P 13
H 6
tof 6 7
tof 6 3
tof 2 6
tof 3 2
T* 6
T* 3
T 2
T* 7
tof 2 6
tof 7 3
tof 7 2
tof 3 2
tof 6 3
tof 3 7
P 7
H 2
tof 2 p
tof 2 3
tof j 2
tof 3 j
T* 2
T* 3
T j
T* p
tof j 2
tof p 3
tof p j
tof 3 j
tof 2 3
tof 3 p
P p
H j
tof j b
tof j 3
tof f j
tof 3 f
T* j
T* 3
T f
T* b
tof f j
tof b 3
tof b f
tof 3 f
tof j 3
tof 3 b
P b
H f
tof f s
tof f 3
tof k f
tof 3 k
T* f
T* 3
T k
T* s
tof k f
tof s 3
tof s k
tof 3 k
tof f 3
tof 3 s
P s
H k
tof 10 3
tof k r
tof k 3
tof 10 k
tof 3 10
T* k
T* 10
T 3
T* r
tof 3 k
tof r 10
tof r 3
tof 10 3
tof k 10
tof 10 r
P r
H 3
tof 3 a
tof 3 h
tof 3 10
tof m 10
tof m 3
tof m 1
tof 3 1
tof 10 m
tof h 1
tof a 1
T* 3
T* m
T 10
T* a
T* h
T 1
tof 10 3
tof a m
tof m 3
tof a 10
tof h a
tof h 1
tof a h
tof m 10
tof 3 m
tof 3 1
tof 3 h
tof 3 a
P h
H 10
tof 10 3
tof 10 a
tof 10 1
tof a 1
tof 3 1
T 10
T* 3
T* a
T 1
tof 3 10
tof a 3
tof a 1
tof 10 3
tof 10 1
tof 3 a
H a
H h
tof r m
tof h m
tof h r
T* r
T m
tof m r
tof h m
tof r m
tof h r
P h
P* r
tof r h
tof h r
H m
tof s k
tof m k
tof m s
T* s
T k
tof k s
tof m k
tof s k
tof m s
P m
P* s
tof s m
tof m s
H k
tof b f
tof k f
tof k b
T* b
T f
tof f b
tof k f
tof b f
tof k b
P k
P* b
tof b k
tof k b
H f
tof p j
tof f j
tof f p
T* p
T j
tof j p
tof f j
tof p j
tof f p
P f
P* p
tof p f
tof f p
H j
tof 7 2
tof j 2
tof j 7
T* 7
T 2
tof 2 7
tof j 2
tof 7 2
tof j 7
P j
P* 7
tof 7 j
tof j 7
H 2
tof 13 6
tof 2 6
tof 2 13
T* 13
T 6
tof 6 13
tof 2 6
tof 13 6
tof 2 13
P 2
P* 13
tof 13 2
tof 2 13
H 6
H t
tof 6 14
tof 1 6
tof 6 1
tof t 6
T 1
T t
T* 14
T* 6
tof 14 t
tof t 1
tof 14 6
tof 1 6
tof t 14
T* 1
T t
tof 14 1
tof 14 t
tof 1 t
H 14
H t
tof t j
tof t 2
tof 13 t
tof 2 13
T* t
T* 2
T 13
T* j
tof 13 t
tof j 2
tof j 13
tof 2 13
tof t 2
tof 2 j
P j
H 13
tof 13 f
tof 13 2
tof 7 13
tof 2 7
T* 13
T* 2
T 7
T* f
tof 7 13
tof f 2
tof f 7
tof 2 7
tof 13 2
tof 2 f
P f
H 7
tof 7 k
tof 7 2
tof p 7
tof 2 p
T* 7
T* 2
T p
T* k
tof p 7
tof k 2
tof k p
tof 2 p
tof 7 2
tof 2 k
P k
H p
tof p m
tof p 2
tof b p
tof 2 b
T* p
T* 2
T b
T* m
tof b p
tof m 2
tof m b
tof 2 b
tof p 2
tof 2 m
P m
H b
tof s 2
tof b h
tof b 2
tof s b
tof 2 s
T* b
T* s
T 2
T* h
tof 2 b
tof h s
tof h 2
tof s 2
tof b s
tof s h
P h
H 2
tof 2 10
tof 2 3
tof 2 s
tof r s
tof r 2
tof r 6
tof 2 6
tof s r
tof 3 6
tof 10 6
T* 2
T* r
T s
T* 10
T* 3
T 6
tof s 2
tof 10 r
tof r 2
tof 10 s
tof 3 10
tof 3 6
tof 10 3
tof r s
tof 2 r
tof 2 6
tof 2 3
tof 2 10
P 3
H s
tof s 2
tof s 10
tof s 6
tof 10 6
tof 2 6
T s
T* 2
T* 10
T 6
tof 2 s
tof 10 2
tof 10 6
tof s 2
tof s 6
tof 2 10
H 10
H 3
tof h r
tof 3 r
tof 3 h
T* h
T r
tof r h
tof 3 r
tof h r
tof 3 h
P 3
P* h
tof h 3
tof 3 h
H r
tof m b
tof r b
tof r m
T* m
T b
tof b m
tof r b
tof m b
tof r m
P r
P* m
tof m r
tof r m
H b
tof k p
tof b p
tof b k
T* k
T p
tof p k
tof b p
tof k p
tof b k
P b
P* k
tof k b
tof b k
H p
tof f 7
tof p 7
tof p f
T* f
T 7
tof 7 f
tof p 7
tof f 7
tof p f
P p
P* f
tof f p
tof p f
H 7
tof j 13
tof 7 13
tof 7 j
T* j
T 13
tof 13 j
tof 7 13
tof j 13
tof 7 j
P 7
P* j
tof j 7
tof 7 j
H 13
H t
tof 13 1
tof 6 13
tof 13 6
tof t 13
T 6
T t
T* 1
T* 13
tof 1 t
tof t 6
tof 1 13
tof 6 13
tof t 1
T* 6
T t
tof 1 6
tof 1 t
tof 6 t
H 1
H t
tof t p
tof t 7
tof j t
tof 7 j
T* t
T* 7
T j
T* p
tof j t
tof p 7
tof p j
tof 7 j
tof t 7
tof 7 p
P p
H j
tof j b
tof j 7
tof f j
tof 7 f
T* j
T* 7
T f
T* b
tof f j
tof b 7
tof b f
tof 7 f
tof j 7
tof 7 b
P b
H f
tof f r
tof f 7
tof k f
tof 7 k
T* f
T* 7
T k
T* r
tof k f
tof r 7
tof r k
tof 7 k
tof f 7
tof 7 r
P r
H k
tof m 7
tof k 3
tof k 7
tof m k
tof 7 m
T* k
T* m
T 7
T* 3
tof 7 k
tof 3 m
tof 3 7
tof m 7
tof k m
tof m 3
P 3
H 7
tof 7 s
tof 7 2
tof 7 m
tof h m
tof 13 7
tof 13 h
tof h 13
tof 7 h
tof m 13
tof 2 7
tof s 7
T* h
T* 13
T m
T* s
T* 2
T 7
tof m h
tof s 13
tof 13 h
tof s m
tof 2 s
tof 2 7
tof s 2
tof 13 m
tof h 13
tof h 7
tof h 2
tof h s
P 2
H m
tof m h
tof m s
tof m 7
tof s 7
tof h 7
T m
T* h
T* s
T 7
tof h m
tof s h
tof s 7
tof m h
tof m 7
tof h s
H s
H 2
tof 3 13
tof 2 13
tof 2 3
T* 3
T 13
tof 13 3
tof 2 13
tof 3 13
tof 2 3
P 2
P* 3
tof 3 2
tof 2 3
H 13
tof r k
tof 13 k
tof 13 r
T* r
T k
tof k r
tof 13 k
tof r k
tof 13 r
P 13
P* r
tof r 13
tof 13 r
H k
tof b f
tof k f
tof k b
T* b
T f
tof f b
tof k f
tof b f
tof k b
P k
P* b
tof b k
tof k b
H f
tof p j
tof f j
tof f p
T* p
T j
tof j p
tof f j
tof p j
tof f p
P f
P* p
tof p f
tof f p
H j
H t
tof j 6
tof 7 j
tof j 7
tof t j
T 7
T t
T* 6
T* j
tof 6 t
tof t 7
tof 6 j
tof 7 j
tof t 6
T* 7
T t
tof 6 7
tof 6 t
tof 7 t
H 6
H t
tof t k
tof t f
tof p t
tof f p
T* t
T* f
T p
T* k
tof p t
tof k f
tof k p
tof f p
tof t f
tof f k
P k
H p
tof p 13
tof p f
tof b p
tof f b
T* p
T* f
T b
T* 13
tof b p
tof 13 f
tof 13 b
tof f b
tof p f
tof f 13
P 13
H b
tof r f
tof b 2
tof b f
tof r b
tof f r
T* b
T* r
T f
T* 2
tof f b
tof 2 r
tof 2 f
tof r f
tof b r
tof r 2
P 2
H f
tof f m
tof f h
tof f r
tof 3 r
tof j f
tof 3 j
tof f j
tof r 3
tof h f
tof m f
T* j
T* 3
T r
T* m
T* h
T f
tof r j
tof m 3
tof 3 j
tof m r
tof h m
tof h f
tof m h
tof 3 r
tof j 3
tof j f
tof j h
tof j m
P h
H r
tof r j
tof r m
tof f r
tof r f
tof m r
tof j r
T f
T* j
T* m
T r
tof j f
tof m j
tof m r
tof f j
tof f r
tof j m
H m
H h
tof 2 3
tof h 3
tof h 2
T* 2
T 3
tof 3 2
tof h 3
tof 2 3
tof h 2
P h
P* 2
tof 2 h
tof h 2
H 3
tof 13 b
tof 3 b
tof 3 13
T* 13
T b
tof b 13
tof 3 b
tof 13 b
tof 3 13
P 3
P* 13
tof 13 3
tof 3 13
H b
tof k p
tof b p
tof b k
T* k
T p
tof p k
tof b p
tof k p
tof b k
P b
P* k
tof k b
tof b k
H p
H t
tof p 7
tof r p
tof p r
tof t p
T r
T t
T* 7
T* p
tof 7 t
tof t r
tof 7 p
tof r p
tof t 7
T* r
T t
tof 7 r
tof 7 t
tof r t
H 7
H t
tof t 3
tof t b
tof k t
tof b k
T* t
T* b
T k
T* 3
tof k t
tof 3 b
tof 3 k
tof b k
tof t b
tof b 3
P 3
H k
tof 13 b
tof k h
tof k b
tof 13 k
tof b 13
T* k
T* 13
T b
T* h
tof b k
tof h 13
tof h b
tof 13 b
tof k 13
tof 13 h
P h
H b
tof b f
tof b j
tof b 13
tof 2 13
tof 2 b
tof 2 p
tof b p
tof 13 2
tof j p
tof f p
T* b
T* 2
T 13
T* f
T* j
T p
tof 13 b
tof f 2
tof 2 b
tof f 13
tof j f
tof j p
tof f j
tof 2 13
tof b 2
tof b p
tof b j
tof b f
P j
H 13
tof 13 b
tof 13 f
tof 13 p
tof f p
tof b p
T 13
T* b
T* f
T p
tof b 13
tof f b
tof f p
tof 13 b
tof 13 p
tof b f
H f
H j
tof h 2
tof j 2
tof j h
T* h
T 2
tof 2 h
tof j 2
tof h 2
tof j h
P j
P* h
tof h j
tof j h
H 2
tof 3 k
tof 2 k
tof 2 3
T* 3
T k
tof k 3
tof 2 k
tof 3 k
tof 2 3
P 2
P* 3
tof 3 2
tof 2 3
H k
H t
tof k r
tof p k
tof k p
tof t k
T p
T t
T* r
T* k
tof r t
tof t p
tof r k
tof p k
tof t r
T* p
T t
tof r p
tof r t
tof p t
H r
H t
tof 3 2
tof t j
tof t 2
tof 3 t
tof 2 3
T* t
T* 3
T 2
T* j
tof 2 t
tof j 3
tof j 2
tof 3 2
tof t 3
tof 3 j
P j
H 2
tof 2 13
tof 2 b
tof 2 3
tof h 3
tof h 2
tof k h
tof h k
tof 2 h
tof 3 k
tof b h
tof 13 h
T* 2
T* k
T 3
T* 13
T* b
T h
tof 3 2
tof 13 k
tof k 2
tof 13 3
tof b 13
tof b h
tof 13 b
tof k 3
tof 2 k
tof 2 h
tof 2 b
tof 2 13
P b
H 3
tof 3 2
tof 3 13
tof 3 h
tof 13 h
tof 2 h
T 3
T* 2
T* 13
T h
tof 2 3
tof 13 2
tof 13 h
tof 3 2
tof 3 h
tof 2 13
H 13
H b
tof j k
tof b k
tof b j
T* j
T k
tof k j
tof b k
tof j k
tof b j
P b
P* j
tof j b
tof b j
H k
H t
tof k t
tof h t
tof k p
tof h k
tof k h
tof t k
T h
T k
T* p
T* t
tof p k
tof k h
tof p t
tof h t
tof k p
T* h
T k
tof p h
tof p k
tof h k
H p
H k
tof k 3
tof k 2
tof k b
tof j b
tof t k
tof j t
tof k t
tof b j
tof 2 k
tof 3 k
T* t
T* j
T b
T* 3
T* 2
T k
tof b t
tof 3 j
tof j t
tof 3 b
tof 2 3
tof 2 k
tof 3 2
tof j b
tof t j
tof t k
tof t 2
tof t 3
P 2
H b
tof b t
tof b 3
tof b k
tof 3 k
tof t k
T b
T* t
T* 3
T k
tof t b
tof 3 t
tof 3 k
tof b t
tof b k
tof t 3
H 3
H 2
H j
tof 2 j
tof k j
tof 2 h
tof k 2
tof 2 k
tof j 2
T k
T 2
T* h
T* j
tof h 2
tof 2 k
tof h j
tof k j
tof 2 h
T* k
T 2
tof h k
tof h 2
tof k 2
H h
H 2
tof 2 b
tof 2 t
tof j 2
tof 2 j
tof t 2
tof b 2
T j
T* b
T* t
T 2
tof b j
tof t b
tof t 2
tof j b
tof j 2
tof b t
H t
H b
tof b j
tof 2 k
tof k 2
tof j k
tof b 2
T b
T* j
T* 2
T k
tof j b
tof 2 j
tof 2 k
tof b j
tof b k
tof j 2
H 2
tof j b
H b
P b
H b
P b
H b
P b
tof 2
Z 2
Y 2
tof h
Z h
Y h
tof t i
tof i t
tof t i
tof i 14
tof 14 i
tof i 14
tof 14 k
tof k 14
tof 14 k
tof k e
tof e k
tof k e
tof e r
tof r e
tof e r
tof r j
tof j r
tof r j
tof j q
tof q j
tof j q
tof q o
tof o q
tof q o
tof o d
tof d o
tof o d
tof d p
tof p d
tof d p
tof p n
tof n p
tof p n
tof n 9
tof 9 n
tof n 9
tof 9 8
tof 8 9
tof 9 8
tof 8 11
tof 11 8
tof 8 11
tof 11 12
tof 12 11
tof 11 12
tof 12 g
tof g 12
tof 12 g
tof g 6
tof 6 g
tof g 6
tof 6 10
tof 10 6
tof 6 10
tof 10 4
tof 4 10
tof 10 4
tof 4 m
tof m 4
tof 4 m
tof m 5
tof 5 m
tof m 5
tof 5 s
tof s 5
tof 5 s
tof s c
tof c s
tof s c
tof c h
tof h c
tof c h
tof h 1
tof 1 h
tof h 1
tof 1 3
tof 3 1
tof 1 3
tof 3 f
tof f 3
tof 3 f
tof f 7
tof 7 f
tof f 7
tof 7 a
tof a 7
tof 7 a
tof a b
tof b a
tof a b
tof b 2
tof 2 b
tof b 2
tof 2 13
tof 13 2
tof 2 13
tof 13 0
tof 0 13
tof 13 0
END
