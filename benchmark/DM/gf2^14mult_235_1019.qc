# (c) Copyright 2004-2010 Dmitri Maslov
# Distributing this file without an express written permission or erasing this note is prohibited.
# Primitive polynomial x^14+x^5+x^3+x+1
.v a0 a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 a11 a12 a13 b0 b1 b2 b3 b4 b5 b6 b7 b8 b9 b10 b11 b12 b13 c0 c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12 c13
.i a0 a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 a11 a12 a13 b0 b1 b2 b3 b4 b5 b6 b7 b8 b9 b10 b11 b12 b13
.o c0 c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12 c13
.c 0 0 0 0 0 0 0 0 0 0 0 0 0 0
BEGIN
H c0
Z a13 b1 c0
H c0
H c0
Z a12 b2 c0
H c0
H c0
Z a11 b3 c0
H c0
H c0
Z a10 b4 c0
H c0
H c0
Z a9 b5 c0
H c0
H c0
Z a8 b6 c0
H c0
H c0
Z a7 b7 c0
H c0
H c0
Z a6 b8 c0
H c0
H c0
Z a5 b9 c0
H c0
H c0
Z a4 b10 c0
H c0
H c0
Z a3 b11 c0
H c0
H c0
Z a2 b12 c0
H c0
H c0
Z a1 b13 c0
H c0
H c1
Z a13 b2 c1
H c1
H c1
Z a12 b3 c1
H c1
H c1
Z a11 b4 c1
H c1
H c1
Z a10 b5 c1
H c1
H c1
Z a9 b6 c1
H c1
H c1
Z a8 b7 c1
H c1
H c1
Z a7 b8 c1
H c1
H c1
Z a6 b9 c1
H c1
H c1
Z a5 b10 c1
H c1
H c1
Z a4 b11 c1
H c1
H c1
Z a3 b12 c1
H c1
H c1
Z a2 b13 c1
H c1
H c2
Z a13 b3 c2
H c2
H c2
Z a12 b4 c2
H c2
H c2
Z a11 b5 c2
H c2
H c2
Z a10 b6 c2
H c2
H c2
Z a9 b7 c2
H c2
H c2
Z a8 b8 c2
H c2
H c2
Z a7 b9 c2
H c2
H c2
Z a6 b10 c2
H c2
H c2
Z a5 b11 c2
H c2
H c2
Z a4 b12 c2
H c2
H c2
Z a3 b13 c2
H c2
H c3
Z a13 b4 c3
H c3
H c3
Z a12 b5 c3
H c3
H c3
Z a11 b6 c3
H c3
H c3
Z a10 b7 c3
H c3
H c3
Z a9 b8 c3
H c3
H c3
Z a8 b9 c3
H c3
H c3
Z a7 b10 c3
H c3
H c3
Z a6 b11 c3
H c3
H c3
Z a5 b12 c3
H c3
H c3
Z a4 b13 c3
H c3
H c4
Z a13 b5 c4
H c4
H c4
Z a12 b6 c4
H c4
H c4
Z a11 b7 c4
H c4
H c4
Z a10 b8 c4
H c4
H c4
Z a9 b9 c4
H c4
H c4
Z a8 b10 c4
H c4
H c4
Z a7 b11 c4
H c4
H c4
Z a6 b12 c4
H c4
H c4
Z a5 b13 c4
H c4
H c5
Z a13 b6 c5
H c5
H c5
Z a12 b7 c5
H c5
H c5
Z a11 b8 c5
H c5
H c5
Z a10 b9 c5
H c5
H c5
Z a9 b10 c5
H c5
H c5
Z a8 b11 c5
H c5
H c5
Z a7 b12 c5
H c5
H c5
Z a6 b13 c5
H c5
H c6
Z a13 b7 c6
H c6
H c6
Z a12 b8 c6
H c6
H c6
Z a11 b9 c6
H c6
H c6
Z a10 b10 c6
H c6
H c6
Z a9 b11 c6
H c6
H c6
Z a8 b12 c6
H c6
H c6
Z a7 b13 c6
H c6
H c7
Z a13 b8 c7
H c7
H c7
Z a12 b9 c7
H c7
H c7
Z a11 b10 c7
H c7
H c7
Z a10 b11 c7
H c7
H c7
Z a9 b12 c7
H c7
H c7
Z a8 b13 c7
H c7
H c8
Z a13 b9 c8
H c8
H c8
Z a12 b10 c8
H c8
H c8
Z a11 b11 c8
H c8
H c8
Z a10 b12 c8
H c8
H c8
Z a9 b13 c8
H c8
H c9
Z a13 b10 c9
H c9
H c9
Z a12 b11 c9
H c9
H c9
Z a11 b12 c9
H c9
H c9
Z a10 b13 c9
H c9
H c10
Z a13 b11 c10
H c10
H c10
Z a12 b12 c10
H c10
H c10
Z a11 b13 c10
H c10
H c11
Z a13 b12 c11
H c11
H c11
Z a12 b13 c11
H c11
H c12
Z a13 b13 c12
H c12
CX c12 c3
CX c12 c1
CX c12 c13
CX c11 c2
CX c11 c0
CX c11 c12
CX c10 c1
CX c10 c13
CX c10 c11
CX c9 c0
CX c9 c12
CX c9 c10
CX c8 c13
CX c8 c11
CX c8 c9
CX c7 c12
CX c7 c10
CX c7 c8
CX c6 c11
CX c6 c9
CX c6 c7
CX c5 c10
CX c5 c8
CX c5 c6
CX c4 c9
CX c4 c7
CX c4 c5
CX c3 c8
CX c3 c6
CX c3 c4
CX c2 c7
CX c2 c5
CX c2 c3
CX c1 c6
CX c1 c4
CX c1 c2
CX c0 c5
CX c0 c3
CX c0 c1
H c13
Z a13 b0 c13
H c13
H c13
Z a12 b1 c13
H c13
H c13
Z a11 b2 c13
H c13
H c13
Z a10 b3 c13
H c13
H c13
Z a9 b4 c13
H c13
H c13
Z a8 b5 c13
H c13
H c13
Z a7 b6 c13
H c13
H c13
Z a6 b7 c13
H c13
H c13
Z a5 b8 c13
H c13
H c13
Z a4 b9 c13
H c13
H c13
Z a3 b10 c13
H c13
H c13
Z a2 b11 c13
H c13
H c13
Z a1 b12 c13
H c13
H c13
Z a0 b13 c13
H c13
H c12
Z a12 b0 c12
H c12
H c12
Z a11 b1 c12
H c12
H c12
Z a10 b2 c12
H c12
H c12
Z a9 b3 c12
H c12
H c12
Z a8 b4 c12
H c12
H c12
Z a7 b5 c12
H c12
H c12
Z a6 b6 c12
H c12
H c12
Z a5 b7 c12
H c12
H c12
Z a4 b8 c12
H c12
H c12
Z a3 b9 c12
H c12
H c12
Z a2 b10 c12
H c12
H c12
Z a1 b11 c12
H c12
H c12
Z a0 b12 c12
H c12
H c11
Z a11 b0 c11
H c11
H c11
Z a10 b1 c11
H c11
H c11
Z a9 b2 c11
H c11
H c11
Z a8 b3 c11
H c11
H c11
Z a7 b4 c11
H c11
H c11
Z a6 b5 c11
H c11
H c11
Z a5 b6 c11
H c11
H c11
Z a4 b7 c11
H c11
H c11
Z a3 b8 c11
H c11
H c11
Z a2 b9 c11
H c11
H c11
Z a1 b10 c11
H c11
H c11
Z a0 b11 c11
H c11
H c10
Z a10 b0 c10
H c10
H c10
Z a9 b1 c10
H c10
H c10
Z a8 b2 c10
H c10
H c10
Z a7 b3 c10
H c10
H c10
Z a6 b4 c10
H c10
H c10
Z a5 b5 c10
H c10
H c10
Z a4 b6 c10
H c10
H c10
Z a3 b7 c10
H c10
H c10
Z a2 b8 c10
H c10
H c10
Z a1 b9 c10
H c10
H c10
Z a0 b10 c10
H c10
H c9
Z a9 b0 c9
H c9
H c9
Z a8 b1 c9
H c9
H c9
Z a7 b2 c9
H c9
H c9
Z a6 b3 c9
H c9
H c9
Z a5 b4 c9
H c9
H c9
Z a4 b5 c9
H c9
H c9
Z a3 b6 c9
H c9
H c9
Z a2 b7 c9
H c9
H c9
Z a1 b8 c9
H c9
H c9
Z a0 b9 c9
H c9
H c8
Z a8 b0 c8
H c8
H c8
Z a7 b1 c8
H c8
H c8
Z a6 b2 c8
H c8
H c8
Z a5 b3 c8
H c8
H c8
Z a4 b4 c8
H c8
H c8
Z a3 b5 c8
H c8
H c8
Z a2 b6 c8
H c8
H c8
Z a1 b7 c8
H c8
H c8
Z a0 b8 c8
H c8
H c7
Z a7 b0 c7
H c7
H c7
Z a6 b1 c7
H c7
H c7
Z a5 b2 c7
H c7
H c7
Z a4 b3 c7
H c7
H c7
Z a3 b4 c7
H c7
H c7
Z a2 b5 c7
H c7
H c7
Z a1 b6 c7
H c7
H c7
Z a0 b7 c7
H c7
H c6
Z a6 b0 c6
H c6
H c6
Z a5 b1 c6
H c6
H c6
Z a4 b2 c6
H c6
H c6
Z a3 b3 c6
H c6
H c6
Z a2 b4 c6
H c6
H c6
Z a1 b5 c6
H c6
H c6
Z a0 b6 c6
H c6
H c5
Z a5 b0 c5
H c5
H c5
Z a4 b1 c5
H c5
H c5
Z a3 b2 c5
H c5
H c5
Z a2 b3 c5
H c5
H c5
Z a1 b4 c5
H c5
H c5
Z a0 b5 c5
H c5
H c4
Z a4 b0 c4
H c4
H c4
Z a3 b1 c4
H c4
H c4
Z a2 b2 c4
H c4
H c4
Z a1 b3 c4
H c4
H c4
Z a0 b4 c4
H c4
H c3
Z a3 b0 c3
H c3
H c3
Z a2 b1 c3
H c3
H c3
Z a1 b2 c3
H c3
H c3
Z a0 b3 c3
H c3
H c2
Z a2 b0 c2
H c2
H c2
Z a1 b1 c2
H c2
H c2
Z a0 b2 c2
H c2
H c1
Z a1 b0 c1
H c1
H c1
Z a0 b1 c1
H c1
H c0
Z a0 b0 c0
H c0
END
