import qecc as q
import gmpy2

def com_bsv(P, Q):
    return gmpy2.popcount((P[0] & Q[1]) ^ (P[1] & Q[0])) & 1

cdef (int, char) MULT_TABLE[383]

cdef int pphash(char x, char y):
    return x << 2 ^ y

for (x, y), (ph, op) in q.PauliClass.MULT_TABLE.items():
    MULT_TABLE[pphash(ord(x), ord(y))] = ph, ord(op)

def mul(self, other):
    if not isinstance(other, q.Pauli):
        return NotImplemented
    op_py = b'\0' * len(self.op)
    op1_py = self.op.encode()
    op2_py = other.op.encode()
    cdef int ph = self.ph + other.ph
    cdef char *op = op_py
    cdef char *op1 = op1_py
    cdef char *op2 = op2_py
    cdef int i
    for i in range(len(self.op)):
        ph_, op[i] = MULT_TABLE[pphash(op1[i], op2[i])]
        ph += ph_
    return q.Pauli(op.decode(), ph % 4)
