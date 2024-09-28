import qecc as q
from gmpy2 import mpz
from collections import Counter

def com(P, Q):
    return sum(1 for P1, Q1 in zip(P.op, Q.op) if P1 != Q1 and 'I' not in [P1, Q1]) & 1

def bsv(P):
    x = z = 0
    for c in P.op:
        x = (x << 1) | (c in 'XY')
        z = (z << 1) | (c in 'YZ')
    return mpz(x), mpz(z)

q.com = com

import pyximport
pyximport.install(language_level = 3)

from . import fastq

def mul(self, other):
    return fastq.mul(self, other)

q.Pauli.__mul__ = mul
com_bsv = fastq.com_bsv

def pauli_init(self, operator, phase=0):
    # don't validate the input string
    if not isinstance(phase, int):
        raise ValueError("Input phase must be an integer, preferably 0 to 3.")

    #If a phase outside of range(4) is input, we use its remainder, mod 4.
    if not( phase > -1 and phase < 4):
        phase= phase % 4

    self.op = operator
    self.ph = phase

q.Pauli.__init__ = pauli_init

def from_QC_file(fin):
    inputs = None
    outputs = None
    gates = []
    for line in fin:
        words = line.split()
        if words != []:
            gate = words[0]
            targets = words[1:]
            if gate == '.v':
                qubits = targets
                qubit_dict = {q: i for i, q in enumerate(qubits)}
            elif gate in ['#', '.c', 'BEGIN', 'END']: pass
            else:
                targets = tuple(qubit_dict.get(q) for q in targets)
                if gate == '.i': inputs = targets
                elif gate == '.o': outputs = targets
                else: gates.append((gate, targets))
    return {'qubits': qubits, 'inputs': inputs, 'outputs': outputs, 'gates': gates}

def to_QC_file(circuit, fout):
    qubits = circuit['qubits']
    fout.write('.v')
    if qubits is not None: 
        for q in qubits: fout.write(' '+q) 
    fout.write('\n')
    inputs = circuit['inputs']
    fout.write('.i')
    if inputs is not None:
        for q in inputs: fout.write(' '+qubits[q]) 
    fout.write('\n')
    outputs = circuit['outputs']
    fout.write('.o')
    if outputs is not None:
        for q in outputs: fout.write(' '+qubits[q]) 
    fout.write('\n')

    fout.write('\n')
    fout.write('BEGIN\n')
    gates = circuit['gates']
    if gates is not None:
        for g in gates:
            name = g[0]
            name = name.lstrip('C')
            if name == 'X': name = 'tof'
            fout.write(name)
            for indexQubit in g[1]: fout.write(' '+qubits[indexQubit]) 
            fout.write('\n')
    fout.write('END\n')


def get_T_paulis(circuit):
    n = len(circuit['qubits'])
    C = q.eye_c(n)
    Xgens, Zgens = q.elem_gens(n)
    res = []
    for gate, targets in circuit['gates']:
        gate = gate.upper()
        if gate in ['T', 'T*']:
            P = Zgens[targets[0]]
            if gate == 'T*': P = -P
            res.append(C(P))
        elif gate == 'Z' and len(targets) == 3:
            P = -q.eye_p(n)
            for k in 0, 1, 0, 2, 0, 1, 0:
                P *= -Zgens[targets[k]]
                res.append(C(P))
        elif gate == 'H':
            t = targets[0]
            C.xout[t], C.zout[t] = C.zout[t], C.xout[t]
        elif gate in ['TOF', 'CX', 'CNOT']:
            if len(targets) ==2:
                C *= q.cnot(n, *targets)
            elif len(targets) ==1:
                C *=Xgens[targets[0]].as_clifford()
            else:
                raise ValueError('parameters error: %s' % gate)
        elif gate in 'XYZ':
            if gate in 'XY': C *= Xgens[targets[0]].as_clifford()
            if gate in 'YZ': C *= Zgens[targets[0]].as_clifford()
        elif gate in ['P', 'P*', 'S', 'S*']:
            C *= q.phase(n, targets[0])
            if gate.endswith('*'): C *= Zgens[targets[0]].as_clifford()
        else: raise ValueError('Unsupported gate: %s' % gate)
    return res, C

def phase_clifford(P):
    Xgens, Zgens = q.elem_gens(P.nq)
    for ncp, gens in ('YZ', Xgens), ('XY', Zgens):
        for i, Q in enumerate(gens):
            if P.op[i] in ncp:
                gens[i] = Q.mul_phase(1) * P
    return q.Clifford(Xgens, Zgens)

def phase_clifford_mult(P, C):
    for gens in C.xout, C.zout:
        for i, Q in enumerate(gens):
            res = Q * P
            if res.ph & 1: gens[i] = res.mul_phase(1)

def remove_duplicates(T_paulis, C0 = None, indices = False):
    if T_paulis == []: return([], []) if indices else ([], C0)
    res = []
    dups = []
    res_cnt = Counter()
    C = q.eye_c(T_paulis[0].nq)
    for j, Q in enumerate(T_paulis):
        Q = C(Q)
        Qb = bsv(Q)
        if res_cnt[Qb]:
            for i, (_, P, Pb) in enumerate(reversed(res), start = 1):
                if Pb == Qb:
                    res.pop(-i)
                    res_cnt[Qb] -= 1
                    if P.ph == Q.ph:
                        phase_clifford_mult(P, C)
                        dups.append(j)
                    Q = None
                    break
                elif com_bsv(Pb, Qb): break
        if Q is not None:
            res.append((j, Q, Qb))
            res_cnt[Qb] += 1
    if indices: return [i for i, P, Pb in res], dups
    if C0 is not None: C *= C0
    return [P for i, P, Pb in res], C

CCZ_gates = [('T', (0,), 0), ('T', (1,), 2), ('CX', (2, 0), None), ('CX', (1, 2), None), ('T*', (0,), 5), ('T*', (2,), 3), ('CX', (1, 0), None), ('CX', (1, 2), None), ('T', (0,), 4), ('CX', (2, 0), None), ('T*', (0,), 1), ('T', (2,), 6), ('CX', (1, 0), None)]

def remove_T_gates(circuit, inds, dups):
    inds = set(inds)
    dups = set(dups)
    res = []
    j = 0
    for gate, targets in circuit['gates']:
        if gate in ['T', 'T*']:
            if j in inds: res.append((gate, targets))
            elif j in dups: res.append(('P' + gate[1:], targets))
            j += 1
        elif gate == 'Z' and len(targets) == 3:
            for gate, tar_i, jj in CCZ_gates:
                tar = tuple(targets[i] for i in tar_i)
                if gate == 'CX' or j + jj in inds: res.append((gate, tar))
                elif j + jj in dups: res.append(('P' + gate[1:], tar))
            j += 7
        else: res.append((gate, targets))
    return {**circuit, 'gates': res}

def from_clifford(C):
    gate_mapping = {'CNOT': 'CX', 'R_pi4': 'P'}
    for kind in 'X', 'Y', 'Z', 'H', 'CZ', 'SWAP': gate_mapping[kind] = kind
    cd = C.circuit_decomposition()
    return [(gate_mapping[gate.kind], gate.qubits) for gate in reversed(cd)]

# returns a gates list
# TODO: synthesize T in layers
def from_T_paulis(T_paulis, C0 = None):
    if T_paulis == []: return from_clifford(C0)
    res = []
    n = T_paulis[0].nq
    C = q.eye_c(n)
    T_layer = {}
    def finish_layer():
        nonlocal T_layer
        res.extend((gate, (target,)) for target, gate in T_layer.items())
        T_layer = {}
    for P in T_paulis:
        P = C(P)
        gate = 'T*' if P.ph else 'T'
        P.set_phase(0)
        targets = []
        target = None
        for i, p in enumerate(P):
            if p != q.I:
                if p == q.Y:
                    res.append(('P', (i,)))
                    C = q.phase(n, i).inv() * C
                if p != q.Z:
                    if i in T_layer: finish_layer()
                    res.append(('H', (i,)))
                    C = q.hadamard(n, i) * C
                targets.append(i)
                if i not in T_layer: target = i
        if targets != []:
            if target is None:
                finish_layer()
                target = targets[-1]
            for i in targets:
                if i != target:
                    res.append(('CX', (i, target)))
                    C = q.cnot(n, i, target) * C
            T_layer[target] = gate
    finish_layer()
    if C0 is not None: C *= C0
    return res + from_clifford(C)

def optimize(T_paulis, C0 = None):
    if T_paulis == []: return [], C0
    res = []
    C = q.eye_c(T_paulis[0].nq)
    for Q in T_paulis:
        Q = C(Q)
        for i, P in enumerate(reversed(res), start = 1):
            if P == Q or P == -Q:
                res.pop(-i)
                if P == Q: C = phase_clifford(P) * C
                Q = None
                break
            elif q.com(P, Q): break
        if Q is not None: res.append(Q)
    if C0 is not None: C *= C0
    return res, C
