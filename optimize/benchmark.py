import numpy as np
import qecc as q
import sys
sys.path.append('..')
import time
import aql
import optimize.T_optimizer as to
import sys
import re

gate_dict = {'X': aql.X, 'Y': aql.Y, 'Z': aql.Z, 'H': aql.H, 'T': aql.T, 'TOF': aql.CX, 'CX': aql.CX, 'CNOT': aql.CX}
gate_dict['CZ'] = aql.Gate(np.diag([1] * 3 + [-1]))
gate_dict['CCZ'] = aql.Gate(np.diag([1] * 7 + [-1]))
gate_dict['SWAP'] = aql.Gate(np.eye(4)[[0, 2, 1, 3]])
gate_dict['T*'] = aql.Gate([[1, 0], [0, aql.sqrt1_2 - aql.sqrt1_2 * 1j]])
gate_dict['P'] = gate_dict['S'] = aql.Gate([[1, 0], [0, 1j]])
gate_dict['P*'] = gate_dict['S*'] = aql.Gate([[1, 0], [0, -1j]])

def eval_circuit(c, state = None):
    n = len(c['qubits'])
    if state is None:
        EPR_pairs = [aql.PureState([1, 0, 0, 1]) for i in range(n)]
        qubits = [pair.qubits[0] for pair in EPR_pairs]
    else:
        qubits = aql.PureState(state).qubits
    for gate, targets in c['gates']:
        gate = gate.upper()
        if gate == 'Z' and len(targets) == 3: gate = 'CCZ'
        elif gate == 'TOF' and len(targets) == 1: gate = 'X'
        qubits[targets[0]].apply_gate(gate_dict[gate], qubits = [qubits[i] for i in targets[1:]])
    if state is None: return qubits[0].container.get_state_vector(qubits = [pair.qubits[1] for pair in EPR_pairs] + qubits).reshape(1 << n, 1 << n)
    return qubits[0].container.get_state_vector(qubits = qubits)

def test(file):
    with open(file, 'r') as fin: c = to.from_QC_file(fin)
    n = len(c['qubits'])
    print("#Input File:"+str(file)+"\n#Number of Qubits:"+str(n))
    T, C = to.get_T_paulis(c)
    print("#T-count : ", str(c['gates']).count('\'T\'')+str(c['gates']).count('\'T*\'')+sum(1 for gate, targets in c['gates'] if gate == 'Z' and len(targets) == 3)*7)
    print("#CNOT-count : ", str(c['gates']).count('\'CX\'')+str(c['gates']).count('\'CZ\'')+sum(1 for gate, targets in c['gates'] if gate.upper() == 'TOF' and len(targets) == 2)+sum(1 for gate, targets in c['gates'] if gate =='Z' and len(targets)==3)*6)

    '''
    TT, CC = to.remove_duplicates(T, C)
    c_new = {**c, 'gates': to.from_T_paulis(TT, CC)}
    print(len(T), len(TT))
    '''
    inds, dups = to.remove_duplicates(T, C, indices = True)
    c_new = to.remove_T_gates(c, inds, dups)
#    print(c_new['gates'])
#    print("#Original T-count", len(T))
    print("#After Optimization \n#T-count : ", str(c_new['gates']).count('\'T\'')+str(c_new['gates']).count('\'T*\''))
    print("#CNOT-count : ", str(c_new['gates']).count('\'CX\'')+str(c_new['gates']).count('\'CZ\'')+sum(1 for gate, targets in c_new['gates'] if gate.upper() == 'TOF' and len(targets) == 2))
    if n <= 20:
        if n <= 10:
            print('Verifying…')
            state = None
        else:
            print('Partially verifying…')
            state = np.random.randn(1 << n) + np.random.randn(1 << n) * 1j
        U = eval_circuit(c, state = state)
        U_new = eval_circuit(c_new, state = state)
        i = abs(U).argmax()
        lam = U_new.reshape(-1)[i] / U.reshape(-1)[i]
        assert np.allclose(abs(lam), 1)
        assert np.allclose(U_new, U * lam)

def check_depth(file):
    test(file)
    '''
    start = time.time()
    with open(file, 'r') as fin: c = to.from_QC_file(fin)
    n = len(c['qubits'])
    print(file, n)
    T, C = to.get_T_paulis(c)
    TT, CC = to.remove_duplicates(T, C)
    depths = [1] * len(TT)
    for i, P in enumerate(TT):
        for j, (Q, d) in enumerate(zip(TT[:i], depths)):
            if depths[i] < d + 1 and q.com(P, Q): depths[i] = d + 1
    c_new = {**c, 'gates': to.from_T_paulis(TT, CC)}
    print("T-count : ", str(c_new['gates']).count('\'T\'')+str(c_new['gates']).count('\'T*\''))
    print("CNOT-count : ", str(c_new['gates']).count('\'CX\'')+str(c_new['gates']).count('\'CZ\''))
    print(len(T), len(TT), max(depths))
    elapsed = time.time()-start
    print(elapsed)
    '''

file=sys.argv[1]
check_depth(file)
