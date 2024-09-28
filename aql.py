import numpy as np

alphabet = (''.join(chr(ord('A') + i) for i in range(26)) +
            ''.join(chr(ord('a') + i) for i in range(26)))
sqrt1_2 = 0.5 ** 0.5

class PureState:
    def __init__(self, mat, qubits = None):
        self.mat = np.array(mat)
        self.mat = self.mat.reshape((2, ) * (self.mat.size.bit_length() - 1))
        if qubits is None:
            self.qubits = [Qubit(container = self) for n in self.mat.shape]
        else:
            assert len(qubits) == len(self.mat.shape)
            self.qubits = qubits

    def __iter__(self):
        return iter(self.qubits)

    def combine_with(self, container):
        if container is not self:
            mat = container.get_state_vector(flatten = False)
            self.mat = np.multiply.outer(self.mat, mat)
            self.qubits.extend(container)
            for q in container:
                q.container = self
            # TODO: destroy the other container?

    def get_indices(self, qubits):
        res = []
        for x in qubits:
            if isinstance(x, Qubit):
                self.combine_with(x.container) # Does nothing if x already in self
                res.append(self.qubits.index(x))
            else: res.append(x)
        return res

    def apply_gate(self, gate, qubits):
        qubits = self.get_indices(qubits)
        self.mat = gate.apply_to_pure_state(self.mat, qubits)

    def get_density_matrix(self, qubits = None, flatten = True):
        if qubits is None: qubits = range(len(self.qubits))
        else: qubits = self.get_indices(qubits)
        n = len(self.qubits)
        m = len(qubits)
        # Example: if n = 5 and qubits = [3, 1], then index_str will be ABCDE, AGCFE -> DBFG
        lhs = ''.join(alphabet[n + qubits.index(i) if i in qubits else i] for i in range(n))
        rhs = ''.join(alphabet[x] for x in qubits) + alphabet[n:n+m]
        index_str = alphabet[:n] + ', ' + lhs + ' -> ' + rhs
        res = np.einsum(index_str, self.mat, self.mat.conj())
        if flatten: return res.reshape((2 ** m, 2 ** m))
        return res 

    def get_state_vector(self, qubits = None, flatten = True):
        if qubits is None: res = self.mat
        else:
            qubits = self.get_indices(qubits)
            assert len(qubits) == len(self.qubits) # TODO
            res = self.mat.transpose(qubits)
        if flatten: return res.flatten()
        return res.copy()

class Qubit:
    named_states = {
            '0': [1, 0],
            '1': [0, 1],
            '+': [sqrt1_2, sqrt1_2],
            '-': [sqrt1_2, -sqrt1_2],
            '+i': [sqrt1_2, sqrt1_2 * 1j],
            '-i': [sqrt1_2, sqrt1_2 * -1j],
    }

    def __init__(self, mat = None, container = None):
        if isinstance(mat, str):
            mat = Qubit.named_states[mat]
        if container is None:
            self.container = PureState(mat, qubits = [self])
        else:
            self.container = container

    def apply_gate(self, gate, qubits = []):
        self.container.apply_gate(gate, [self] + qubits)

    def get_density_matrix(self):
        return self.container.get_density_matrix(qubits = [self])

    def get_state_vector(self):
        return self.container.get_state_vector(qubits = [self])

class Gate:
    def __init__(self, mat):
        self.mat = np.array(mat)
        self.mat = self.mat.reshape((2, ) * (self.mat.size.bit_length() - 1))
        assert len(self.mat.shape) % 2 == 0
        self.fan_in = len(self.mat.shape) // 2

    def apply_to_pure_state(self, state, qubits):
        n = len(state.shape)
        m = self.fan_in
        # Example: if n = 5 and qubits = [3, 1], then index_str will be ABCDE, DBFG -> AGCFE
        lhs = ''.join(alphabet[x] for x in qubits) + alphabet[n:n+m]
        rhs = ''.join(alphabet[n + qubits.index(i) if i in qubits else i] for i in range(n))
        index_str = alphabet[:n] + ', ' + lhs + ' -> ' + rhs
        return np.einsum(index_str, state, self.mat)

X = Gate([[0, 1], [1, 0]])
Y = Gate([[0, -1j], [1j, 0]])
Z = Gate([[1, 0], [0, -1]])
H = Gate([[sqrt1_2, sqrt1_2], [sqrt1_2, -sqrt1_2]])
T = Gate([[1, 0], [0, sqrt1_2 + sqrt1_2 * 1j]])
CX = Gate(np.eye(4)[[0, 1, 3, 2]])

class Circuit:
    pass
