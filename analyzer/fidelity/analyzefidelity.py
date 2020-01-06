import numpy as np
import matplotlib.pyplot as plt
import copy

from tqdm import trange, tqdm
from qiskit import QuantumRegister, QuantumCircuit, ClassicalRegister
from qiskit import Aer, execute
from qiskit.quantum_info import state_fidelity
from qiskit.providers.aer.noise import NoiseModel
from qiskit.providers.aer.noise.errors import depolarizing_error
from qiskit import transpile


QASM = Aer.get_backend('qasm_simulator')
STATEVEC = Aer.get_backend('statevector_simulator')


class FidelityAnalyzer:
    '''
    This analyzer is for analyzing the fidelity
    performance under the noisy situation.
    Errors are very naive.
    If you put quantum circuit, this module analyze the
    fidelity decreasing automatically.
    '''

    def __init__(self, one_error, two_error, measure_qubit,
                 extime=10, shots=5000, drawing=True):
        if not isinstance(one_error, (float, int, np.ndarray, list)):
            raise ValueError('one error is must be float, int, array or list.')
        else:
            # error rate of u3
            self.one_error = one_error

        if not isinstance(two_error, (float, int, np.ndarray, list)):
            raise ValueError('one error is must be float, int, array or list.')
        else:
            # error rate of cx
            self.two_error = two_error
        # TODO: make 3D plot when one error and two error are array
        self.shots = shots
        self.mes_qbt = measure_qubit
        self.extime = extime
        self.drawing = drawing

    def fidelity_drop(self, qc):
        nqc = transpile(qc, optimization_level=0, basis_gates=['cx', 'u3'])
        # HACK efficient ways?
        if isinstance(self.one_error,
                      (float, int)) and isinstance(self.two_error,
                                                   (float, int)):
            fidelitis = self.fixed_fidelity(nqc)
        # FIXME more seeable
        elif isinstance(self.one_error,
                        (float, int)) and isinstance(self.two_error,
                                                     (np.ndarray, list)):
            fidelities, std = self._u3fix(nqc)
            print('fid!!', fidelities, std)
        elif isinstance(self.two_error,
                        (float, int)) and isinstance(self.one_error,
                                                     (np.ndarray, list)):
            cxerror = depolarizing_error(self.two_error, 2)
            fidelities = self._cxfix(nqc)
        else:
            fidelities = self._nofix(nqc)
        if self.drawing:
            self._draw()
        return fidelities

    def _u3fix(self, qc):
        print(qc.count_ops())
        nst = 2**len(self.mes_qbt)
        bins = [format(i, '0%db' % len(self.mes_qbt))
                for i in range(nst)]

        # ideal result of this circuit
        ideal = execute(qc, backend=QASM, shots=self.shots*10)
        idealcounts = ideal.result().get_counts()
        idealst = np.array([idealcounts.get(i, 0)/(self.shots*10)
                           for i in bins])
        # making noise model with error rate
        u3error = depolarizing_error(self.one_error, 1)
        # start simulations
        mean_fid = []
        std_fid = []
        for cxerror in tqdm(self.two_error):
            mid = []
            # HACK: might be efficient in top layer
            noise_model = NoiseModel()
            noise_model.add_all_qubit_quantum_error(u3error, ['u3'])
            cxerror = depolarizing_error(cxerror, 2)
            noise_model.add_all_qubit_quantum_error(cxerror, ['cx'])
            for t in range(self.extime):
                # execute!
                job = execute(qc, backend=QASM, noise_model=noise_model,
                              shots=self.shots)
                counts = job.result().get_counts()
                stvec = [counts.get(i, 0)/self.shots for i in bins]
                stf = state_fidelity(idealst, stvec)
                mid.append(stf)
            mean_fid.append(np.mean(mid))
            std_fid.append(np.std(mid))
        return mean_fid, std_fid

    def _cxfix(self):
        pass

    def _nofix(self):
        pass

    def _fixed_fidelity(self, qc):
        fidelity = 0
        return fidelity

    def _draw(self, errorbar=True, **kwargs):
        '''
        drawing fidelity dropping
        '''
        pass


if __name__ == '__main__':
    q = QuantumRegister(4)
    c = ClassicalRegister(2)
    qc = QuantumCircuit(q, c)
    qc.x(q[0])
    qc.cx(q[0], q[1])
    qc.measure(q[0], c[0])
    qc.measure(q[1], c[1])
    analyzer = FidelityAnalyzer(0.1, np.arange(0, 1, 0.01), [0, 1])
    result = analyzer.fidelity_drop(qc)
