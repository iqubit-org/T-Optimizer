# Original circuit
#   qubits: 9
#   qubits used: 9
#   H: 14
#   cnot: 55
#   X: 8
#   T: 49
#   P: 0
#   Z: 0
#   tdepth (by partitions): 21
#   depth  (by critical paths): 54
#   tdepth (by critical paths): 15
# Optimized circuit
#   qubits: 9
#   qubits used: 9
#   H: 14
#   cnot: 104
#   X: 13
#   T: 35
#   P: 2
#   Z: 0
#   tdepth (by partitions): 7
#   depth  (by critical paths): 74
#   tdepth (by critical paths): 7
#   Time: 0.001 s
.v x0 x1 x2 y0 y1 y2 y3 y4 y5
.i x0 x1 x2
.o x0 x1 x2 y0 y1 y2 y3 y4 y5

BEGIN
H y4
tof y4 y0
tof y4 y2
tof y4 y3
tof x2 y2
tof x2 y1
tof x0 y3
tof x0 y1
tof x0 y4
tof x2
tof y1
tof y2
tof y3
T* x0
T y0
T* y4
T x2
P y1
T y1
T y2
T* y3
tof y4 y0
tof x2 y3
tof x2 y2
tof x2 y1
tof y4 y3
tof y4 y2
tof x0 y0
tof x0 y2
tof x0 y1
tof x0 y4
H y4
tof x0
H y3
H y5
tof y5 x0
tof y4 x1
tof y5 x2
tof y3 y0
tof y3 y1
tof y4 y1
tof y4 y2
tof y4 y3
tof y3 y4
tof x1 y4
tof x1 y3
tof x0
T y2
T y0
T* x1
T* y3
T* y1
T y4
T y5
T* x2
T x0
tof y3 y2
tof y4 y2
tof x1 y0
tof y4 y0
tof y4 x1
tof y4 y3
tof x2 y5
tof y5 x0
tof y5 x2
tof y3 y4
tof y3 y1
tof x1 y4
tof x1 y1
T* x0
tof y5 x0
tof x2 x0
H x2
H x1
tof x1 y0
tof x1 y2
tof x2 y3
tof y4
H y0
H y2
tof y0 y4
tof y3 y0
tof x2 y0
tof y3 y1
tof x2 y3
tof y0 y3
tof y4 y1
tof y4 y0
tof y0
tof y1
T x2
T y3
T y4
T y0
T y1
tof y0 x2
tof y2 x2
tof y0 y3
tof y1 y3
tof y2 y3
tof y2 y4
tof y2 y0
tof y0 y1
tof y3 y1
tof x2 y4
T y2
T* y3
T x2
T* y0
T y4
tof x2 y2
tof y0 y3
tof y0 x2
tof y1 y4
tof x2 y0
tof y2 y3
tof y2 y0
P* y4
H y0
H y1
tof y0 y3
tof x1 x2
tof x0
H x1
H x2
tof x2 x0
tof x1 y2
tof x2 y5
tof x1 y3
tof y0 y3
tof x1 y0
tof y3 y0
tof y2 y4
tof x1
tof y4
tof x0
T y3
T* y2
T x2
T* y5
T x1
T* y4
T x0
tof y2 y3
tof x1 y3
tof y4 y3
tof y4 y2
tof y5 x2
tof x1 y4
tof x2 x0
tof x2 y5
tof y2 x1
T* x0
tof x2 x0
tof y5 x0
tof y0 y3
H y5
H x1
tof x1 y5
tof y4 y2
tof y2 y4
tof y4 y2
tof y2 x1
tof x1 y2
tof y2 x1
tof y1 y0
tof y0 y1
tof y1 y0
END
