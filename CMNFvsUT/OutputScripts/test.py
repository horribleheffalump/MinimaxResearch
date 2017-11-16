import numpy

x = numpy.array([3,2,4,0])
y = numpy.array([4,1,3,1])

th0 = -2
th1 = 0.5
def h(x):
    return th0 + th1*x
print(h(6))

print(x.size)
print(1/2/x.size * numpy.sum(numpy.power((h(x)-y),2)))

X = numpy.array([[1, 89, 7921, 96], [1, 72, 5184, 74], [1, 94, 8836, 87], [1, 69, 4761, 78]])
print(X)
print(X[:,2])

av = sum(X[:,2])/X[:,2].size
range = max(X[:,2])-min(X[:,2])
print(av, range,X[3,2] )
print((X[3,2] - av)/range)