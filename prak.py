from enum import Enum
import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation

class Grid:
    h = 0
    tau = 0
    l = 0
    def __init__(self, h, tau, l):
       self.h, self.tau, self.l = h, tau, l

    def getH(self):
        return self.h

    def getL(self):
        return self.l

    def getTau(self):
        return self.tau

    def getStepCount(self):
        return math.ceil(self.l / self.h)


class ConditionMode(Enum):
    NormalFunction = 0,
    DeriviativeFunction = 1

class Condition:
    def __init__(self, Condition, mode = ConditionMode):
        self.Condition, self.mode = Condition, mode

class Task:
    def __init__(self, l, a, StartCondition, BoundaryConditions):
        assert isinstance(StartCondition, Condition)
        self.l, self.a  = l, a
        self.StartCondition, self.BoundaryConditions = StartCondition, BoundaryConditions
    def getL(self):
        return self.l

    def getA(self):
        return self.a

    def getStartConditions(self):
        return self.StartCondition

    def getBoundaryConditions(self):
        return self.BoundaryConditions

def buildGenerator(condition, stepCount, step):
    assert isinstance(condition, Condition)
    prevElement = 0
    for i in range(stepCount):
        if (condition.mode == ConditionMode.DeriviativeFunction):
            currentElement =  prevElement + step * condition.Condition(i * step)
            yield currentElement
            prevElement = currentElement
        else:
            yield condition.Condition(i * step)


def primeCalc(prevU, h, tau, a):
    currentU = [0] * len(prevU)
    for j in range(1, len(prevU) - 1):
        currentU[j] = prevU[j] + a * tau * (prevU[j - 1] - 2 * prevU[j] + prevU[j + 1])/(h * h)

    return currentU


def Solve(task, grid, iterationCount, recalcFunc):
    assert isinstance(task, Task)
    assert (task.getL() == grid.getL())

    prevU = [x for x in buildGenerator(task.getStartConditions(), grid.getStepCount(), grid.getH())]

    xx = []
    for i in range(grid.getStepCount()):
        xx.append(i * grid.getH())


    currentU = prevU
    for i in range(iterationCount):
        print(currentU)
        yield xx, currentU
        prevU = [x for x in currentU]
        h = grid.getH()
        tau = grid.getTau()
        a = task.getA()
        #ux,t+1 - ux,t/tau = a*(ux-1,t -2ux,t + ux+1,t)/h**2
        ##curentU[grid.getStepCount() - 1] =
        currentU = recalcFunc(prevU, h, tau, a)

        if (task.getBoundaryConditions()[0].mode == ConditionMode.DeriviativeFunction):
            currentU[0] = currentU[1] - grid.getH() * task.getBoundaryConditions()[0].Condition(i * tau)
        else:
            currentU[0] = task.getBoundaryConditions()[0].Condition(i * tau)
        sc = grid.getStepCount()
        if (task.getBoundaryConditions()[1].mode == ConditionMode.DeriviativeFunction):
            currentU[sc - 1] = currentU[sc - 2] + task.getBoundaryConditions()[1].Condition(i * tau) * grid.getH()
        else:
            currentU[sc - 1] = task.getBoundaryConditions()[1].Condition(i * tau)


def shuttleCalc(prevU, h, tau, a):
    #ux,t+1 - ux,t/tau = a*(ux-1,t+1 -2ux,t+1 + ux+1,t+1)/h**2
    #-a/h**2(ux-1,t+1 + ux+1,t+1) + (1/tau +2a/h**2) * ux,t+1 = 1/tau * ux,t

    ## -a/h**2 (1/tau +2a/h**2) -a/h**2
    a_orig = -a/h**2
    b_orig = (1/tau + 2 * a / h**2)
    c_orig = a_orig
    n = len(prevU)

    currentU = [0] * n
    a_c = [0] * n
    b_c = [0] * n
    c_c = [0] * n
    d_c = [0] * n

    e = a_orig
    c = a_orig
    d = b_orig

    a_c[1] = -e / d
    b_c[1] = prevU[0] / d
    for i in range(2, n):
        a_c[i] = -e/(d + c  * a_c[i - 1])
        b_c[i] = (-c* b_c[i - 1] + prevU[i - 1]) / (d + c * a_c[i - 1])


    currentU[n - 1] = (-c * b_c[n - 1] + prevU[n - 1]) / (d + c * a_c[n - 1])
    for i in range(n - 2, -1, -1):
        currentU[i] = a_c[i + 1] * currentU[i + 1] + b_c[i + 1]

    return currentU

def startCondition(x):
    return x

def boundaryCondition0(x):
    return 0

def boundaryConditionL(x):
    return 0



fig, ax = plt.subplots()
line, = ax.plot([], [], lw=2)
xdata, ydata = [], []
ax.set_ylim(-1.1, 1.1)
ax.set_xlim(0, 5)
ax.grid()
ax.set_xlabel("X")
ax.set_ylabel("Temperature")
def run(data):
    # update the data
    t,y = data
    xdata = t
    ydata = y
    xmin, xmax = ax.get_xlim()

    if t[-1] >= xmax:
         ax.set_xlim(xmin, 2*xmax)
         ax.figure.canvas.draw()


    ymin, ymax = ax.get_ylim()
    y = list(y)
    try:
        curYmin = y[0]
        curYmax = y[0]
    except:
        curYmin = 0
        curYmax = 0.5
    for i in range(len(y)):
        curYmax = max(curYmax, y[i])
        curYmin = min(curYmin, y[i])
    if curYmax >= ymax or curYmin <= ymin:
        if (curYmin < 0):
            curYmin *= 2
        else:
            curYmin /= 2

        if (curYmax < 0):
            curYmax /= 2
        else:
            curYmax *= 2

        ax.set_ylim(curYmin, curYmax)
        ax.figure.canvas.draw()

    line.set_data(xdata, ydata)

def main():
    l = 5
    a = 1
    startConditionWrapper = Condition(startCondition,ConditionMode.NormalFunction)
    boundCondition0Wrapper = Condition(boundaryCondition0, ConditionMode.NormalFunction)
    boundConditionLWrapper = Condition(boundaryConditionL, ConditionMode.NormalFunction)


    ## if we use primeCalc (a * tau) /h < 1/2
    task = Task(l, a, startConditionWrapper , [boundCondition0Wrapper, boundConditionLWrapper])
    h = 0.7
    tau = 0.2
    grid = Grid(h, tau, task.getL())

    ani = animation.FuncAnimation(fig, run, Solve(task, grid, 100, shuttleCalc), blit=False, interval=1000, repeat=False)
    plt.show()



if __name__ == "__main__" :
    main()

