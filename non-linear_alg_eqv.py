class Alg_eqv:
    def __init__(self, function):
        self.function = function

    def set_derivative(self, derivative):
        self.derivative = derivative

    def bisection(self, left, right, eps):
        print("\n\nCalling bisection method for interval [" + str(left) + ", " + str(right) +"].\n")
        f = self.function
        iter = 1
        while right - left > eps or abs(f((left + right)/2)) > eps:
            mid = (right + left)/2
            if f(mid) == 0:
                break
            if f(left)*f(mid) < 0:
                right = mid
            else:
                left = mid
            print("Iteration " + str(iter) +".\t Sector: [" + "%.6f" % left + ", " + "%.6f" % right + "],\t \t f(" + "%.6f" % left +") = "  + "%.6f" % f(left) +", \t f(" + "%.6f" % right + ") = " + "%.6f" % f(right))
            iter += 1
        return (right + left)/2

    def newton(self, x_0, eps):
        print("\n\nCalling Newton method with first approximation " + str(x_0) + ".\n")
        f, df = self.function, self.derivative
        approximations = [x_0]
        iter = 1
        while True:
            x_k = approximations[-1]
            approximations.append(x_k - f(x_k)/df(x_k))
            if abs(x_k - approximations[-1]) < eps and abs(f(approximations[-1])) < eps:
                break
            print("Iteration " + str(iter) +".\t Approximation:  = "  + "%.6f" % approximations[-1] +", \t\t f(" + "%.6f" % approximations[-1] +  ") = %.6f" % f(approximations[-1]))
            iter += 1
        return approximations[-1]

    def chord(self, left, right, eps):
        print("\n\nCalling chord method for interval [" + str(left) + ", " + str(right) +"].\n")
        f = self.function
        approximations = [left]
        iter = 1
        while True:
            x_k = (left * f(right) - right* f(left))/(f(right) - f(left))
            if abs(approximations[-1] - x_k) < eps and abs(f(x_k)) < eps:
                break
            if f(left)*f(x_k) < 0:
                right = x_k
            else:
                left = x_k
            print("Iteration " + str(iter) +".\t Sector: [" + "%.6f" % left + ", " + "%.6f" % right + "],\t \t f(" + "%.6f" % left +") = "  + "%.6f" % f(left) +", \t f(" + "%.6f" % right + ") = " + "%.6f" % f(right))
            approximations.append(x_k)
            iter += 1
        return x_k

def polynom(x):
    return x**5-3*x**4+7*x**2-9
def dpoly(x):
    return 5*x**4-12*x**3+14*x

f = Alg_eqv(polynom)
f.set_derivative(dpoly)

print("Bisection method result:", "%.6f" % f.bisection(1, 2, 0.00001))
print("Chord method result:", "%.6f" % f.chord(1, 2, 0.00001))
print("Newton method result:", "%.6f" % f.newton(2, 0.00001))
