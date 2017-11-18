import numpy as np



'''
full duplex
'''

class Result(object):
    def __init__(self, timeSlots, throughputRate):
        self.timeSlots = timeSlots
        self.throughputRate = throughputRate

class Solution(object):
    def __init__(self, Dj, numofue, eh):
        self.numOfUser = numofue
        self.H = eh
        self.minD = Dj
        self.lamb = 1.0
        self.mu = np.ones(self.numOfUser, dtype=float)
        self.z = np.ones(self.numOfUser, dtype=float)
        self.step = 0.006
        self.error = 0.005
        self.tau = [0.0] * (self.numOfUser + 1)
        self.throughput = np.zeros(self.numOfUser, dtype=float)

    def findMu(self):
        count = 0
        grad = []
        while count < 3000:

            self.findLambda()
            tmp = [1] * (self.numOfUser + 1)
            for i in range(1, self.numOfUser + 1):
                tmp[i] = sum(tmp[:i]) * self.H[i - 1] / self.z[i - 1]
            self.tau[0] = 1 / (sum(tmp))
            for i in xrange(1, self.numOfUser + 1):
                self.tau[i] = self.tau[0] * tmp[i]

            t = np.array(self.tau[1:])
            self.throughput = t * np.log(1 + self.z)
            gradient = self.minD - self.throughput
            if (np.sqrt(sum(self.mu * (gradient ** 2))) < self.error):
                break

            grad.append(np.sqrt(sum(self.mu * gradient ** 2)))
            self.mu += self.step * gradient
            for i in range(len(self.mu)):
                if self.mu[i] < 0:
                    self.mu[i] = 0
            count += 1
        self.grad = grad
        # plt.plot(grad)
        # plt.show()
        # print self.mu
        # print self.throughput
        # print self.tau
        if grad[-1] > 0.5:
            return None
        return Result(self.tau, self.throughput)

    def findLambda(self):
        la_max = 1000
        la_min = 0
        while la_max - la_min > 0.05:
            self.lamb = (la_max + la_min) * 0.5
            self.findZ()
            #print (1 + self.mu) * self.H / (1 +self.z)

            if sum((1 + self.mu) * self.H / (1 + self.z)) > self.lamb:
                la_min = self.lamb
            else:
                la_max = self.lamb

    def findZ(self):
        for i in range(self.numOfUser):
            z_max = 1000000
            z_min = 0
            while z_max - z_min > 0.005:
                z_mid = (z_max + z_min) * 0.5
                if np.log(1 + z_mid) - z_mid / (1 + z_mid) > self.lamb / (1 + self.mu[i]):
                    z_max = z_mid
                else:
                    z_min = z_mid

            self.z[i] = z_mid

if __name__ == "__main__":
    s = Solution([0.1,0.2,0.1,0.4],4,[2.4,2.6,2.1,2.7])
    res = s.findMu()
