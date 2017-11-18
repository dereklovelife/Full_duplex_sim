import numpy as np
import numpy.linalg as la
import scipy.stats as stats
import cvxpy as cvx

import matplotlib.pyplot as plt

class Beam(object):
    def __init__(self, H1, H2, HL):
        self.H1 = H1
        self.H2 = H2
        self.HL = HL
        self.w1 = None
        self.w2 = None
        self.w3 = None
        self.w4 = None
        self.ret = []
        self.ret2 = []
        self.ret3 = []
        self.vt = None
        self.preIndex = -1


    def getW1(self):
        ## get beamforming at transmitter

        u, sigma, vt = la.svd(self.H1)
        h = self.H1.dot(vt.T)
        h = h ** 2

        s = np.zeros((len(vt), len(vt)))
        s[0,0] = 1
        self.s = vt.conj().T.dot(s).dot(vt)
        ret = self.H1.dot(s).dot(self.H1.conj().T)
        self.energyReceive = [ret[i,i].real for i in xrange(len(ret))]
        self.H2 = np.diag(np.array(np.sqrt(self.energyReceive))).dot(self.H2)
        self.w1 = self.s
        return self.s


    '''
    full-duplex (no self-interference) Receiving beamforming.
    '''
    def getW2(self):
        self.ret = []
        u, simga, vt = la.svd(self.H2)
        s = np.zeros((len(vt), len(vt)))
        s[0,0] = 1
        tmp = self.H2.dot(vt.conj().T).dot(s).dot(vt).dot(self.H2.conj().T)
        self.ret = [tmp[i,i].real for i in xrange(len(tmp))]
        self.w2 = vt.conj().T.dot(s).dot(vt)
        return self.w2

    def getRetval(self):
        if self.w1 == None:
            self.getW1()
        if self.w2 == None:
            self.getW2()
        return self.ret


    def getRetval2(self):
        if self.w1 == None:
            self.getW1()
        if self.w3 == None:
            self.getW3()
        return self.ret2

    def getRetval3(self):
        if self.w1 == None:
            self.getW1()
        if self.w4 == None:
            self.getW4()
        return self.ret3

    '''
    full-duplex (self-interference) Receiving beamforming.
    '''
    def getW3(self):
        G = self.HL.dot(self.s)
        G = G.conj().T
        _, sigma1, vt = la.svd(G)
        sigma = np.zeros(len(self.HL))
        for i in xrange(len(sigma1)):
            sigma[i] = sigma1[i]
        sigma1 = sigma
        sigma1 = sigma1 ** 2
        h = vt.dot(self.H2.conj().T)
        hh = h.dot(h.conj().T)
        sigma2 = np.array([hh[i,i].real for i in xrange(len(hh))])
        x = cvx.Variable(len(sigma2))
        obj = cvx.Maximize(sigma2 * x)
        cons = [sigma1 * x == 0, x >= 0, cvx.sum_entries(x) <= 1]
        prob = cvx.Problem(obj, cons)
        prob.solve()
        x = x.value
        x = [i[0,0] for i in x]
        self.w3 = np.diag(np.array(x))
        self.w3 = vt.conj().T.dot(self.w3).dot(vt)
        tmp = self.H2.dot(self.w3).dot(self.H2.conj().T)
        self.ret2 = np.array([tmp[i,i].real for i in xrange(len(tmp))])

    def getW33(self, t):
        receivingTime = [0] * (len(t) - 1)
        base = 0.0
        for i in xrange(0, len(t) - 1):
            receivingTime[i] = t[i] + base
            base += t[i]
        receivingTime = np.array(receivingTime)
        transTime = np.array(t[1:])
        G = self.HL.dot(self.s)
        G = G.conj().T
        _, sigma1, vt = la.svd(G)
        h = vt.dot(self.H2.conj().T)
        hh = h.dot(h.conj().T)
        sigma2 = np.array([hh[i,i].real for i in xrange(len(hh))])
        curIndex = -1
        througputMax = 0.0
        ret22max = None
        for i in xrange(1, len(sigma2)):
            x = np.zeros(len(sigma2))
            x[i] = 1
            self.w33 = np.diag(np.array(x))
            self.w33 = vt.conj().T.dot(self.w33).dot(vt)
            tmp = self.H2.dot(self.w33).dot(self.H2.conj().T)
            self.ret22 = np.array([tmp[i,i].real for i in xrange(len(tmp))])
            throughput = np.sum(transTime * np.log(1 + self.ret22 * receivingTime / transTime))
            if throughput > througputMax:
                througputMax = throughput
                curIndex = i
                ret22max = self.ret22

        if self.preIndex == -1:
            if curIndex == -1:
                print "error occurs when choose channel"
            self.preIndex = curIndex

        return ret22max, curIndex

    def getRetval22(self, t):
        if self.w1 == None:
            self.getW1()

        return self.getW33(t)


'''
    def getW4(self):
        G = self.HL.dot(self.s)
        G = G.T
        _, sigma1, vt = la.svd(G)
        sigma = np.zeros(len(self.HL))
        for i in xrange(len(sigma1)):
            sigma[i] = sigma1[i]
        sigma1 = sigma
        sigma1 = sigma1 ** 2
        for ah in self.H2:
            ah = ah.dot(vt.T)
            ah = ah ** 2
            x = cvx.Variable(len(ah))
            obj = cvx.Maximize(ah * x)
            cons = [sigma1 * x == 0, x >= 0, cvx.sum_entries(x) <= 1]
            prob = cvx.Problem(obj, cons)
            prob.solve()
            x = x.value
            self.ret3.append(prob.value)
        self.w4 = np.ones(1)
'''

if __name__ == "__main__":
    numOfUE = 4
    numofTransmitter = 4
    numOfReciver = 4
    ##H = stats.rayleigh.rvs(size = numofTransmitter)
    ##HSI = np.zeros((numOfReciver, numofTransmitter))
    ##for i in xrange(numOfReciver):
        ##HSI[i] = H
    HSI = stats.norm.rvs(size = (numOfReciver, numofTransmitter)) + stats.norm.rvs(size = (numOfReciver, numofTransmitter)) * np.complex(0,1)

    H1 = stats.norm.rvs(size = (numOfUE, numofTransmitter)) + stats.norm.rvs(size = (numOfUE, numofTransmitter)) * np.complex(0,1)
    H2 = stats.norm.rvs(size = (numOfUE, numOfReciver)) + stats.norm.rvs(size = (numOfUE, numOfReciver)) * np.complex(0,1)
    b = Beam(H1, H2, HSI)
    print b.getRetval()
    print sum(b.getRetval())
    print b.getRetval2()
    print sum(b.getRetval2())
    print "Testing: "
    print "func 2:"
    print "Self Interference:"
    print b.H2.dot(b.w2).dot(b.H2.conj().T).real
    print (b.w3.dot(b.HL).dot(b.w1).dot(b.HL.conj().T).dot(b.w3.conj().T)).real
    print b.H2.dot(b.w3).dot(b.H2.conj().T).real
