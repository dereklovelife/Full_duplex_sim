import numpy as np
import numpy.linalg as la
import scipy.stats as stats

class BeamOptimization(object):
    def __init__(self):
        self.H1 = None ## BS to UEs K*Nt
        self.H2 = None ## UEs to BS Nr*K
        self.HSI = None ## self-interference Nr*Nt

        self.sendBeam = None  ## matrix St
        self.recvBeam = None  ## matrix Sr

        self.sendBeamVector = None ## vector wt
        self.recvBeamVector = None ## vector wr

        self.sumEH = 0 ## value sum_EH
        self.EH = None ## matrix Energy harvested by UEs

    def setH1(self, h):
        self.H1 = np.mat(h)

    def setH2(self, h):
        self.H2 = np.mat(h)

    def setHSI(self, h):
        self.HSI = np.mat(h)

    def getSendBeam(self):
        _, sigma, vh = la.svd(self.H1)
        self.sendBeam = vh[0].H * vh[0]
        self.sendBeamVector = vh[0].H
        self.sumEH = sigma[0] ** 2
        tmp = self.H1 * self.sendBeam * self.H1.H
        tmp = tmp.real
        self.EH = np.mat(np.zeros((len(tmp), len(tmp))))
        for i in xrange(len(tmp)):
            self.EH[i,i] = tmp[i,i] ** 0.5

    def getRecvBeam(self):
        G = self.sendBeam.H * self.HSI.H
        _, _, vGh = la.svd(G)

        Hu = self.EH.H * self.H2.H * vGh.H
        D = Hu[1:, 1:]
        b = Hu[0, 1:]

        newHu = b.T * b + D.T * D
        _,sigma,V = la.svd(newHu)







        self.recvBeam = np.mat(np.zeros((len(self.H2), len(self.H2)), dtype = np.complex))
        self.recvBeam[1:, 1:] = tmp
        ##self.recvBeamVector = np.mat(np.zeros(len(self.H2)), dtype = np.complex)
        ##self.recvBeamVector = self.recvBeamVector.T
        ##for i in xrange(1, len(self.H2)):
            ##self.recvBeamVector[i,0] = (vHg[0].H)[i-1, 0]

        self.recvBeam = vGh.H * self.recvBeam * vGh
        ##self.recvBeamVector = vGh.H * self.recvBeamVector
        h2 = self.EH * self.H2.H

        ret = (h2 * self.recvBeam * h2.H).real
        print np.trace(ret)
        print sigma[0]




if __name__ == "__main__":
    a = stats.norm.rvs(size = (3,3))
    b = stats.norm.rvs(size = (3,3))
    H1 = a + b * np.complex(0,1)
    a = stats.norm.rvs(size = (3,3))
    b = stats.norm.rvs(size = (3,3))
    H2 = a + b * np.complex(0,1)
    a = stats.norm.rvs(size = (3,3))
    b = stats.norm.rvs(size = (3,3))
    HSI = a + b * np.complex(0,1)

    b = BeamOptimization()
    b.setH1(H1)
    b.setH2(H2)
    b.setHSI(HSI)
    b.getSendBeam()
    b.getRecvBeam()

