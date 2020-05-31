class BorderArray():
    def __init__(self,t):
        self.comparisons=0
        self.S = "$"+t  # throw-away char at zero
        self.B = []   # border
        self.BB=[]    # strong border by linear algorithm
        self.NL=[]    # strong border by easy non-linear algorithm
        for x in self.S:
            self.B.append(0)  # one zero for each letter in t
            self.BB.append(0)  # one zero for each letter in t
            self.NL.append(0)  # one zero for each letter in t
        self.B[0]= -1  # throw-away value at zero
        self.BB[0]= -1  # required value at zero
        self.NL[0]= -1  # required value at zero
        self.computeB()  # fast, standard border array
        self.computeC()  # modified to count comparisons (overwrites B)
        self.computeStrongNL()
        self.computeStrongL()

    def computeStrongL(self):
        pass

    def computeB(self):
        S=self.S
        B=self.B   # alias
        n = self.getLength()
        B[1]=0
        for i in range(1,n): # inclusive/exclusive i.e. 1 to n-1
            b = B[i]
            while b>0 and S[i+1]!=S[b+1]:
                b=B[b]
            if S[i+1]==S[b+1]:
                B[i+1]=b+1
            else:
                B[i+1]=0

    def compare(self,i,b):
        S=self.S
        self.comparisons += 1
        result=(S[i+1]==S[b+1])
        return result

    def computeC(self):
        S=self.S
        B=self.B  # alias
        n = self.getLength()
        B[1]=0
        for i in range(1,n):  # inclusive/exclusive i.e. 1 to n-1
            b = B[i]
            while b>0 and not self.compare(i,b):
                b=B[b]
            if self.compare(i,b):
                B[i+1]=b+1
            else:
                B[i+1]=0

    def getString(self):
        return self.S[1:]

    def getLength(self):
        return len(self.S)-1 # minus 1 to overlook the added $ sign

    def arrayToString(ary):
        conv=""
        for a in ary:
            if a<0:
                conv = conv + 'n'
            else:
                conv = conv + (str(a))
        return conv

    def getBorderString(self):
        conv=BorderArray.arrayToString(self.B)
        disp = conv[1:] # exclude the throw-away value at zero
        return disp

    def getStrongLString(self):
        conv=BorderArray.arrayToString(self.BB)
        return conv

    def getStrongNLString(self):
        conv=BorderArray.arrayToString(self.NL)
        return conv

    def getComparisons(self):
        return self.comparisons

    def borderDemo(question):
        t=""
        n=0
        for i in range(10):
            if question==3:
                tt=t+"b"  # question AA iii
            elif question==2:
                tt="b"+t  # question AA ii
            else:
                tt=t+"a"  # question AA i
            n=n+1
            ba = BorderArray(tt)
            print("n=%d, S=%s, B=%s, C=%d"%(n,ba.getString(),ba.getBorderString(),ba.getComparisons()))
            t=t+"a"

    def demo():
        T= "abaababaabaab"
        print ("T= "+T)
        print ("Expect B= 0011232345645")
        ba = BorderArray(T)
        print("B=%s"%ba.getBorderString())

    def strongDemo1():
        tt='abaab'
        sba=BorderArray(tt)
        modified='n'+sba.getBorderString()
        print("NonLinear")
        print("P= %s, B= %s, BB=%s"%(sba.getString(),modified,sba.getStrongNLString()))
        tt='ababaccababb'
        sba=BorderArray(tt)
        modified='n'+sba.getBorderString()
        print("NonLinear")
        print("P= %s, B= %s, BB=%s"%(sba.getString(),modified,sba.getStrongNLString()))
        tt='abaaaaaaaa'
        sba=BorderArray(tt)
        modified='n'+sba.getBorderString()
        print("NonLinear")
        print("P= %s, B= %s, BB=%s"%(sba.getString(),modified,sba.getStrongNLString()))

        #print("Linear")
        #print("P= %s, B= %s, BB=%s"%(sba.getString(),modified,sba.getStrongLString()))

    def computeStrongX(self):
        n = self.getLength()
        S=self.S
        B=self.B   # alias to border array
        BB=self.BB   # alias to strong border array
        B[0]=-1
        B[1]=0
        BB[1]=-1
        BB[1]=0
        for i in range(1,n): # inclusive/exclusive i.e. 1 to n-1
            b = B[i]
            while b>0 and S[i+1]!=S[b+1]:
                b=B[b]
            if S[i+1]==S[b+1]:
                B[i+1]=b+1
            else:
                B[i+1]=0
        BB[m]=B[m]

    def computeStrongY(self):
        # Linear-time algorithm
        # When computing B, add the check that subsequent characters don't match
        m = self.getLength()
        P=self.S
        B=self.B   # alias to border array
        BB=self.BB   # alias to strong border array
        BB[0]=-1
        BB[1]=0
        BB[m]=B[m]    # requires border array!
        for i in range(2,m): # inclusive/exclusive
            b = B[i]
            while b>0 and P[i+1]==P[b+1]:
                b=B[b]
            if b>=0 and P[i+1]!=P[b+1]:
                BB[i]=b
            else:
                BB[i]=-1

    def computeStrongZ(self):
        # Linear-time algorithm
        # When computing B, add the check that subsequent characters don't match
        # Suspect this is wrong. Need to use border array more.
        P=self.S
        BB=self.BB   # alias
        m = self.getLength()
        BB[0]= -1
        BB[1]=0
        for i in range(1,m-1): # inclusive/exclusive i.e. 1 to m-2
            b = max(0,BB[i])
            while b>0 and (P[i+1]!=P[b+1] or P[i+2]==P[b+2]):
                b=BB[b]
            if P[i+1]==P[b+1] and P[i+2]!=P[b+2]:
                BB[i+1]=b+1
            else:
                BB[i+1]=-1
        B=self.B
        BB[m]=B[m]    # requires border array!

    def computeStrongNL(self):   # first try was correct but non-linear
        P=self.S
        m = self.getLength()
        B=self.B
        BB=self.NL  # alias
        BB[0]= -1
        BB[1]= 0
        BB[m]=B[m]    # requires border array
        for j in range(2,m):  # exclusive of m
            b = B[j]   # requires border array
            max = b + 1
            min = max
            while b>=0:
                if P[j+1]!=P[b+1]: # strong?
                    min=b
                b = B[b]
            if min < max: # found a strong border
                BB[j]=min
            else:
                BB[j]= -1

if __name__ == '__main__':
    #print("--- Basic Demo ---")
    #BorderArray.demo()
    #print("--- Demo 1 ---")
    #BorderArray.borderDemo(1)
    #print("--- Demo 2 ---")
    #BorderArray.borderDemo(2)
    #print("--- Demo 3 ---")
    #BorderArray.borderDemo(3)
    print("--- Strong Demo ---")
    BorderArray.strongDemo1()
