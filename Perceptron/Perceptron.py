'''
Perceptron.
See Charu Aggarwal chapter 1.
'''
import argparse
import traceback
import sys

class OutputLayer ():
    def __init__(self):
        self.output=0
    def set_weights(self,w1,w2):
        self.w1=w1
        self.w2=w2
    def get_weights(self):
        return (self.w1,self.w2)
    def adjust_weights(self,w1_delta,w2_delta):
        self.w1 += w1_delta
        self.w2 += w2_delta
class Perceptron ():
    def __init__(self):
        self.out = OutputLayer()
        self.alpha = 0
    def setup(self,w1,w2,alpha):
        self.alpha=alpha
        self.out.set_weights(w1,w2)
    def classify_one(self,x1,x2):
        (w1,w2)=self.out.get_weights()
        sum = (w1*x1)+(w2*x2)
        yhat = self.activation(sum)
        return yhat
    def activation(self,sum):
        sign = -1
        if (sum>0):
            sign = 1
        return sign
    def train_one(self,x1,x2,y):
        #say("initial W = "+str(self.out.get_weights()))
        yhat = self.classify_one(x1,x2)
        say("yhat = "+str(yhat))
        loss = y-yhat
        adj_loss = self.alpha*loss
        say("adj loss = "+str(adj_loss))
        w1_delta = adj_loss*x1
        w2_delta = adj_loss*x2
        say("delta = "+str((w1_delta,w2_delta)))
        self.out.adjust_weights(w1_delta,w2_delta)
        say("adjusted W = "+str(self.out.get_weights()))

def say(statement):
    try:
        if (args.debug):
            print(statement, file=sys.stderr, end="\n")
    except Exception:
        pass

def args_parse():
    global args
    parser = argparse.ArgumentParser(description='Linear classifier.')
    parser.add_argument('--train', help = 'Learn parameters', action = 'store_true')
    parser.add_argument('--run', help = 'Classify data', action = 'store_true')
    parser.add_argument('--alpha', help= 'Learn rate (1)', type=float, default=float(1))
    parser.add_argument('--debug', action = 'store_true')
    args = parser.parse_args()  # on error, program exits

if __name__ == '__main__':
    """
    Command line invocation:
    $ python3 Perceptron.py --help
    """
    try:
        args_parse()
        p = Perceptron()
        p.setup(0,0,args.alpha)
        if args.train:
            p.train_one(1,1,1)
            p.train_one(1,-1,1)
            p.train_one(-1,1,-1)
            p.train_one(-1,-1,-1)
        if args.run:
            result1=p.classify_one(1,1)
            result2=p.classify_one(1,-1)
            result3=p.classify_one(-1,1)
            result4=p.classify_one(-1,-1)
            print(result1,result2,result3,result4)
    except Exception as e:
        print("\nThere was an error.")
        if args.debug:
            print(traceback.format_exc())
        else:
            print ('Consider running with --debug')
