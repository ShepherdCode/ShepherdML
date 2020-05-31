'''
Softmax.
Performs multinomial logistic regression.
Implemented as multiclass shallow NN.
Ouputs are probabilities 0 to 1.
One output per class r; these sum to 1.
Output is negative log-likelihood based on posterior probability.
Loss function is cross-entropy: neg log of true class.
One hidden layer, one weight vector Wr per node, linear activation.
One output node per r with softmax activation.
Softmax: P(r|Xi)=e^(Wr*Xi)/sum(e^(WXi)) where W represents all Wr.
See Charu Aggarwal chapter 2 exercise 6, and figure 2.5(c).
'''
import argparse
import traceback
import sys
import csv

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
    def setup(self,w1,w2,alpha,epochs):
        self.alpha=alpha
        self.epochs=epochs
        self.out.set_weights(w1,w2)
    def show_weights(self):
        return self.out.get_weights()
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
        #say("yhat = "+str(yhat))
        loss = y-yhat
        adj_loss = self.alpha*loss
        #say("adj loss = "+str(adj_loss))
        w1_delta = adj_loss*x1
        w2_delta = adj_loss*x2
        #say("delta = "+str((w1_delta,w2_delta)))
        self.out.adjust_weights(w1_delta,w2_delta)
        #say("adjusted W = "+str(self.out.get_weights()))
    def train(self,csvfilename):
        for i in range(0,self.epochs):
            # TO DO: load the file rather than re-open it
            with open (csvfilename,"r") as csvfile:
                reader = csv.reader(csvfile,delimiter=',')
                for line in reader:
                    (x1,x2,y)=line
                    x1=float(x1)
                    x2=float(x2)
                    y=int(y)
                    self.train_one(x1,x2,y)
    def classify(self,csvfilename):
        with open (csvfilename,"r") as csvfile:
            reader = csv.reader(csvfile,delimiter=',')
            for line in reader:
                (x1,x2)=line
                x1=float(x1)
                x2=float(x2)
                yhat=self.classify_one(x1,x2)
                say("Point (%d,%d) classified %d"%(x1,x2,yhat))

def say(statement):
    try:
        if (args.debug):
            print(statement, file=sys.stderr, end="\n")
    except Exception:
        pass

def args_parse():
    global args
    parser = argparse.ArgumentParser(description='Linear classifier.')
    parser.add_argument('--train', help = 'Labeled set filename', type=str)
    parser.add_argument('--classify', help = 'Unlabeled set filename', type=str)
    parser.add_argument('--alpha', help= 'Learn rate (1)', type=float, default=float(1))
    parser.add_argument('--debug', action = 'store_true')
    args = parser.parse_args()  # on error, program exits

if __name__ == '__main__':
    """
    Command line invocation:
    $ python3 Softmax.py --help
    """
    try:
        args_parse()
        p = Perceptron()
        epochs=3 # TO DO: user parameter
        p.setup(0,0,args.alpha,epochs)
        if args.train is not None:
            p.train(args.train)
            print("Trained weights: "+str(p.show_weights()))
        if args.classify is not None:
            p.classify(args.classify)
    except Exception as e:
        print("\nThere was an error.")
        if args.debug:
            print(traceback.format_exc())
        else:
            print ('Consider running with --debug')
