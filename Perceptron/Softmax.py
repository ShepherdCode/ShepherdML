'''
Softmax.
Performs multinomial logistic regression.
Implemented as multiclass shallow NN.
Ouputs are probabilities 0 to 1.
One output per class r; these sum to 1.
Output is negative log-likelihood based on posterior probability.
Loss function is cross-entropy: neg log of true class.
Udate is Wr=Wr+alpha{r is correct ? Xi * (1-P(r|Xi)) : -Xi * P(r|Xi)}.
One hidden layer, one weight vector Wr per node, linear activation.
One output node per r with softmax activation.
Softmax: P(r|Xi)=e^(Wr*Xi)/sum(e^(WXi)) where W represents all Wr.
See Charu Aggarwal chapter 2 exercise 6, and figure 2.5(c).
'''
import argparse
import traceback
import sys
import csv

class Weight_Vector():
    def __init__(self):
        self.vec = []
    def set_weight(self,i,w):
        self.vec[i]=w
class Output_Node ():
    def __init__(self):
        pass
class Hidden_Node ():
    def __init__(self):
        pass
class Multinomial_Logistic_Regression ():
    def __init__(self,epochs=3, alpha=1, dimensions=5, classes=3):
        self.epochs=3
        self.alpha=1
        self.dimensions=dimensions
        self.classes=classes

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
    parser.add_argument('--debug', action = 'store_true')
    args = parser.parse_args()  # on error, program exits

if __name__ == '__main__':
    """
    Command line invocation:
    $ python3 Softmax.py --help
    """
    try:
        args_parse()
        nn = Multinomial_Logistic_Regression (3,1,5,3)
        if args.train is not None:
            nn.train(args.train)
        if args.classify is not None:
            nn.classify(args.classify)
    except Exception as e:
        print("\nThere was an error.")
        if args.debug:
            print(traceback.format_exc())
        else:
            print ('Consider running with --debug')
