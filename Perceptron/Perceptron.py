'''
Perceptron.
See Charu Aggarwal chapter 1.
'''
import argparse
import traceback

class OutputLayer ():
    def __init__(self):
        self.output=0
    def set_weights(self,w1,w2):
        self.w1=w1
        self.w2=w2
    def get_weights(self):
        return (self.w1,self.w2)
class Perceptron ():
    def __init__(self):
        self.out = OutputLayer()
        self.alpha = 0
    def setup(self,w1,w2,alpha):
        self.alpha=alpha
        self.out.set_weights(w1,w2)
    def classify(self,x1,x2):
        (w1,w2)=self.out.get_weights()
        sum = (w1*x1)+(w2*x2)
        return sum

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
        if args.train:
            pass
        elif args.run:
            p.setup(0,0,1)
            result=p.classify(1,1)
            print(result)
    except Exception as e:
        print("\nThere was an error.")
        if args.debug:
            print(traceback.format_exc())
        else:
            print ('Consider running with --debug')
