'''
Softmax.
Performs multinomial logistic regression.
Implemented as multiclass shallow NN.
Ouputs are probabilities 0 to 1.
One output per class r; these sum to 1.
Output is negative log-likelihood based on posterior probability.
Loss function is cross-entropy: neg log prob of true class given Xi.
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

class Vector():
    def __init__(self,dimensions):
        self.dimensions=dimensions
        self.vec = []
        for i in range(0,dimensions):
            self.vec.append(1)  # initial guess
    def set(self,i,w):
        if i>=0 and i<self.dimensions:
            self.vec[i]=w
    def get(self,i):
        if i>=0 and i<self.dimensions:
            return self.vec[i]
    def get_dimension(self):
        return self.dimensions
class Output_Node ():
    def __init__(self):
        pass
class Hidden_Node ():
    def __init__(self):
        self.sum=0
    def add(self,value):
        self.sum += value
    def get_output(self):
        return activation(self.sum)
    def activation(self,value):
        return value # linear/identity activation
class Multinomial_Logistic_Regression ():
    def __init__(self,epochs=3, alpha=1, dimensions=5, classes=3):
        self.epochs=3
        self.alpha=1
        self.dimensions=dimensions
        self.classes=classes
        self.input_layer=Vector(dimensions)
        self.weight_vectors=[]
        for r in range(0,classes):
            self.weight_vectors.append(Vector(dimensions))
        self.hidden_nodes=[]
        for r in range(0,classes):
            self.hidden_nodes.append(Hidden_Node())
        self.output_nodes=[]
        for r in range(0,classes):
            self.output_nodes.append(Output_Node())
    def train(self,filename):
        with open (csvfilename,"r") as csvfile:
            reader = csv.reader(csvfile,delimiter=',')
            for epoch in range(0,self.epochs):
                print("epoch %d"%epoch)
    def classify(self,filename):
        print("Assume weights are trained or initialized")
        instances=[]
        with open (filename,"r") as csvfile:
            reader = csv.reader(csvfile,delimiter=',')
            reader.__next__() # skip header
            for instance in reader:
                instances.append(instance)
        for instance in instances:
            print("Classify instance %d"%int(instance[0]))
            for i in range(0,self.dimensions):
                xi = float(instance[i+1])
                self.input_layer.set(i,xi)
            for r in range(0,self.classes):
                hidden_node = self.hidden_nodes[r]
                weight_vector = self.weight_vectors[r]
                for d in range(0,self.dimensions):
                    xi = self.input_layer.get(d)
                    wi = weight_vector.get(d)
                    wx = wi * xi
                    hidden_node.add(wx)

def say(statement):
    try:
        if (args.debug):
            print(statement, file=sys.stderr, end="\n")
    except Exception:
        pass

def create_sample_data(prefix):
    training_file=prefix+".labeled.csv"
    with open (training_file,"w") as csvfile:
        writer = csv.writer(csvfile,delimiter=',')
        writer.writerow(["n","x1","x2","x3","x4","x5","class"])
        writer.writerow([1,1.3,1.0,1.2,1.1,1.4,"one"])
        writer.writerow([1,1.3,1.0,1.2,1.1,1.4,"one"])
        writer.writerow([2,2.3,2.0,2.2,2.1,2.4,"two"])
        writer.writerow([2,2.3,2.0,2.2,2.1,2.4,"two"])
        writer.writerow([3,3.3,3.0,3.2,3.1,3.4,"three"])
        writer.writerow([3,3.3,3.0,3.2,3.1,3.4,"three"])
    testing_file=prefix+".unlabeled.csv"
    with open (testing_file,"w") as csvfile:
        writer = csv.writer(csvfile,delimiter=',')
        writer.writerow(["n","x1","x2","x3","x4","x5"])
        writer.writerow([1,2.5,2.1,2.1,2.1,2.0])

def args_parse():
    global args
    parser = argparse.ArgumentParser(description='Linear classifier.')
    parser.add_argument('--example', help = 'Prefix for generated files', type=str)
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
        if args.example is not None:
            create_sample_data(args.example)
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
