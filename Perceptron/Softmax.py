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
#from numpy import exp
from math import exp

class Vector():
    def __init__(self,dimensions):
        self.dimensions=dimensions
        self.vec = []
        for i in range(0,dimensions):
            self.vec.append(0.1)  # initial guess, get overflow if set >=1
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
        self.numerator=0.0
        self.denominator=0.0
        self.name=""
    def add(self,value,same_class):
        try:
            if same_class:
                self.numerator=exp(value)
            self.denominator += exp(value)
        except:
            print("Exception")
            print("Given value = %f",value)
            raise Exception("Math error")
    def get_output(self):
        prob = self.numerator / self.denominator
        return prob
    def set_class(self,name):
        self.name = name
    def get_class(self):
        return self.name
class Hidden_Node ():
    def __init__(self):
        self.sum=0
    def add(self,value):
        self.sum += value
    def get_output(self):
        activation = 1*self.sum   # linear activation
        return activation
class Multinomial_Logistic_Regression ():
    def __init__(self,epochs=3, alpha=1, dimensions=5, classes=3):
        self.epochs=epochs
        self.alpha=alpha
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
    def set_class_names(self,instances):
        c = -1
        current_class_name = ""
        for instance in instances:
            given_class_name = instance[self.dimensions+1]
            if (given_class_name != current_class_name):
                c += 1
                output_node = self.output_nodes[c]
                output_node.set_class(given_class_name)
                current_class_name = given_class_name
    def load_instances(self,filename):
        instances=[]
        with open (filename,"r") as csvfile:
            reader = csv.reader(csvfile,delimiter=',')
            reader.__next__() # skip header
            for instance in reader:
                instances.append(instance)
        return instances
    def train_file(self,filename):
        say("TRAIN")
        instances=self.load_instances(filename)
        say("Load class names for %d instances"%len(instances))
        self.set_class_names(instances)
        for epoch in range(0,self.epochs):
            say("epoch %d..."%epoch)
            for instance in instances:
                say("instance ..."+str(instance))
                yhat = self.classify_instance(instance)
                self.update_weights(instance)
    def update_weights(self,instance):
        alpha = self.alpha
        y = instance[self.dimensions+1]  # true label = last field in file
        for r in range(0,self.classes):
            weight_vector = self.weight_vectors[r]
            output_node = self.output_nodes[r]
            prob_of_r = output_node.get_output()
            nodename = output_node.get_class()
            for d in range(0,self.dimensions):
                wr = weight_vector.get(d)
                xd = float(instance[d+1])  # TO DO: create a class for instance!
                if (y==nodename):
                    wr = wr + alpha*xd*(1.0-prob_of_r)
                else:
                    wr = wr - alpha*xd*(prob_of_r)
                weight_vector.set(d,wr)
    def classify_file(self,filename):
        say("CLASSIFY")
        say("Assume weights are trained or initialized")
        instances = self.load_instances(filename)
        for instance in instances:
            prediction=self.classify_instance(instance)
            print("input: "+str(instance)+" prediction: "+prediction)
            if (args.debug):
                self.show_all_output_values()
    def show_all_output_values(self):
        for c in range(0,self.classes):
            output_node = self.output_nodes[c]
            prob = output_node.get_output()
            classname = output_node.get_class()
            print("%10.7f %s"%(prob,classname))
    def classify_instance(self,instance):
        say("Classify instance "+instance[0])
        for i in range(0,self.dimensions):
            xi = float(instance[i+1])
            self.input_layer.set(i,xi)
        say("Compute hidden layer values ")
        for r in range(0,self.classes):
            hidden_node = self.hidden_nodes[r]
            weight_vector = self.weight_vectors[r]
            for d in range(0,self.dimensions):
                xi = self.input_layer.get(d)
                wi = weight_vector.get(d)
                wx = wi * xi
                hidden_node.add(wx)
        say("Compute output layer values ")
        for c in range(0,self.classes):
            output_node = self.output_nodes[c]
            for h in range (0,self.classes):
                hidden_node = self.hidden_nodes[h]
                value = hidden_node.get_output()
                same_class=(c==h)
                output_node.add(value,same_class)
        best_class = -1
        max_prob = 0
        for c in range(0,self.classes):
            output_node = self.output_nodes[c]
            prob = output_node.get_output()
            if (prob > max_prob):
                best_class = c
        best_class_name = self.get_class_name(best_class)
        return best_class_name
    def get_class_name(self,classnum):
        output_node = self.output_nodes[classnum]
        classname = output_node.get_class()
        return classname

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
        writer.writerow([2,1.3,1.0,1.2,1.1,1.4,"one"])
        writer.writerow([3,2.3,2.0,2.2,2.1,2.4,"two"])
        writer.writerow([4,2.3,2.0,2.2,2.1,2.4,"two"])
        writer.writerow([5,3.3,3.0,3.2,3.1,3.4,"three"])
        writer.writerow([6,3.3,3.0,3.2,3.1,3.4,"three"])
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
        # TO DO: get dimensions and classes from input file
        # TO DO: get epochs and learn rate from user
        nn = Multinomial_Logistic_Regression (1,1,5,3)
        if args.example is not None:
            create_sample_data(args.example)
        if args.train is not None:
            nn.train_file(args.train)
        if args.classify is not None:
            nn.classify_file(args.classify)
    except Exception as e:
        print("\nThere was an error.")
        if args.debug:
            print(traceback.format_exc())
        else:
            print ('Consider running with --debug')
