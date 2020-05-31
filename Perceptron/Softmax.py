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
    def __init__(self,dimension):
        self.dimension=dimension
        self.vec = []
        for i in range(0,dimension):
            self.vec.append(0.1)  # initial guess, get overflow if set >=1
    def set(self,i,w):
        if i>=0 and i<self.dimension:
            self.vec[i]=w
    def get(self,i):
        if i>=0 and i<self.dimension:
            return self.vec[i]
    def get_dimension(self):
        return self.dimension
    def __str__(self):
        return str(self.vec)

class Instance():   # TO DO: extend vector
    def __init__(self,dimension,aslist):
        self.id=aslist[0]
        self.label=None
        self.x=[]
        for d in range(1,dimension+1):
            value=float(aslist[d])
            self.x.append(value)
        if len(aslist)>=dimension+2:
            self.label=aslist[dimension+1]
    def get_label(self):
        return self.label
    def get_id(self):
        return self.id
    def get_value(self,i):
        return self.x[i]
    def __str__(self):
        return str(self.x)

class Output_Node ():
    def __init__(self):
        self.numerator=0.0
        self.denominator=0.0
        self.name=""
    def reset(self):
        self.numerator=0.0
        self.denominator=0.0
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
    def __init__(self,epochs=3, alpha=1.0, dimension=5, classes=3):
        self.epochs=epochs
        self.alpha=alpha
        self.dimension=dimension
        self.classes=classes
        self.input_layer=Vector(dimension)
        self.weight_vectors=[]
        for r in range(0,classes):  # TO DO: enable lambda functions
            self.weight_vectors.append(Vector(dimension))
        self.hidden_nodes=[]
        for r in range(0,classes):
            self.hidden_nodes.append(Hidden_Node())
        self.output_nodes=[]
        for r in range(0,classes):
            self.output_nodes.append(Output_Node())
        self.data_instances=[]  # TO DO: rename because it also holds an unlabeled instance
    def set_class_labels(self):
        all_labels=[]
        for instance in self.data_instances:
            given_label = instance.get_label() # last column in file
            all_labels.append(given_label)
        uniq_labels=sorted(set(all_labels))
        num_labels=len(uniq_labels)
        num_classes=self.classes
        if num_labels != num_classes:
            print("Mismatch in training data: %d labels for %d classes"%
            (num_labels,num_classes))
            raise Exception ("Mismatch in training data")
        for c in range(0,num_classes):
            output_node = self.output_nodes[c]
            this_label = uniq_labels[c]
            output_node.set_class(this_label)
    def load_instances(self,filename):
        self.data_instances=[]
        with open (filename,"r") as csvfile:
            reader = csv.reader(csvfile,delimiter=',')
            reader.__next__() # skip header
            for oneline in reader:
                instance = Instance(self.dimension,oneline)
                self.data_instances.append(instance)
    def train_file(self,filename):
        say("TRAIN")
        if (args.debug):
            self.show_all_weights()
        self.load_instances(filename)
        say("Load class names for %d instances"%len(self.data_instances))
        self.set_class_labels()
        for epoch in range(0,self.epochs):
            say("EPOCH %d..."%epoch)
            for instance in self.data_instances:
                yhat = self.classify_instance(instance)
                say("instance "+str(instance)+" classified as "+yhat)
                if (args.debug):
                    self.show_all_output_values()
                say("update weights...")
                self.update_weights(instance)
                if (args.debug):
                    self.show_all_weights()
    def show_all_weights(self):
        print ("WEIGHTS")
        for r in range(0,self.classes):
            weight_vector = self.weight_vectors[r]
            print(weight_vector)
    def update_weights(self,instance):
        # Assume instance was just classified so output nodes have probabilities
        alpha = self.alpha
        y = instance.get_label()
        say("true class is "+y)
        for r in range(0,self.classes):
            weight_vector = self.weight_vectors[r]
            output_node = self.output_nodes[r]
            prob_of_r = output_node.get_output()
            nodename = output_node.get_class()
            for d in range(0,self.dimension):
                wr = weight_vector.get(d)
                xd = instance.get_value(d)
                if (y==nodename):
                    wr = wr + alpha*xd*(1.0-prob_of_r)
                else:
                    wr = wr - alpha*xd*(prob_of_r)
                weight_vector.set(d,wr)
    def classify_file(self,filename):
        say("CLASSIFY")
        say("Assume weights are trained or initialized")
        self.load_instances(filename)
        for instance in self.data_instances:
            prediction=self.classify_instance(instance)
            print("input: "+str(instance)+" prediction: "+prediction)
            if (args.debug):
                self.show_all_output_values()
    def show_all_output_values(self):
        for c in range(0,self.classes):
            output_node = self.output_nodes[c]
            prob = output_node.get_output()
            classname = output_node.get_class()
            print("%f %s"%(prob,classname))
    def classify_instance(self,instance):
        say("Classify instance "+instance.get_id())
        for c in range(0,self.classes):
            output_node = self.output_nodes[c]
            output_node.reset()
        for i in range(0,self.dimension):
            xi = instance.get_value(i)
            self.input_layer.set(i,xi)
        for r in range(0,self.classes):
            hidden_node = self.hidden_nodes[r]
            weight_vector = self.weight_vectors[r]
            for d in range(0,self.dimension):
                xi = self.input_layer.get(d)
                wi = weight_vector.get(d)
                wx = wi * xi
                hidden_node.add(wx)
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
                max_prob = prob
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
        writer.writerow([1,1.0,1.0,1.1,1.1,1.1,"one"])
        writer.writerow([2,2.1,2.1,2.4,2.1,2.2,"two"])
        writer.writerow([3,3.0,3.1,3.0,3.1,3.2,"three"])
        writer.writerow([4,1.3,1.0,1.2,1.1,1.4,"one"])
        writer.writerow([5,2.3,2.0,2.2,2.1,2.4,"two"])
        writer.writerow([6,3.3,3.0,3.2,3.1,3.4,"three"])
    testing_file=prefix+".unlabeled.csv"
    with open (testing_file,"w") as csvfile:
        writer = csv.writer(csvfile,delimiter=',')
        writer.writerow(["n","x1","x2","x3","x4","x5"])
        writer.writerow([7,2.0,2.1,2.1,2.1,2.0])

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
    $ python3 Softmax.py --example regression
    $ python3 Softmax.py --train regression.labeled.csv --classify regression.unlabeled.csv
    # As yet, no means of saving the trained weights, so train & classify.
    """
    try:
        args_parse()
        # TO DO: get dimension and classes from input file
        # TO DO: get epochs and learn rate from user
        nn = Multinomial_Logistic_Regression (10,0.1,5,3)
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
