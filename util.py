import pandas
import networkx as nx
import random

def load_original_expression():
    dataframe = pandas.read_csv('data/original/ecoli_10_4/Ecoli_subnet_10_4-1_dream4_timeseries.tsv', sep = '\t', header=0)
    return dataframe

def extract_training_data(original_expression, percent):
    normalized_percent = percent * 0.01
    data_size = original_expression.Time.count() * normalized_percent
    dataframe = original_expression[:int(data_size)]
    return dataframe


def extract_testing_data(original_expression, percent):
    normalized_percent = percent * 0.01
    data_size = original_expression.Time.count() * normalized_percent
    dataframe = original_expression.tail(int(data_size))
    return dataframe


def load_gene_list(expression_data):
    return list(expression_data.columns.values)[1:]


def train_grn(grn, training_expression):
    return grn


def predict_expression(trained_grn):
    return 1


def calculate_mse(predicted_expression, testing_expression):
    return 10**-4


class GrnSpace:

    def __init__(self,gene_list,max_regulator_count,space_size = 1000):
        self.gene_list = gene_list
        self.max_regulators_count = max_regulator_count
        self.space_size = space_size
        self.generated_count = 0

    def has_next(self):
        if self.generated_count < self.space_size:
            return True
        else:
            return False

    def get(self):

        graph = self._generate_graph()

        while not nx.is_directed_acyclic_graph(graph):
            print("is cyclic")
            graph = self._generate_graph()

        self.generated_count = self.generated_count + 1
        return graph

    def _generate_graph(self):
        in_degree = []

        for _gene in self.gene_list:
            in_degree.append(random.randint(0, 4))

        out_degree = in_degree[:]
        random.shuffle(out_degree)
        graph = nx.directed_configuration_model(in_degree, out_degree)
        return graph