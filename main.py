import util
import matplotlib.pyplot as plt
import networkx as nx
import tqdm
from multiprocessing import Pool
import numpy

# load the expression data
original_expression = util.load_original_expression()

# split original expression in training and testing data
training_expression = util.extract_training_data(original_expression, 80)
testing_expression = util.extract_testing_data(original_expression, 20)

# set the max allowed regulators to 4
regulators = 4

# maximum iteration for bapso
maxiter = 1

# initialize the grn space
grn_space_size = 1

# least mse
least_mse = 0.9



def train_grn(grn):

    # train the grn
    trained_grn = util.train_grn(grn, training_expression, maxiter)

    # simulate the network to predict expression
    predicted_expression = util.predict_expression(trained_grn, training_expression, testing_expression.shape[1])

    # calculate the mse
    mse = util.calculate_mse(predicted_expression, testing_expression)

    return (mse,trained_grn)


if __name__ == "__main__":

    # load the gene list
    gene_list = util.load_gene_list(original_expression)

    grn_space = util.GrnSpace(gene_list, regulators,space_size=grn_space_size)

    # init selected grn to null
    selected_grn = None


    grns = []

    pbar = tqdm.tqdm(total=grn_space_size)

    # loop over all the grns
    while grn_space.has_next():

        # get the grn
        grn = grn_space.get()

        grns.append(grn)

        pbar.update(1)

    pbar.close()


    # do the task
    trained_grns = []
    pool = Pool(processes=8)
    for value in tqdm.tqdm(pool.imap_unordered(train_grn,grns), total=len(grns)):
        trained_grns.append(value)

    # loop over the trained grns
    pbar = tqdm.tqdm(total=grn_space_size)
    for value in trained_grns:
        mse,trained_grn = value
        if mse < least_mse:
            least_mse = mse
            selected_grn = trained_grn
        pbar.update(1)
    pbar.close()

    if selected_grn:
        matrix = (nx.to_numpy_matrix(selected_grn.graph))

        data = util.serialize_predictions(matrix,selected_grn.genes)

        data.to_csv('data\data.csv',header=False,sep='\t',index=False)

        nx.draw(selected_grn.graph,with_labels=True)

        print("yeah we found a grn")
        print("mse is : ",least_mse)

        predicted_expression = util.predict_expression(selected_grn,training_expression,41)
        predicted_expression = predicted_expression.T
        predicted_expression.plot()

        testing_expression = testing_expression.T
        testing_expression.plot()

        plt.draw()
        plt.show()

    else:
        print("Sorry could not find a grn")









