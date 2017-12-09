import util
import matplotlib.pyplot as plt
import networkx as nx
import tqdm
from multiprocessing import Pool
import time



def train_grn(training_config):

    grn = training_config.grn
    training_expression = training_config.training_expression
    testing_expression = training_config.testing_expression
    maxiter = training_config.maxiter
    pertubations = training_config.pertubation

    # train the grn
    trained_grn = util.train_grn(grn, training_expression, maxiter,pertubations)

    # simulate the network to predict expression
    predicted_expression = util.predict_expression(trained_grn, training_expression, testing_expression.shape[1])

    # calculate the mse
    mse = util.calculate_mse(predicted_expression, testing_expression)

    return (mse,trained_grn)


if __name__ == "__main__":

    # grn search space
    grn_space_size = 1000

    # maximum iteration for bapso
    maxiter = 20

    # mse to use as cutoff
    least_mse = 0.1

    # set the max allowed regulators to 4,
    regulators = 4                              # will change according to selection
    num_time_series = 25                        # will change according to selection
    time_points_per_series = 21                 # constant

    selected_option = 0
    while not 1 <= selected_option <= 10 :
        print("1. Ecoli 5 genes 2 regulators 4 interactions 25 timepoints")
        print("2. Ecoli 5 genes 2 regulators 4 interactions 50 timepoints")
        try:
            selected_option = int(input("Select the dataset to use [1] : "))
        except:
            selected_option = 1

        try:
            maxiter = int(input("Enter the maximum iteration you want in BAPSO [20] : "))
        except:
            maxiter = 20

        try:
            grn_space_size = int(input("Enter the size of GRN search space [1000] : "))
        except:
            grn_space_size = 1000

        try:
            least_mse = float(input("Enter the mse you want to use as cutoff [0.1] : "))
        except:
            least_mse = 0.1


    # get the data
    original_expression,regulators,num_time_series,pertubations = util.load_original_expression(selected_option)
    training_expression = util.extract_training_data(original_expression,num_time_series,time_points_per_series)
    testing_expression = util.extract_testing_data(original_expression,num_time_series,time_points_per_series)

    # load the gene list
    gene_list = util.load_gene_list(original_expression)

    # create the grn space
    grn_space = util.GrnSpace(gene_list, regulators,space_size=grn_space_size)

    # init selected grn to null
    selected_grn = None

    # list to hold all the training configs
    grn_training_configs = []

    pbar = tqdm.tqdm(total=grn_space_size)

    # loop over all the grns
    while grn_space.has_next():

        # get the grn
        grn = grn_space.get()

        training_config = util.TrainingConfig(grn,training_expression,testing_expression,maxiter,pertubations)

        grn_training_configs.append(training_config)

        pbar.update(1)

    pbar.close()

    # do the task
    trained_grns = []
    pool = Pool(processes=8)
    for value in tqdm.tqdm(pool.imap_unordered(train_grn, grn_training_configs), total=len(grn_training_configs)):
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

        # get the path to save this file based on the option selected
        path_to_save = util.get_path_to_save(selected_option)

        # save the data
        data.to_csv(path_to_save, header=False, sep='\t', index=False)

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
        time.sleep(0.1)
        print("Sorry could not find a grn")









