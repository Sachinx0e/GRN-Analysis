import util
import matplotlib.pyplot as plt
import networkx as nx
import tqdm
from multiprocessing import Pool
import time
import multiprocessing



def train_grn(training_config):

    grn = training_config.grn
    training_expression = training_config.training_expression
    testing_expression = training_config.testing_expression
    maxiter = training_config.maxiter
    pertubations = training_config.pertubation

    # train the grn
    trained_grn = util.train_grn(grn, training_expression, maxiter)

    # simulate the network to predict expression
    predicted_expression = util.predict_expression(trained_grn, testing_expression,pertubations)

    # calculate the mse
    mse = util.calculate_mse(predicted_expression, testing_expression)

    return (mse,trained_grn)

def run_with_option(selected_options):

    if selected_options.option == 1:
        print("1. Ecoli 5 genes 2 regulators 4 interactions 25 timepoints")
    elif selected_options.option == 2:
        print("2. Ecoli 5 genes 2 regulators 4 interactions 50 timepoints")
    elif selected_options.option == 3:
        print("3. Ecoli 5 genes 2 regulators 7 interactions 25 timepoints")
    elif selected_options.option == 4:
        print("4. Ecoli 5 genes 2 regulators 7 interactions 50 timepoints")
    elif selected_options.option == 5:
        print("5. Ecoli 10 genes 4 regulators 14 interactions 25 timepoints")
    elif selected_options.option == 6:
        print("6. Ecoli 10 genes 4 regulators 14 interactions 50 timepoints")
    elif selected_options.option == 7:
        print("7. Ecoli 10 genes 4 regulators 20 interactions 25 timepoints")
    elif selected_options.option == 8:
        print("8. Ecoli 10 genes 4 regulators 20 interactions 50 timepoints")
    elif selected_options.option == 9:
        print("9. Ecoli 20 genes 4 regulators 24 interactions 25 timepoints")
    elif selected_options.option == 10:
        print("10. Ecoli 20 genes 4 regulators 24 interactions 50 timepoints")
    elif selected_options.option == 11:
        print("11. Ecoli 20 genes 4 regulators 35 interactions 25 timepoints")
    elif selected_options.option == 12:
        print("12. Ecoli 20 genes 4 regulators 35 interactions 50 timepoints")

    option = selected_options.option

    # grn search space
    grn_space_size = selected_options.grn_space

    # maximum iteration for bapso
    maxiter = selected_options.maxiter

    # set the max allowed regulators to 4,
    regulators = 4 # will change according to selection
    num_time_series = 25  # will change according to selection
    time_points_per_series = 21  # constant


    # get the data
    original_expression, regulators, num_time_series, pertubations = util.load_original_expression(option)
    training_expression = util.extract_training_data(original_expression, num_time_series, time_points_per_series)
    testing_expression, testing_pertubations = util.extract_testing_data(original_expression, num_time_series,
                                                                         time_points_per_series, pertubations)

    # load the gene list
    gene_list = util.load_gene_list(original_expression)

    # create the grn space
    grn_space = util.GrnSpace(gene_list, regulators, space_size=grn_space_size)

    # init selected grn to null
    selected_grn = None

    # list to hold all the training configs
    grn_training_configs = []

    pbar = tqdm.tqdm(total=grn_space_size)

    # loop over all the grns
    while grn_space.has_next():
        # get the grn
        grn = grn_space.get()

        training_config = util.TrainingConfig(grn, training_expression, testing_expression, maxiter,
                                              testing_pertubations)

        grn_training_configs.append(training_config)

        pbar.update(1)

    pbar.close()

    # do the task
    trained_grns = []

    cpu_count = multiprocessing.cpu_count()
    pool = Pool(processes=cpu_count)
    for value in tqdm.tqdm(pool.imap_unordered(train_grn, grn_training_configs), total=len(grn_training_configs)):
        trained_grns.append(value)

    # loop over the trained grns
    least_mse = 10
    pbar = tqdm.tqdm(total=grn_space_size)
    for value in trained_grns:
        mse, trained_grn = value
        if mse < least_mse:
            least_mse = mse
            selected_grn = trained_grn
        pbar.update(1)
    pbar.close()

    matrix = (nx.to_numpy_matrix(selected_grn.graph))
    network_data = util.serialize_predictions(matrix, selected_grn.genes)

    # get the path to save this file based on the option selected
    path_to_save_network = util.get_path_to_save(option)
    # save the data
    network_data.to_csv(path_to_save_network, header=False, sep='\t', index=False)

    # save the least mse
    path_to_save_least_mse = util.get_path_to_save_mse(option)
    f = open(path_to_save_least_mse, 'w')
    f.write(str(least_mse))
    f.close()


if __name__ == "__main__":

    selected_options = None
    choice = 0
    while not 1 <= choice <= 13:
        print("1. Ecoli 5 genes 2 regulators 4 interactions 25 timepoints")
        print("2. Ecoli 5 genes 2 regulators 4 interactions 50 timepoints")
        print("3. Ecoli 5 genes 2 regulators 7 interactions 25 timepoints")
        print("4. Ecoli 5 genes 2 regulators 7 interactions 50 timepoints")
        print("5. Ecoli 10 genes 4 regulators 14 interactions 25 timepoints")
        print("6. Ecoli 10 genes 4 regulators 14 interactions 50 timepoints")
        print("7. Ecoli 10 genes 4 regulators 20 interactions 25 timepoints")
        print("8. Ecoli 10 genes 4 regulators 20 interactions 50 timepoints")
        print("9. Ecoli 20 genes 4 regulators 24 interactions 25 timepoints")
        print("10. Ecoli 20 genes 4 regulators 24 interactions 50 timepoints")
        print("11. Ecoli 20 genes 4 regulators 35 interactions 25 timepoints")
        print("12. Ecoli 20 genes 4 regulators 35 interactions 50 timepoints")
        print("13. Run all with maxiter = 20 and grn search space = 1000")


        try:
            choice = int(input("Select the dataset to use [1] : "))
        except:
            choice = 1


        try:
            maxiter = int(input("Enter the maximum iteration you want in BAPSO [20] : "))
        except:
            maxiter = 20

        try:
            grn_space_size = int(input("Enter the size of GRN search space [1000] : "))
        except:
            grn_space_size = 1000

        selected_options = util.SelectedOptions(choice, grn_space_size, maxiter)

    if selected_options.option != 13 :
        run_with_option(selected_options)
    else:
        for i in range(1,12):
            selected_options.option = i
            run_with_option(selected_options)






