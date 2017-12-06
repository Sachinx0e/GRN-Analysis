import util
import matplotlib.pyplot as plt
import networkx as nx

if __name__ == "__main__":
    # load the expression data
    original_expression = util.load_original_expression()

    # split original expression in training and testing data
    training_expression = util.extract_training_data(original_expression,80)
    testing_expression = util.extract_testing_data(original_expression,20)

    # load the gene list
    gene_list = util.load_gene_list(original_expression)

    # set the max allowed regulators to 4
    regulators = 4

    # initialize the grn space
    grn_space = util.GrnSpace(gene_list, regulators)

    # init selected grn to null
    selected_grn = None

    # least mse
    least_mse = 10**-3

    # loop over all the grns
    if grn_space.has_next():

        # get the grn
        grn = grn_space.get()

        # train the grn
        trained_grn = util.train_grn(grn,training_expression)

        # simulate the network to predict expression
        predicted_expression = util.predict_expression(trained_grn)

        # calculate the mse
        mse = util.calculate_mse(predicted_expression,testing_expression)

        # check if mse is lower than or equal to previous value
        if mse <= least_mse:
            # select this grn
            selected_grn = trained_grn

            # update the mse value
            least_mse = mse

    if selected_grn:
        nx.draw(selected_grn)
        plt.draw()
        plt.show()
        print("yeah we found a grn")

    else:
        print("Sorry could not find a grn")





