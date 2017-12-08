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
    grn_space = util.GrnSpace(gene_list, regulators,space_size=20)

    # init selected grn to null
    selected_grn = None

    # least mse
    least_mse = 0.5

    # loop over all the grns
    while grn_space.has_next():

        # get the grn
        grn = grn_space.get()

        # train the grn
        trained_grn = util.train_grn(grn,training_expression)

        # simulate the network to predict expression
        predicted_expression = util.predict_expression(trained_grn,training_expression,testing_expression.shape[1])

        # calculate the mse
        mse = util.calculate_mse(predicted_expression,testing_expression)
        print("MSE is : ",mse)

        # check if mse is lower than or equal to previous value
        if mse <= least_mse:
            # select this grn
            selected_grn = trained_grn

            # update the mse value
            least_mse = mse

    if selected_grn:
        matrix = (nx.to_numpy_matrix(selected_grn.graph))
        nx.draw(selected_grn.graph,with_labels=True)

        print("yeah we found a grn")

        predicted_expression = util.predict_expression(selected_grn,training_expression,41)
        predicted_expression = predicted_expression.T
        predicted_expression.plot()

        testing_expression = testing_expression.T
        testing_expression.plot()

        plt.draw()
        plt.show()

    else:
        print("Sorry could not find a grn")





