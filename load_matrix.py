import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def largest_indices(ary, n):
    flat = ary.flatten()
    indices = np.argpartition(flat, -n)[-n:]
    indices = indices[np.argsort(-flat[indices])]
    return np.unravel_index(indices, ary.shape)


def the_middle(pair):
    if pair[0] in range(pair[1], pair[1]+15):
        return True
    if pair[0] in range(pair[1]-15, pair[1]):
        return True
    else:
        return False


def find_indices(protein_file, how_many):
    with open(protein_file) as json_file:
        data = json.load(json_file)
    probs = data['contact_probs']
    probs_df = pd.DataFrame(probs)
    probs_np = np.array(probs)
    idx = largest_indices(probs_np, how_many)
    idx_list = list(zip(idx[0], idx[1]))
    idx_list = [(int(idx[0]), int(idx[1])) for idx in idx_list]

    new_list = []
    for ind in idx_list:
        if not the_middle(ind):
            new_list.append(ind)

    return probs_df, new_list


def plot(probs_df, new_list):
    xs = [ind[0] for ind in new_list]
    ys = [ind[1] for ind in new_list]

    plt.imshow(probs_df)
    plt.scatter(xs, ys, s=2.5, c="red")
    plt.show()
