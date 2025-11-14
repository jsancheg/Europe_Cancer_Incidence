import matplotlib.pyplot as plt
def plot_serie(ts_data, title, xlabel, ylabel):
    plt.figure(figsize=(10, 4))
    ts_data.plot(title=title, linewidth=2, marker='o')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.show()



def descriptive_statistics(ts_data):
    ds_dict = {'Mean': [ts_data.mean()],
               'Median': [ts_data.median()],
               'SD': [ts_data.std()],
               'Min': [ts_data.min()],
               'Max': [ts_data.max()],
               'Range': [ts_data.max() - ts_data.min()],
               'CV': [(ts_data.std() / ts_data.mean() ) *100]
    }
    print ("==== Descriptive Statistics ===")

    return ds_dict

