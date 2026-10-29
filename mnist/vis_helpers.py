## Visualizations for grid plots
import pandas as pd
import numpy as np
import plotnine as p9
# from plotnine import ggplot, aes, geom_tile, facet_wrap, scale_fill_cmap, theme_void, theme, element_text

def build_tidy_grid(images_array, cluster_label, padding=2, grid_width=3):
    data_list = []
    H = images_array.shape[1]
    W = images_array.shape[2]
    for idx in range(len(images_array)):
        # Calculate 3x3 grid position (0 to 2)
        grid_row = idx // grid_width
        grid_col = idx % grid_width
        
        # Calculate coordinate offsets for this image slot
        x_offset = grid_col * (W + padding)
        y_offset = (2 - grid_row) * (H + padding) # 2 - grid_row keeps row 0 at the top
        
        for r in range(H):
            for c in range(W):
                # Map pixels to absolute canvas coordinates
                # (H - 1 - r) flips the Y-axis so images aren't upside down
                data_list.append({
                    'X': x_offset + c,
                    'Y': y_offset + (H - 1 - r),
                    'Intensity': images_array[idx, r, c],
                    'Cluster': cluster_label
                })
                
    return pd.DataFrame(data_list)


def produce_grid_plot(X, masks, titles, grid_width=7, seed=1, ncol=None, nrow=1):
    if ncol is None:
        ncol = len(masks)

    plt_dfs = []
    np.random.seed(seed)
    for idx in range(len(masks)):
        mask = masks[idx]
        title = titles[idx]
        sampled = np.random.choice([x for x in range(np.sum(mask))], size=grid_width**2, replace=False)
        imgs = X[mask][sampled, 0]
        plt_dfs.append(build_tidy_grid(imgs, title, grid_width=grid_width))

    plot_df = pd.concat(plt_dfs, ignore_index=True)
    plot_df["Cluster"] = pd.Categorical(plot_df["Cluster"], 
                                       categories=titles, ordered=True)
    return (
        p9.ggplot(plot_df, p9.aes(x='X', y='Y', fill='Intensity'))
        + p9.geom_tile()
        + p9.facet_wrap('~Cluster', ncol=ncol, nrow=nrow)
        + p9.scale_fill_cmap('gray')
        + p9.theme_void()
        + p9.theme(
            strip_text=p9.element_text(size=14, face='bold'),
            legend_position="none"
        )
    )

def map_clusters_to_label_mask(posterior_matrix, class_labels, target_class, top_n=2):
    posterior_matrix = np.asarray(posterior_matrix)
    class_labels = np.asarray(class_labels)
    class_mask = (class_labels == target_class)
    mean_posteriors = np.mean(posterior_matrix[class_mask], axis=0)
    indices = np.argsort(mean_posteriors)[-top_n:][::-1]
    output = []
    for idx in indices:
        output.append((np.argmax(posterior_matrix, axis=1) == idx))

    return output