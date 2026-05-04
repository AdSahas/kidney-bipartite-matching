import pandas as pd
import matplotlib.pyplot as plt
from statsmodels.nonparametric.smoothers_lowess import lowess


def plot_metric_lowess(metric, ylabel, output_plot, output_csv, frac=0.2):
    greedy_df = pd.read_csv("greedy.csv")
    kuhn_df   = pd.read_csv("kuhn.csv")
    jv_df     = pd.read_csv("jv.csv")

    plt.figure(figsize=(10, 6))

    def process(df, label, marker):
        x = df["size"].values
        y = df[metric].values

        # LOWESS smoothing
        smoothed = lowess(y, x, frac=frac, return_sorted=False)

        # Plot raw (faint)
        plt.plot(x, y, marker=marker, alpha=0.25)

        # Plot smoothed (bold)
        plt.plot(x, smoothed, label=label)

        return smoothed

    greedy_smooth = process(greedy_df, "Greedy", 'o')
    kuhn_smooth   = process(kuhn_df, "Kuhn", 's')
    jv_smooth     = process(jv_df, "Jonker-Volgenant", '^')

    plt.xlabel("Graph Size")
    plt.ylabel(ylabel)
    plt.title(f"{ylabel} vs Graph Size (LOWESS Smoothed)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(output_plot)
    plt.close()

    print(f"Saved {output_plot}")

    # Save smoothed data to CSV
    smoothed_df = pd.DataFrame({
        "size": greedy_df["size"],
        "greedy_smoothed": greedy_smooth,
        "kuhn_smoothed": kuhn_smooth,
        "jv_smoothed": jv_smooth
    })

    smoothed_df.to_csv(output_csv, index=False)
    print(f"Saved {output_csv}")


def graphAllLowess():
    plot_metric_lowess("time", "Execution Time (seconds)", "time_lowess.png", "time_lowess.csv")
    plot_metric_lowess("matching", "Average Matching Size", "matching_lowess.png", "matching_lowess.csv")
    plot_metric_lowess("hla", "Average HLA Score", "hla_lowess.png", "hla_lowess.csv")
    plot_metric_lowess("pra", "Average PRA Score", "pra_lowess.png", "pra_lowess.csv")


if __name__ == "__main__":
    graphAllLowess()