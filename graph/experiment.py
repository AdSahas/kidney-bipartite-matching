# Driver for running experiments. This code was
# Modified by Claude Sonnet 4.6 to enhance ease of use.


from compatibility import *
from builder import *
import timeit
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ─────────────────────────────────────────────
# Experiment 1: Scalability (varying graph size)
# ─────────────────────────────────────────────

def runAll(rej_chance):
    sizes = range(100, 510, 5)

    greedyTime, greedyMatching, greedyHLA, greedyPRA = [], [], [], []
    kuhnTime, kuhnMatching, kuhnHLA, kuhnPRA = [], [], [], []
    jvTime, jvMatching, jvHLA, jvPRA = [], [], [], []

    for size in sizes:
        print(f"Running all algorithms on size {size}")

        gTimes, gMatch, gHLA, gPRA = [], [], [], []
        kTimes, kMatch, kHLA, kPRA = [], [], [], []
        jTimes, jMatch, jHLA, jPRA = [], [], [], []

        for _ in range(3):
            donors, recipients = build_graph(size, size)

            # Greedy
            start = timeit.default_timer()
            gMatches = GreedyMatching(donors, recipients, rej_chance)
            stop = timeit.default_timer()

            gTimes.append(stop - start)
            gPairs = [(m, i) for i, m in enumerate(gMatches) if m != -1]
            gMatch.append(len(gPairs))

            if gPairs:
                gHLA.append(np.mean([HLAMatchScore(donors[d], recipients[r]) for d, r in gPairs]))
                gPRA.append(np.mean([recipients[r]['pra'] for d, r in gPairs]))
            else:
                gHLA.append(0)
                gPRA.append(0)

            # Kuhn
            start = timeit.default_timer()
            kMatches = kuhn(donors, recipients, rej_chance)
            stop = timeit.default_timer()

            kTimes.append(stop - start)
            kPairs = [(m, i) for i, m in enumerate(kMatches) if m != -1]
            kMatch.append(len(kPairs))

            if kPairs:
                kHLA.append(np.mean([HLAMatchScore(donors[d], recipients[r]) for d, r in kPairs]))
                kPRA.append(np.mean([recipients[r]['pra'] for d, r in kPairs]))
            else:
                kHLA.append(0)
                kPRA.append(0)

            # Jonker-Volgenant
            start = timeit.default_timer()
            jMatches = jonkerVolgenant(donors, recipients, rej_chance)
            stop = timeit.default_timer()

            jTimes.append(stop - start)
            jMatch.append(len(jMatches))

            if jMatches:
                jHLA.append(np.mean([HLAMatchScore(donors[d], recipients[r]) for d, r in jMatches]))
                jPRA.append(np.mean([recipients[r]['pra'] for d, r in jMatches]))
            else:
                jHLA.append(0)
                jPRA.append(0)

        greedyTime.append(np.mean(gTimes))
        greedyMatching.append(np.mean(gMatch))
        greedyHLA.append(np.mean(gHLA))
        greedyPRA.append(np.mean(gPRA))

        kuhnTime.append(np.mean(kTimes))
        kuhnMatching.append(np.mean(kMatch))
        kuhnHLA.append(np.mean(kHLA))
        kuhnPRA.append(np.mean(kPRA))

        jvTime.append(np.mean(jTimes))
        jvMatching.append(np.mean(jMatch))
        jvHLA.append(np.mean(jHLA))
        jvPRA.append(np.mean(jPRA))

    pd.DataFrame({
        'size': list(sizes),
        'time': greedyTime,
        'matching': greedyMatching,
        'hla': greedyHLA,
        'pra': greedyPRA
    }).to_csv('greedy.csv', index=False)

    pd.DataFrame({
        'size': list(sizes),
        'time': kuhnTime,
        'matching': kuhnMatching,
        'hla': kuhnHLA,
        'pra': kuhnPRA
    }).to_csv('kuhn.csv', index=False)

    pd.DataFrame({
        'size': list(sizes),
        'time': jvTime,
        'matching': jvMatching,
        'hla': jvHLA,
        'pra': jvPRA
    }).to_csv('jv.csv', index=False)

    print("Saved greedy.csv, kuhn.csv, jv.csv")


def plot_metric(metric, ylabel, output_file):
    greedy_df = pd.read_csv("greedy.csv")
    kuhn_df = pd.read_csv("kuhn.csv")
    jv_df = pd.read_csv("jv.csv")

    plt.figure(figsize=(10, 6))
    plt.plot(greedy_df["size"], greedy_df[metric], marker='o', label="Greedy")
    plt.plot(kuhn_df["size"], kuhn_df[metric], marker='s', label="Kuhn")
    plt.plot(jv_df["size"], jv_df[metric], marker='^', label="Jonker-Volgenant")

    plt.xlabel("Graph Size")
    plt.ylabel(ylabel)
    plt.title(f"{ylabel} vs Graph Size")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

    print(f"Saved {output_file}")


def graphAll():
    plot_metric("time", "Execution Time (seconds)", "time_plot.png")
    plot_metric("matching", "Average Matching Size", "matching_plot.png")
    plot_metric("hla", "Average HLA Score", "hla_plot.png")
    plot_metric("pra", "Average PRA Score", "pra_plot.png")


# ─────────────────────────────────────────────
# Experiment 2: Donor : Recipient Ratio
# ─────────────────────────────────────────────

def runRatios(rej_chance=0.2, n_recipients=500, trials=5):
    ratios = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    rows = []

    for ratio in ratios:
        n_donors = int(n_recipients * ratio)
        print(f"Running ratio {ratio:.1f} with donors={n_donors}, recipients={n_recipients}")

        gTimes, gMatch, gHLA, gPRA = [], [], [], []
        kTimes, kMatch, kHLA, kPRA = [], [], [], []
        jTimes, jMatch, jHLA, jPRA = [], [], [], []

        for _ in range(trials):
            donors, recipients = build_graph(n_donors, n_recipients)

            # Greedy
            start = timeit.default_timer()
            gMatches = GreedyMatching(donors, recipients, rej_chance)
            stop = timeit.default_timer()

            gTimes.append(stop - start)
            gPairs = [(m, i) for i, m in enumerate(gMatches) if m != -1]
            gMatch.append(len(gPairs))

            if gPairs:
                gHLA.append(np.mean([HLAMatchScore(donors[d], recipients[r]) for d, r in gPairs]))
                gPRA.append(np.mean([recipients[r]['pra'] for d, r in gPairs]))
            else:
                gHLA.append(0)
                gPRA.append(0)

            # Kuhn
            start = timeit.default_timer()
            kMatches = kuhn(donors, recipients, rej_chance)
            stop = timeit.default_timer()

            kTimes.append(stop - start)
            kPairs = [(m, i) for i, m in enumerate(kMatches) if m != -1]
            kMatch.append(len(kPairs))

            if kPairs:
                kHLA.append(np.mean([HLAMatchScore(donors[d], recipients[r]) for d, r in kPairs]))
                kPRA.append(np.mean([recipients[r]['pra'] for d, r in kPairs]))
            else:
                kHLA.append(0)
                kPRA.append(0)

            # Jonker-Volgenant
            start = timeit.default_timer()
            jMatches = jonkerVolgenant(donors, recipients, rej_chance)
            stop = timeit.default_timer()

            jTimes.append(stop - start)
            jMatch.append(len(jMatches))

            if jMatches:
                jHLA.append(np.mean([HLAMatchScore(donors[d], recipients[r]) for d, r in jMatches]))
                jPRA.append(np.mean([recipients[r]['pra'] for d, r in jMatches]))
            else:
                jHLA.append(0)
                jPRA.append(0)

        rows.append({
            'algorithm': 'Greedy',
            'ratio': ratio,
            'donors': n_donors,
            'recipients': n_recipients,
            'time': np.mean(gTimes),
            'matching': np.mean(gMatch),
            'match_pct': np.mean(gMatch) / n_recipients,
            'hla': np.mean(gHLA),
            'pra': np.mean(gPRA)
        })

        rows.append({
            'algorithm': 'Kuhn',
            'ratio': ratio,
            'donors': n_donors,
            'recipients': n_recipients,
            'time': np.mean(kTimes),
            'matching': np.mean(kMatch),
            'match_pct': np.mean(kMatch) / n_recipients,
            'hla': np.mean(kHLA),
            'pra': np.mean(kPRA)
        })

        rows.append({
            'algorithm': 'Jonker-Volgenant',
            'ratio': ratio,
            'donors': n_donors,
            'recipients': n_recipients,
            'time': np.mean(jTimes),
            'matching': np.mean(jMatch),
            'match_pct': np.mean(jMatch) / n_recipients,
            'hla': np.mean(jHLA),
            'pra': np.mean(jPRA)
        })

    pd.DataFrame(rows).to_csv("ratios.csv", index=False)
    print("Saved ratios.csv")


def plotRatios(metric, ylabel, output_file):
    df = pd.read_csv("ratios.csv")

    plt.figure(figsize=(10, 6))

    for algo, marker in [('Greedy', 'o'), ('Kuhn', 's'), ('Jonker-Volgenant', '^')]:
        temp = df[df['algorithm'] == algo]
        plt.plot(temp['ratio'], temp[metric], marker=marker, label=algo)

    plt.xlabel("Donor : Recipient Ratio")
    plt.ylabel(ylabel)
    plt.title(f"{ylabel} vs Donor : Recipient Ratio")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

    print(f"Saved {output_file}")


def graphRatioResults():
    plotRatios("time", "Execution Time (seconds)", "ratio_time.png")
    plotRatios("matching", "Average Matching Size", "ratio_matching.png")
    plotRatios("match_pct", "Match Percentage", "ratio_match_pct.png")
    plotRatios("hla", "Average HLA Score", "ratio_hla.png")
    plotRatios("pra", "Average PRA Score", "ratio_pra.png")


# ─────────────────────────────────────────────
# Experiment 3: Varying Rejection Threshold
# ─────────────────────────────────────────────

def runRejectionChance(size=500, trials=5):
    rej_chances = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5]
    rows = []

    for rej_chance in rej_chances:
        print(f"Running rejection chance {rej_chance} with size {size}")

        gTimes, gMatch, gHLA, gPRA = [], [], [], []
        kTimes, kMatch, kHLA, kPRA = [], [], [], []
        jTimes, jMatch, jHLA, jPRA = [], [], [], []

        for _ in range(trials):
            donors, recipients = build_graph(size, size)

            # Greedy
            start = timeit.default_timer()
            gMatches = GreedyMatching(donors, recipients, rej_chance)
            stop = timeit.default_timer()

            gTimes.append(stop - start)
            gPairs = [(m, i) for i, m in enumerate(gMatches) if m != -1]
            gMatch.append(len(gPairs))

            if gPairs:
                gHLA.append(np.mean([HLAMatchScore(donors[d], recipients[r]) for d, r in gPairs]))
                gPRA.append(np.mean([recipients[r]['pra'] for d, r in gPairs]))
            else:
                gHLA.append(0)
                gPRA.append(0)

            # Kuhn
            start = timeit.default_timer()
            kMatches = kuhn(donors, recipients, rej_chance)
            stop = timeit.default_timer()

            kTimes.append(stop - start)
            kPairs = [(m, i) for i, m in enumerate(kMatches) if m != -1]
            kMatch.append(len(kPairs))

            if kPairs:
                kHLA.append(np.mean([HLAMatchScore(donors[d], recipients[r]) for d, r in kPairs]))
                kPRA.append(np.mean([recipients[r]['pra'] for d, r in kPairs]))
            else:
                kHLA.append(0)
                kPRA.append(0)

            # Jonker-Volgenant
            start = timeit.default_timer()
            jMatches = jonkerVolgenant(donors, recipients, rej_chance)
            stop = timeit.default_timer()

            jTimes.append(stop - start)
            jMatch.append(len(jMatches))

            if jMatches:
                jHLA.append(np.mean([HLAMatchScore(donors[d], recipients[r]) for d, r in jMatches]))
                jPRA.append(np.mean([recipients[r]['pra'] for d, r in jMatches]))
            else:
                jHLA.append(0)
                jPRA.append(0)

        rows.append({
            'algorithm': 'Greedy',
            'rej_chance': rej_chance,
            'size': size,
            'time': np.mean(gTimes),
            'matching': np.mean(gMatch),
            'match_pct': np.mean(gMatch) / size,
            'hla': np.mean(gHLA),
            'pra': np.mean(gPRA)
        })

        rows.append({
            'algorithm': 'Kuhn',
            'rej_chance': rej_chance,
            'size': size,
            'time': np.mean(kTimes),
            'matching': np.mean(kMatch),
            'match_pct': np.mean(kMatch) / size,
            'hla': np.mean(kHLA),
            'pra': np.mean(kPRA)
        })

        rows.append({
            'algorithm': 'Jonker-Volgenant',
            'rej_chance': rej_chance,
            'size': size,
            'time': np.mean(jTimes),
            'matching': np.mean(jMatch),
            'match_pct': np.mean(jMatch) / size,
            'hla': np.mean(jHLA),
            'pra': np.mean(jPRA)
        })

    pd.DataFrame(rows).to_csv("rejection.csv", index=False)
    print("Saved rejection.csv")


def plotRejection(metric, ylabel, output_file):
    df = pd.read_csv("rejection.csv")

    plt.figure(figsize=(10, 6))

    for algo, marker in [('Greedy', 'o'), ('Kuhn', 's'), ('Jonker-Volgenant', '^')]:
        temp = df[df['algorithm'] == algo]
        plt.plot(temp['rej_chance'], temp[metric], marker=marker, label=algo)

    plt.xlabel("Rejection Chance Threshold")
    plt.ylabel(ylabel)
    plt.title(f"{ylabel} vs Rejection Chance Threshold")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

    print(f"Saved {output_file}")


def graphRejectionResults():
    plotRejection("time", "Execution Time (seconds)", "rejection_time.png")
    plotRejection("matching", "Average Matching Size", "rejection_matching.png")
    plotRejection("match_pct", "Match Percentage", "rejection_match_pct.png")
    plotRejection("hla", "Average HLA Score", "rejection_hla.png")
    plotRejection("pra", "Average PRA Score", "rejection_pra.png")


# ─────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────

if __name__ == "__main__":
    # Experiment 1: Scalability
    print("\n=== Experiment 1: Scalability ===")
    runAll(rej_chance=0.2)
    graphAll()

    # Experiment 2: Donor : Recipient Ratio
    print("\n=== Experiment 2: Donor : Recipient Ratio ===")
    runRatios(rej_chance=0.2, n_recipients=500, trials=5)
    graphRatioResults()

    # Experiment 3: Rejection Threshold
    print("\n=== Experiment 3: Rejection Threshold ===")
    runRejectionChance(size=500, trials=5)
    graphRejectionResults()

    print("\nDone. Saved all CSVs and plots.")