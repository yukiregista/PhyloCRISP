# <img width="35" alt="logo" src="https://github.com/user-attachments/assets/b05ea043-676b-406d-a2be-8cc45c9cc315"/> PhyloCRISP: Phylogenetic Consensus Resolution via Improved Split Proximity

PhyloCRISP is a software package for constructing more resolved, meaningful consensus trees from a set of phylogenetic trees than the conventional majority-rule trees. It implements methods proposed in the paper

> Yuki Takazawa, Atsushi Takeda, Momoko Hayamizu, and Olivier Gascuel. **Outperforming the Majority-Rule Consensus Tree Using Fine-Grained Dissimilarity Measures**, _submitted_.

The software supports both a graphical user interface (GUI) and a command-line interface (CLI), as demonstrated in the following tutorial.

<!-- While the entire datasets analyzed in the paper is available at DRYAD, some of the data are also found in `example_data` for the tutorial purpose.  -->

Some of the data can be found in `example_data` for the tutorial purpose.

## Install (GUI)

You can download and install the GUI version of PhyloCRISP for Windows, macOS, and Linux:

- [Windows](https://github.com/yukiregista/PhyloCRISP/releases/latest/download/PhyloCRISP-Setup-0.0.1-win-x64.exe)
- [Mac](https://github.com/yukiregista/PhyloCRISP/releases/latest/download/PhyloCRISP-0.0.1-macos.dmg)
- [Linux](https://github.com/yukiregista/PhyloCRISP/releases/latest/download/PhyloCRISP-0.0.1-linux-x64.AppImage)

## Tutorial (GUI)

This tutorial demonstrates how to reproduce the consensus tree from the **mammals dataset** (Figure 6 in the paper). We also illustrate how to compare it with the majority rule consensus tree.

### Step 1: Prepare Input Files

You will need the following two files from the Mammals dataset:

- [bootstrap_trees_full.standardized.nw](https://github.com/yukiregista/PhyloCRISP/releases/latest/download/bootstrap_trees_full.standardized.nw) — the set of bootstrap trees (input trees, ~30 MB)
- [reference_tree.standardized.nw](https://github.com/yukiregista/PhyloCRISP/releases/latest/download/reference_tree.standardized.nw) — the reference tree used as the starting tree (~90 KB)

<!-- These files are also available at DRYAD. -->

### Step 2: Launch the Application

Run the GUI application.

On macOS and Windows, you may see security warnings when running software from outside the official app store or package manager. On macOS, allow the app to run via **System Settings → Privacy & Security** ([Apple support](https://support.apple.com/guide/mac-help/open-a-mac-app-from-an-unknown-developer-mh40616/mac)). On Windows, click **More info → Run anyway** in the SmartScreen warning ([Microsoft SmartScreen overview](https://learn.microsoft.com/en-us/windows/security/operating-system-security/virus-and-threat-protection/microsoft-defender-smartscreen)).

### Step 3: Configure the Settings

In the GUI, configure the following options as shown in the screenshot below:

1. **Input Trees**: Click _Browse Input Trees File_ and select `bootstrap_trees_full.standardized.nw`.
2. **Consensus Method**: Here we select **Unscaled Transfer** from the dropdown menu to measure the dissimilarity between trees using this type of distance.
3. **Optimization Strategy**: Here we select **Prune Only** for computing a median tree, which requires us specify an initial tree in the next step.
4. **Starting Tree**: Click _Browse Starting Tree File_ and select `reference_tree.standardized.nw` as an initial tree. _(This option is required when using Prune Only.)_
5. **Output Filename**: Set to a name of your choice. We also recommend tick the checkbox "Auto-save output".
6. **Output Directory**: Choose your desired output directory.

<img width="605" height="1064" alt="image" src="https://github.com/user-attachments/assets/115e8d4e-5939-464f-b9c8-1185770c2253" />

### Step 4: Run and Save

Click the **Generate Consensus Tree** button to start the computation. For this dataset, the computation is expected to end within 10 minutes. Once complete, the consensus tree will be displayed in the _Consensus Tree_ panel at the bottom, and the Newick file (.nwk), together with the transfer support values (.nwk.tsupp) and the frequency values (.nwk.fsupp), will be saved in your specified directly. (If you forgot to tick "Auto-save output", you can still save the result to your specified output directory by clicking the **Save Output** button.

### Step 5: Visualization of the Consensus Tree

You can visualize the consensus tree using your preferred phylogenetic tool. You can easily reproduce the image in Figure 6(b) of the paper by uploading the newick file to [ITOL v7](https://itol.embl.de/).

<img width="1302" height="957" alt="image" src="https://github.com/user-attachments/assets/36a37c7c-7f68-4b03-9fee-730f16016786" />

### Step 6: Comparison against the Majority Rule Consensus Tree

You can easily compare this consensus tree with the majority rule consensus tree. To get the majority rule tree for the same dataset, select **Majority Rule** from the dropdown menu for Consensus Methods in Step 3. You can reproduce Figure 6(a) of the paper by visualizing the output Newick file with [ITOL v7](https://itol.embl.de/).

## Install (CLI)

You can install and use the CLI version of PhyloCRISP as follows.
CLI support via `pip install` is currently Linux/macOS only.

1. Clone the repository:

```sh
git clone git@github.com:yukiregista/PhyloCRISP.git
```

2. Install the package:

```sh
pip install .
```

3. Build `booster` (for support-value outputs and transfer prune-only workflows):

```sh
cd src/booster/src
make
cd ../../..
```

On macOS with Apple clang, install OpenMP first if needed:

```sh
brew install libomp
```

4. Ensure `tqdist` is available when using `--method quartet`:

- Install `tqdist` by following the official instructions: https://www.birc.au.dk/~cstorm/software/tqdist/
- Ensure `quartet_dist` and `pairs_quartet_dist` are available on `PATH`

5. Confirm installation:

```sh
phylocrisp-cli --help
```

## Usage (CLI)

```sh
phylocrisp-cli \
  --input-trees <path/to/input_trees.nw> \
  --method <scaled_transfer|unscaled_transfer|quartet|mr> \
  [--strategy <add_and_prune|prune>] \
  [--starting-tree <path/to/starting_tree.nw>] \
  [--output-filename <name_or_path.nwk>] \
  [--output-dir <output_directory>]
```

File format information:

- `--input-trees`: Newick file containing input trees.
- `--starting-tree`: Newick file containing one starting tree.
- Main output is written as Newick (`*.nwk`).
- When `booster` is available, additional outputs are `*.fsupp` (frequency support) and `*.tsupp` (transfer support).

Argument notes:

- `--input-trees` and `--method` are required.
- `--strategy` defaults to `add_and_prune`.
- `--method quartet` supports `--strategy prune` only.
- `--starting-tree` is required for prune workflows.
