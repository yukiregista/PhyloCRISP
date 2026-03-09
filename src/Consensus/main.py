import argparse
import os
import platform
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from ._consensus import Tree_with_support, TreeList_with_support
from ._greedy import STDGreedyConsensus, SUTDGreedyConsensus
from ._c_pruning import c_prune
from ._app_options import (
    effective_strategy_token,
    normalize_method_token,
    resolve_output_path,
)


class ConsensusError(Exception):
    """Custom exception for errors in the consensus process."""
    pass

def _strip_branch_lengths(tree):
    """Remove branch lengths from all edges of the tree."""
    for edge in tree.postorder_edge_iter():
        edge.length = None
    return tree


def _clear_internal_labels(tree):
    """Remove internal node labels (used for plain output tree)."""
    for node in tree.postorder_internal_node_iter(exclude_seed_node=True):
        node.label = None
    return tree


def _write_plain_newick(tree, output_path):
    out_tree = tree.clone(depth=1)
    out_tree.encode_bipartitions()
    _strip_branch_lengths(out_tree)
    _clear_internal_labels(out_tree)
    out_tree.write(path=output_path, schema="newick", suppress_rooting=True)


def _plain_newick_string(tree):
    out_tree = tree.clone(depth=1)
    out_tree.encode_bipartitions()
    _strip_branch_lengths(out_tree)
    _clear_internal_labels(out_tree)
    return out_tree.as_string(schema="newick", suppress_rooting=True).strip()


def _write_support_newick(tree, support_map, output_path):
    out_tree = tree.clone(depth=1)
    out_tree.encode_bipartitions()
    _strip_branch_lengths(out_tree)
    for node in out_tree.postorder_internal_node_iter(exclude_seed_node=True):
        splitint = node.bipartition.split_as_int()
        node.label = f"{support_map.get(splitint, 0.0):.6f}"
    out_tree.write(path=output_path, schema="newick", suppress_rooting=True)


def _platform_folder():
    if sys.platform == "win32":
        return "win"
    if sys.platform == "darwin":
        return "mac"
    return "linux"


def _arch_folder():
    machine = platform.machine().lower()
    if machine in ("x86_64", "amd64"):
        return "x64"
    if machine in ("arm64", "aarch64"):
        return "arm64"
    if machine in ("i386", "i686", "x86"):
        return "x86"
    return machine or "unknown"


def _iter_booster_exec_candidates():
    seen = set()

    def _yield(path_like):
        if not path_like:
            return
        path = os.path.abspath(os.fspath(path_like))
        norm = os.path.normcase(path)
        if norm in seen:
            return
        seen.add(norm)
        if os.path.isfile(path):
            yield path

    exe_name = "booster.exe" if os.name == "nt" else "booster"

    env_path = os.environ.get("CONSENSUS_BOOSTER_EXECUTABLE")
    if env_path:
        yield from _yield(env_path)

    which_path = shutil.which("booster")
    if which_path:
        yield from _yield(which_path)

    module_dir = Path(__file__).resolve().parent
    repo_src_dir = module_dir.parent
    platform_dir = Path(_platform_folder())
    arch_dir = Path(_arch_folder())
    universal_dir = platform_dir / "universal2"

    candidates = [
        module_dir / "bin" / platform_dir / arch_dir / exe_name,
        module_dir / "bin" / universal_dir / exe_name,
        module_dir.parent / "bin" / platform_dir / arch_dir / exe_name,
        module_dir.parent / "bin" / universal_dir / exe_name,
        repo_src_dir / "booster" / "src" / exe_name,
    ]

    if hasattr(sys, "_MEIPASS"):
        meipass = Path(sys._MEIPASS)
        candidates.extend(
            [
                meipass / exe_name,
                meipass / "_internal" / exe_name,
                meipass / platform_dir / arch_dir / exe_name,
                meipass / "_internal" / platform_dir / arch_dir / exe_name,
                meipass / universal_dir / exe_name,
                meipass / "_internal" / universal_dir / exe_name,
            ]
        )

    exe_dir = Path(sys.executable).resolve().parent if sys.executable else None
    if exe_dir:
        candidates.extend(
            [
                exe_dir / exe_name,
                exe_dir / "_internal" / exe_name,
                exe_dir / "bin" / platform_dir / arch_dir / exe_name,
                exe_dir / "_internal" / "bin" / platform_dir / arch_dir / exe_name,
                exe_dir / "bin" / universal_dir / exe_name,
                exe_dir / "_internal" / "bin" / universal_dir / exe_name,
            ]
        )

    for candidate in candidates:
        yield from _yield(candidate)


def _find_booster_executable():
    candidates = list(_iter_booster_exec_candidates())
    if candidates:
        return candidates[0]
    raise ConsensusError(
        "Cannot find the compiled 'booster' executable. "
        "Install/build booster and make it available on PATH, or set CONSENSUS_BOOSTER_EXECUTABLE."
    )


def _run_booster_cli(booster_executable, algo, reference_tree_file, bootstrap_trees_file, output_file):
    cmd = [
        booster_executable,
        "-a",
        algo,
        "-i",
        reference_tree_file,
        "-b",
        bootstrap_trees_file,
        "-o",
        output_file,
        "-q",
    ]
    run_kwargs = {"capture_output": True, "text": True}
    if os.name == "nt":
        # Hide console windows when invoking booster.exe from the GUI build.
        creationflags = getattr(subprocess, "CREATE_NO_WINDOW", 0)
        startupinfo = subprocess.STARTUPINFO()
        startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW
        run_kwargs["creationflags"] = creationflags
        run_kwargs["startupinfo"] = startupinfo

    proc = subprocess.run(cmd, **run_kwargs)
    if proc.returncode != 0:
        raise ConsensusError(
            f"booster failed with algorithm '{algo}' (exit code {proc.returncode}).\n"
            f"Command: {' '.join(cmd)}\n"
            f"stderr: {proc.stderr.strip()}\n"
            f"stdout: {proc.stdout.strip()}"
        )


def _read_support_tree_to_map(support_tree_file, taxon_namespace):
    support_tree = Tree_with_support.get(
        path=support_tree_file,
        schema="newick",
        taxon_namespace=taxon_namespace,
    )
    support_tree.encode_bipartitions()
    support_map = {}
    for node in support_tree.postorder_internal_node_iter(exclude_seed_node=True):
        splitint = node.bipartition.split_as_int()
        try:
            support_value = float(node.label) if node.label is not None else 0.0
        except ValueError:
            support_value = 0.0
        support_map[splitint] = support_value
    return support_map


def _compute_support_maps_with_booster(consensus_tree, input_trees):
    support_maps = {}
    errors = []
    try:
        booster_executable = _find_booster_executable()
    except ConsensusError as e:
        return support_maps, [str(e)]

    with tempfile.TemporaryDirectory(prefix="consensus_booster_") as tmpdir:
        reference_tree_file = os.path.join(tmpdir, "reference.nwk")
        bootstrap_trees_file = os.path.join(tmpdir, "bootstrap_trees.nwk")
        fbp_tree_file = os.path.join(tmpdir, "fbp_output.nwk")
        tbe_tree_file = os.path.join(tmpdir, "tbe_output.nwk")

        with open(reference_tree_file, "w") as f:
            f.write(_plain_newick_string(consensus_tree) + "\n")
        with open(bootstrap_trees_file, "w") as f:
            for tree in input_trees:
                f.write(_plain_newick_string(tree) + "\n")

        try:
            _run_booster_cli(booster_executable, "fbp", reference_tree_file, bootstrap_trees_file, fbp_tree_file)
            support_maps["fsupp"] = _read_support_tree_to_map(fbp_tree_file, consensus_tree.taxon_namespace)
        except Exception as e:
            errors.append(f"booster fbp failed: {e}")

        try:
            _run_booster_cli(booster_executable, "rtbe", reference_tree_file, bootstrap_trees_file, tbe_tree_file)
            support_maps["tsupp"] = _read_support_tree_to_map(tbe_tree_file, consensus_tree.taxon_namespace)
        except Exception as e:
            errors.append(f"booster rtbe failed: {e}")

    return support_maps, errors


def _write_all_outputs(consensus_tree, input_trees, output_file):
    _write_plain_newick(consensus_tree, output_file)
    support_maps, errors = _compute_support_maps_with_booster(consensus_tree, input_trees)
    support_outputs = []

    if "fsupp" in support_maps:
        fsupp_path = output_file + ".fsupp"
        _write_support_newick(consensus_tree, support_maps["fsupp"], fsupp_path)
        support_outputs.append((fsupp_path, "frequency support"))
    if "tsupp" in support_maps:
        tsupp_path = output_file + ".tsupp"
        _write_support_newick(consensus_tree, support_maps["tsupp"], tsupp_path)
        support_outputs.append((tsupp_path, "transfer support"))

    for err in errors:
        print(f"Warning: {err}", file=sys.stderr)

    return output_file, support_outputs


def main(args=None):
    parser = argparse.ArgumentParser(
        prog="PhyloCRISP CLI",
        description=(
            "Build a consensus tree from input trees and write Newick output.\n"
            "The CLI writes one main tree file (*.nwk) and two support-value trees:\n"
            "*.fsupp (frequency support; fbp in the case of bootstrap trees) and *.tsupp (transfer support; tbe in the case of bootstrap trees)."
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--input-trees",
        dest="input_trees",
        type=str,
        help="Path to input tree file in Newick format.",
    )
    parser.add_argument(
        "--method",
        type=str,
        help=(
            "Consensus method:\n"
            "  scaled_transfer   Scaled transfer distance greedy consensus\n"
            "  unscaled_transfer Unscaled transfer distance greedy consensus\n"
            "  quartet           Quartet-distance prune-only consensus\n"
            "  mr                Majority-rule consensus"
        ),
    )
    parser.add_argument(
        "--strategy",
        type=str,
        default="add_and_prune",
        help=(
            "Optimization strategy:\n"
            "  add_and_prune  Build from initial/majority tree (default)\n"
            "  prune          Prune-only optimization (requires --starting-tree except for quartet)"
        ),
    )
    parser.add_argument(
        "--starting-tree",
        dest="starting_tree",
        type=str,
        help="Path to a starting tree in Newick format.",
    )
    parser.add_argument(
        "--output-filename",
        dest="output_filename",
        type=str,
        help=(
            "Main output filename/path (for the .nwk tree).\n"
            "If omitted, a default name is generated from method/strategy."
        ),
    )
    parser.add_argument(
        "--output-dir",
        dest="output_dir",
        type=str,
        help="Output directory (used with --output-filename or default naming).",
    )

    if args is None:
        args = sys.argv[1:]

    if not args:
        parser.print_help()
        return

    args = parser.parse_args(args)

    if args.input_trees is None:
        raise ConsensusError("Please specify the path to input trees with --input-trees.")

    if args.method is None:
        raise ConsensusError("Please specify the consensus method with --method.")
    try:
        method_token = normalize_method_token(args.method)
    except ValueError as e:
        raise ConsensusError(str(e)) from e

    try:
        strategy_token = effective_strategy_token(method_token, args.strategy)
    except ValueError as e:
        raise ConsensusError(str(e)) from e

    output_file = resolve_output_path(
        input_trees_path=args.input_trees,
        method_token=method_token,
        strategy_token=strategy_token,
        output_filename=args.output_filename,
        output_dir=args.output_dir,
        starting_tree_path=args.starting_tree,
    )

    dir_name = os.path.dirname(output_file)
    if dir_name and not os.path.exists(dir_name):
        raise ConsensusError(f"Directory '{dir_name}' does not exist.")

    if method_token == "quartet" and strategy_token != "prune":
        raise ConsensusError("Method 'quartet' currently supports prune-only mode only.")

    try:
        input_trees = TreeList_with_support.get(path=args.input_trees, schema="newick")
    except Exception as e:
        raise ConsensusError(f"Error in reading input trees (Newick expected): {e}")

    consensus_tree = None

    if method_token == "mr":
        consensus_tree = input_trees.majority_rule_consensus()
    elif strategy_token == "prune" and method_token in ["scaled_transfer", "unscaled_transfer"]:
        if args.starting_tree is None:
            raise ConsensusError("You need to specify the path to the starting tree with --starting-tree.")
        try:
            starting_tree = Tree_with_support.get(
                path=args.starting_tree,
                schema="newick",
                taxon_namespace=input_trees.taxon_namespace,
            )
        except Exception as e:
            raise ConsensusError(f"Error in reading starting tree: {e}")

        n_taxa = len(starting_tree.taxon_namespace)
        if n_taxa < 4:
            raise ConsensusError("Starting tree must have at least four taxa for pruning.")

        if method_token == "scaled_transfer":
            gre_nwk = c_prune(args.starting_tree, args.input_trees, K=min(30, n_taxa - 3), scaled=True)
        else:
            gre_nwk = c_prune(args.starting_tree, args.input_trees, K=min(30, n_taxa - 3), scaled=False)

        consensus_tree = Tree_with_support.get(
            data=gre_nwk,
            schema="newick",
            taxon_namespace=input_trees.taxon_namespace,
        )
    elif method_token == "quartet":
        if args.starting_tree is None:
            raise ConsensusError("You need to specify the path to (nearly) fully-resolved starting tree with SQDG approach.")
        try:
            starting_tree = Tree_with_support.get(
                path=args.starting_tree,
                schema="newick",
                taxon_namespace=input_trees.taxon_namespace,
            )
        except Exception as e:
            raise ConsensusError(f"Error in reading starting tree: {e}")
        consensus_tree = starting_tree.SQD_greedy_pruning(input_trees)
    else:
        if args.starting_tree is None:
            starting_tree = input_trees.majority_rule_consensus()
        else:
            try:
                starting_tree = Tree_with_support.get(
                    path=args.starting_tree,
                    schema="newick",
                    taxon_namespace=input_trees.taxon_namespace,
                )
            except Exception as e:
                raise ConsensusError(f"Error in reading starting tree: {e}")

        method_dict = {
            "scaled_transfer": STDGreedyConsensus,
            "unscaled_transfer": SUTDGreedyConsensus,
        }

        gre = method_dict[method_token](input_trees)
        gre.specify_initial_tree(starting_tree)
        gre.greedy(method="most")
        consensus_tree = gre.return_current_tree()
        if not isinstance(consensus_tree, Tree_with_support):
            consensus_tree = Tree_with_support(consensus_tree)

    output_path, support_paths = _write_all_outputs(consensus_tree, input_trees, output_file)
    print(f"Output saved to: {output_path}")
    if support_paths:
        print("\nTree files with support values:")
        for path, label in support_paths:
            print(f"- {path} ({label})")
            
if __name__ == "__main__":
    main()
