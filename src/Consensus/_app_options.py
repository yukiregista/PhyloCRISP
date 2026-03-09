import os
from typing import Optional


METHOD_TOKEN_TO_GUI_LABEL = {
    "scaled_transfer": "Scaled Transfer",
    "unscaled_transfer": "Unscaled Transfer",
    "quartet": "Quartet",
    "mr": "Majority Rule",
}

GUI_LABEL_TO_METHOD_TOKEN = {v: k for k, v in METHOD_TOKEN_TO_GUI_LABEL.items()}

METHOD_ALIASES = {
    "scaled_transfer": "scaled_transfer",
    "unscaled_transfer": "unscaled_transfer",
    "quartet": "quartet",
    "mr": "mr",
}

STRATEGY_ALIASES = {
    "add_and_prune": "add_and_prune",
    "prune": "prune",
}


def slugify_token(text: str) -> str:
    token = "".join(ch.lower() if ch.isalnum() else "_" for ch in text.strip())
    while "__" in token:
        token = token.replace("__", "_")
    return token.strip("_") or "value"


def normalize_method_token(value: str) -> str:
    key = value.strip().lower()
    if key in METHOD_ALIASES:
        return METHOD_ALIASES[key]
    raise ValueError(
        "Unknown method. Use one of: scaled_transfer, unscaled_transfer, quartet, mr."
    )


def normalize_strategy_token(value: Optional[str]) -> str:
    if value is None or not str(value).strip():
        return "add_and_prune"
    key = value.strip().lower()
    if key in STRATEGY_ALIASES:
        return STRATEGY_ALIASES[key]
    raise ValueError("Unknown strategy. Use one of: add_and_prune, prune.")


def effective_strategy_token(method_token: str, strategy_token: Optional[str]) -> str:
    if method_token == "mr":
        return ""
    normalized_strategy = normalize_strategy_token(strategy_token)
    if method_token == "quartet":
        return "prune"
    return normalized_strategy


def method_token_from_gui_label(label: str) -> str:
    return GUI_LABEL_TO_METHOD_TOKEN.get(label, "scaled_transfer")


def default_output_filename(
    method_token: str,
    strategy_token: Optional[str],
    include_starting_tree: bool = False,
    starting_tree_name: Optional[str] = None,
) -> str:
    parts = [slugify_token(method_token)]

    effective_strategy = effective_strategy_token(method_token, strategy_token)
    if effective_strategy:
        parts.append("addprune" if effective_strategy == "add_and_prune" else "prune")

    if include_starting_tree and method_token != "mr":
        init_name = starting_tree_name or "starting_tree"
        parts.append(f"init-{slugify_token(init_name)}")

    return "_".join(parts) + ".nwk"


def default_output_directory(input_trees_path: str, output_dir: Optional[str] = None) -> str:
    if output_dir:
        return os.path.abspath(output_dir)
    if input_trees_path:
        input_dir = os.path.dirname(os.path.abspath(input_trees_path))
        if input_dir:
            return input_dir
    return os.getcwd()


def resolve_output_path(
    input_trees_path: str,
    method_token: str,
    strategy_token: Optional[str],
    output_filename: Optional[str] = None,
    output_dir: Optional[str] = None,
    starting_tree_path: Optional[str] = None,
) -> str:
    if output_filename:
        if os.path.isabs(output_filename):
            return output_filename
        if os.path.dirname(output_filename):
            return os.path.abspath(output_filename)
        base_dir = default_output_directory(input_trees_path, output_dir)
        return os.path.join(base_dir, output_filename)

    starting_name = None
    include_starting = False
    if method_token != "mr" and starting_tree_path:
        include_starting = True
        starting_name = os.path.splitext(os.path.basename(starting_tree_path))[0]

    default_name = default_output_filename(
        method_token=method_token,
        strategy_token=strategy_token,
        include_starting_tree=include_starting,
        starting_tree_name=starting_name,
    )
    return os.path.join(default_output_directory(input_trees_path, output_dir), default_name)
