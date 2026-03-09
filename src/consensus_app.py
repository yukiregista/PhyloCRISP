import os
import logging
import sys
import tempfile
import subprocess
import shutil
import traceback
from pathlib import Path
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QVBoxLayout, QWidget, QLabel, QTextEdit, QPushButton, QComboBox, QMessageBox, 
    QHBoxLayout, QFileDialog, QRadioButton, QButtonGroup, QLineEdit, QCheckBox, QProgressBar
)
from PyQt6.QtCore import QTimer, QThread, pyqtSignal, QStandardPaths
from PyQt6.QtGui import QIcon
import Consensus
from Consensus._app_options import (
    default_output_filename,
    effective_strategy_token,
    method_token_from_gui_label,
)

_ = Consensus.__file__  # Ensure the Consensus package is imported

APP_NAME = "PhyloCRISP"
APP_ICON_RELATIVE_PATH = os.path.join("assets", "icons", "image1.png")


def _startup_log_file() -> str:
    """Return a writable startup log path across dev and installed builds."""
    candidates = []

    if os.name == "nt":
        local_app_data = os.environ.get("LOCALAPPDATA")
        if local_app_data:
            candidates.append(os.path.join(local_app_data, APP_NAME, "logs"))

    qt_app_data = QStandardPaths.writableLocation(QStandardPaths.StandardLocation.AppDataLocation)
    if qt_app_data:
        candidates.append(os.path.join(qt_app_data, "logs"))

    candidates.append(os.path.join(tempfile.gettempdir(), APP_NAME, "logs"))
    candidates.append(os.getcwd())

    for directory in candidates:
        try:
            os.makedirs(directory, exist_ok=True)
            return os.path.join(directory, "fatal_startup.log")
        except OSError:
            continue

    return "fatal_startup.log"


def _is_writable_directory(path: str) -> bool:
    """Check whether a directory exists (or can be created) and is writable."""
    if not path:
        return False
    try:
        os.makedirs(path, exist_ok=True)
        fd, probe_path = tempfile.mkstemp(prefix=".__phylo_write_test_", dir=path)
        os.close(fd)
        os.remove(probe_path)
        return True
    except OSError:
        return False


def _find_app_icon_path() -> str | None:
    """Locate app icon path in dev and PyInstaller layouts."""
    candidates: list[Path] = []

    if getattr(sys, "frozen", False):
        meipass = getattr(sys, "_MEIPASS", "")
        if meipass:
            candidates.append(Path(meipass) / APP_ICON_RELATIVE_PATH)
        exe_dir = Path(sys.executable).resolve().parent
        candidates.append(exe_dir / APP_ICON_RELATIVE_PATH)
        candidates.append(exe_dir / "_internal" / APP_ICON_RELATIVE_PATH)

    script_dir = Path(__file__).resolve().parent
    candidates.append(script_dir.parent / APP_ICON_RELATIVE_PATH)
    candidates.append(script_dir / APP_ICON_RELATIVE_PATH)

    for candidate in candidates:
        if candidate.exists():
            return os.fspath(candidate)
    return None


# Setup early logging (must never crash app startup)
try:
    logging.basicConfig(
        filename=_startup_log_file(),
        level=logging.DEBUG,
        format="%(asctime)s - %(levelname)s - %(message)s"
    )
except OSError:
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s - %(levelname)s - %(message)s"
    )
logging.debug("===== Application startup initiated =====")

# File size threshold for auto-hiding content (in bytes)
LARGE_FILE_THRESHOLD = 10 * 1024  # 10KB
LARGE_OUTPUT_PREVIEW_THRESHOLD = 100 * 1024  # 100KB


def _maybe_run_embedded_cli_mode():
    """Run consensus CLI logic when launched with internal CLI flag."""
    flag = "--run-consensus-cli"
    if flag not in sys.argv:
        return
    idx = sys.argv.index(flag)
    cli_args = sys.argv[idx + 1:]
    from Consensus.main import main as consensus_main
    try:
        consensus_main(cli_args)
        raise SystemExit(0)
    except SystemExit:
        raise
    except Exception:
        traceback.print_exc()
        raise SystemExit(1)


def _build_consensus_command(args):
    """Build command for running CLI logic from GUI worker."""
    if getattr(sys, "frozen", False):
        # Re-enter this same executable in internal CLI mode.
        return [sys.executable, "--run-consensus-cli"] + args
    # Development mode: run the module directly from the current interpreter.
    return [sys.executable, "-m", "Consensus.main"] + args


def _cli_popen_window_kwargs():
    """Return subprocess kwargs to keep the CLI hidden on Windows."""
    if os.name != "nt":
        return {}
    startupinfo = subprocess.STARTUPINFO()
    startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW
    startupinfo.wShowWindow = 0
    return {
        "creationflags": subprocess.CREATE_NO_WINDOW,
        "startupinfo": startupinfo,
    }


_maybe_run_embedded_cli_mode()

class ConsensusWorkerThread(QThread):
    """Worker thread for running consensus calculations via CLI command."""
    finished = pyqtSignal(str)  # Signal emitted when done
    error = pyqtSignal(str)     # Signal emitted on error
    
    def __init__(self, args, temp_dir):
        super().__init__()
        self.args = args
        self.temp_dir = temp_dir
        self._stop_requested = False
        self.process = None
    
    def stop(self):
        """Request the thread to stop."""
        self._stop_requested = True
        if self.process:
            self.process.terminate()
        self.requestInterruption()
    
    def run(self):
        try:
            if self._stop_requested:
                return
            
            cmd = _build_consensus_command(self.args)
            
            logging.info(f"Running command: {' '.join(cmd)}")
            
            
            self.process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                **_cli_popen_window_kwargs(),
            )
            
            stdout, stderr = self.process.communicate()
            
            if self._stop_requested:
                return
                
            if self.process.returncode == 0:
                logging.info(f"Command completed successfully. stdout: {stdout}")
                self.finished.emit("success")
            else:
                error_msg = f"Command failed with return code {self.process.returncode}"
                if stderr:
                    error_msg += f"\nError: {stderr}"
                if stdout:
                    error_msg += f"\nOutput: {stdout}"
                logging.error(error_msg)
                self.error.emit(error_msg)
                
        except FileNotFoundError:
            error_msg = "Consensus execution command not found. Reinstall PhyloCRISP."
            logging.error(error_msg)
            self.error.emit(error_msg)
        except Exception as e:
            if not self._stop_requested:
                error_msg = f"Error running consensus: {str(e)}"
                logging.error(error_msg)
                self.error.emit(error_msg)

class ConsensusApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("PhyloCRISP")
        icon_path = _find_app_icon_path()
        if icon_path:
            self.setWindowIcon(QIcon(icon_path))
        self.setGeometry(100, 100, 700, 850)

        # Store file paths and content
        self.input_trees_content = ""
        self.input_trees_filename = ""
        self.input_trees_file_path = ""
        self.starting_tree_content = ""
        self.starting_tree_filename = ""
        self.starting_tree_file_path = ""

        # Output state
        self.output_content = ""
        self.output_temp_file_path = ""
        self.output_size_bytes = 0
        self.output_preview_hidden = False
        self.output_saved = False
        self.last_saved_output_path = ""
        self._last_default_output_filename = ""
        self.last_browse_dir = ""
        
        # Worker thread
        self.worker_thread = None
        self.temp_dir_path = None

        self.setup_ui()

    def setup_ui(self):
        """Set up the user interface."""
        # Main layout
        layout = QVBoxLayout()

        # Input trees section
        self.input_label = QLabel("Input Trees:")
        layout.addWidget(self.input_label)

        # Input trees display option
        input_display_layout = QHBoxLayout()
        self.input_display_checkbox = QCheckBox("Display content")
        self.input_display_checkbox.setChecked(True)
        self.input_display_checkbox.toggled.connect(self.toggle_input_display)
        input_display_layout.addWidget(self.input_display_checkbox)
        
        # Add file size info label
        self.input_size_label = QLabel("")
        input_display_layout.addWidget(self.input_size_label)
        input_display_layout.addStretch()
        layout.addLayout(input_display_layout)

        # Input trees text area
        self.input_text = QTextEdit()
        self.input_text.textChanged.connect(self.on_input_text_changed)
        layout.addWidget(self.input_text)

        # File browsing layout for input trees
        input_file_layout = QHBoxLayout()
        self.input_file_button = QPushButton("Browse Input Trees File")
        self.input_file_button.clicked.connect(self.browse_input_file)
        input_file_layout.addWidget(self.input_file_button)
        layout.addLayout(input_file_layout)

        # Method selection
        self.method_label = QLabel("Consensus Method:")
        layout.addWidget(self.method_label)
        self.method_combo = QComboBox()
        self.method_combo.addItems(["Scaled Transfer", "Unscaled Transfer", "Quartet", "Majority Rule"])
        self.method_combo.currentTextChanged.connect(self.on_method_changed)
        layout.addWidget(self.method_combo)

        # Optimization strategy selection
        self.strategy_label = QLabel("Optimization Strategy:")
        layout.addWidget(self.strategy_label)

        self.add_and_prune_radio = QRadioButton("Add and Prune")
        self.add_and_prune_radio.setChecked(True)
        self.prune_radio = QRadioButton("Prune Only")
        self.strategy_group = QButtonGroup()
        self.strategy_group.addButton(self.add_and_prune_radio)
        self.strategy_group.addButton(self.prune_radio)
        self.add_and_prune_radio.toggled.connect(self.on_strategy_changed)
        self.prune_radio.toggled.connect(self.on_strategy_changed)

        layout.addWidget(self.add_and_prune_radio)
        layout.addWidget(self.prune_radio)

        # Starting tree section
        self.starting_tree_label = QLabel("Starting Tree (Optional for Add and Prune, Required for Prune Only):")
        layout.addWidget(self.starting_tree_label)

        # Starting tree display option
        starting_display_layout = QHBoxLayout()
        self.starting_display_checkbox = QCheckBox("Display content")
        self.starting_display_checkbox.setChecked(True)
        self.starting_display_checkbox.toggled.connect(self.toggle_starting_display)
        starting_display_layout.addWidget(self.starting_display_checkbox)
        
        # Add file size info label
        self.starting_size_label = QLabel("")
        starting_display_layout.addWidget(self.starting_size_label)
        starting_display_layout.addStretch()
        layout.addLayout(starting_display_layout)

        # Starting tree text area
        self.starting_tree_text = QTextEdit()
        self.starting_tree_text.textChanged.connect(self.on_starting_text_changed)
        layout.addWidget(self.starting_tree_text)

        # File browsing layout for starting tree
        starting_file_layout = QHBoxLayout()
        self.starting_file_button = QPushButton("Browse Starting Tree File")
        self.starting_file_button.clicked.connect(self.browse_starting_file)
        starting_file_layout.addWidget(self.starting_file_button)
        layout.addLayout(starting_file_layout)

        # Output filename section
        output_filename_layout = QHBoxLayout()
        self.output_filename_label = QLabel("Output Filename:")
        self.output_filename_input = QLineEdit()
        self.output_filename_input.setPlaceholderText("Enter output filename (optional)")
        self.output_filename_input.textChanged.connect(self.on_output_filename_changed)
        output_filename_layout.addWidget(self.output_filename_label)
        output_filename_layout.addWidget(self.output_filename_input)
        self.auto_save_output_checkbox = QCheckBox("Auto-save output")
        self.auto_save_output_checkbox.setToolTip(
            "Automatically save generated output to the input file folder (or Downloads if no input file path is available)."
        )
        output_filename_layout.addWidget(self.auto_save_output_checkbox)
        layout.addLayout(output_filename_layout)

        # Output directory section
        output_dir_layout = QHBoxLayout()
        self.output_dir_label = QLabel("Output Directory:")
        self.output_dir_input = QLineEdit()
        self.output_dir_input.setPlaceholderText("Optional (default: input folder / Downloads)")
        self.output_dir_input.textChanged.connect(self.on_output_directory_changed)
        self.output_dir_button = QPushButton("Browse Folder")
        self.output_dir_button.clicked.connect(self.browse_output_directory)
        output_dir_layout.addWidget(self.output_dir_label)
        output_dir_layout.addWidget(self.output_dir_input)
        output_dir_layout.addWidget(self.output_dir_button)
        layout.addLayout(output_dir_layout)

        # Generate and Stop button layout
        button_layout = QHBoxLayout()
        self.generate_button = QPushButton("Generate Consensus Tree")
        self.generate_button.clicked.connect(self.generate_consensus_tree)
        button_layout.addWidget(self.generate_button)
        
        self.stop_button = QPushButton("Stop")
        self.stop_button.clicked.connect(self.stop_consensus_generation)
        self.stop_button.setEnabled(False)
        self.stop_button.setStyleSheet("QPushButton { background-color: #ff6b6b; color: white; }")
        button_layout.addWidget(self.stop_button)
        
        layout.addLayout(button_layout)

        # Progress bar and status
        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        layout.addWidget(self.progress_bar)

        self.status_label = QLabel("")
        self.status_label.setVisible(False)
        layout.addWidget(self.status_label)

        # Output section
        self.output_label = QLabel("Consensus Tree:")
        layout.addWidget(self.output_label)

        # Output display option
        output_display_layout = QHBoxLayout()
        self.output_display_checkbox = QCheckBox("Display content")
        self.output_display_checkbox.setChecked(True)
        self.output_display_checkbox.toggled.connect(self.toggle_output_display)
        output_display_layout.addWidget(self.output_display_checkbox)
        
        # Save output button
        self.save_output_button = QPushButton("Save Output")
        self.save_output_button.clicked.connect(self.save_output_file)
        self.save_output_button.setEnabled(False)
        output_display_layout.addWidget(self.save_output_button)
        output_display_layout.addStretch()
        layout.addLayout(output_display_layout)

        self.output_path_label = QLabel("Output location: not saved yet")
        self.output_path_label.setWordWrap(True)
        layout.addWidget(self.output_path_label)

        # Output text area
        self.output_text = QTextEdit()
        self.output_text.setReadOnly(True)
        layout.addWidget(self.output_text)

        # Set central widget
        container = QWidget()
        container.setLayout(layout)
        self.setCentralWidget(container)
        self.on_method_changed(self.method_combo.currentText())
        self.refresh_output_filename_placeholder()

    def on_input_text_changed(self):
        """Called when input text is edited."""
        # Only update content if we're showing actual content (not filename)
        if self.input_display_checkbox.isChecked():
            self.input_trees_content = self.input_text.toPlainText()

    def on_starting_text_changed(self):
        """Called when starting tree text is edited."""
        # Only update content if we're showing actual content (not filename)
        if self.starting_display_checkbox.isChecked():
            self.starting_tree_content = self.starting_tree_text.toPlainText()
        self.refresh_output_filename_placeholder()

    def format_file_size(self, size_bytes):
        """Format file size in human-readable format."""
        if size_bytes < 1024:
            return f"{size_bytes} B"
        elif size_bytes < 1024 * 1024:
            return f"{size_bytes / 1024:.1f} KB"
        else:
            return f"{size_bytes / (1024 * 1024):.1f} MB"

    def on_output_filename_changed(self, _text=None):
        """Refresh output status text when the filename preset changes."""
        if self.has_output_generated():
            self.refresh_output_view()

    def on_output_directory_changed(self, _text=None):
        """Refresh output status text when the output directory changes."""
        if self.has_output_generated():
            self.refresh_output_view()

    def _current_method_token(self):
        return method_token_from_gui_label(self.method_combo.currentText())

    def _current_strategy_token(self):
        requested = "prune" if self.prune_radio.isChecked() else "add_and_prune"
        return effective_strategy_token(self._current_method_token(), requested)

    def get_default_output_filename(self):
        """Build the default output filename from method/strategy/starting-tree settings."""
        method_token = self._current_method_token()
        strategy_token = self._current_strategy_token()
        include_starting = method_token != "mr" and bool(self.starting_tree_content.strip())
        if self.starting_tree_filename:
            init_name = os.path.splitext(self.starting_tree_filename)[0]
        elif self.starting_tree_file_path:
            init_name = os.path.splitext(os.path.basename(self.starting_tree_file_path))[0]
        else:
            init_name = "starting_tree"
        return default_output_filename(
            method_token=method_token,
            strategy_token=strategy_token,
            include_starting_tree=include_starting,
            starting_tree_name=init_name,
        )

    def refresh_output_filename_placeholder(self):
        default_name = self.get_default_output_filename()
        self.output_filename_input.setPlaceholderText(f"Optional (default: {default_name})")
        current_name = self.output_filename_input.text().strip()
        if not current_name or current_name == self._last_default_output_filename:
            if current_name != default_name:
                self.output_filename_input.setText(default_name)
        self._last_default_output_filename = default_name
        if self.has_output_generated():
            self.refresh_output_view()

    def has_output_generated(self):
        """Return True if a consensus output has been generated in this session."""
        return bool(self.output_temp_file_path or self.output_content)

    def on_method_changed(self, method):
        """Update strategy availability according to selected consensus method."""
        if method == "Quartet":
            self.add_and_prune_radio.setEnabled(False)
            self.prune_radio.setEnabled(True)
            self.prune_radio.setChecked(True)
        elif method == "Majority Rule":
            # Majority-rule has no optimization strategy choice.
            self.add_and_prune_radio.setChecked(True)
            self.add_and_prune_radio.setEnabled(False)
            self.prune_radio.setEnabled(False)
        else:
            self.add_and_prune_radio.setEnabled(True)
            self.prune_radio.setEnabled(True)
        self.refresh_output_filename_placeholder()

    def on_strategy_changed(self, _checked=False):
        self.refresh_output_filename_placeholder()

    def cleanup_temp_dir(self):
        """Remove the previous temporary directory, if present."""
        if self.temp_dir_path and os.path.isdir(self.temp_dir_path):
            try:
                shutil.rmtree(self.temp_dir_path)
            except Exception as e:
                logging.warning(f"Failed to clean temporary directory {self.temp_dir_path}: {e}")
        self.temp_dir_path = None

    def reset_output_state(self):
        """Clear output state for a new run."""
        self.output_content = ""
        self.output_temp_file_path = ""
        self.output_size_bytes = 0
        self.output_preview_hidden = False
        self.output_saved = False
        self.last_saved_output_path = ""
        if hasattr(self, "save_output_button"):
            self.save_output_button.setEnabled(False)
        if hasattr(self, "output_path_label"):
            self.output_path_label.setText("Output location: not saved yet")
        if hasattr(self, "output_text"):
            self.refresh_output_view()

    def get_default_output_directory(self):
        """Choose a default directory for saving output files."""
        candidates = []

        user_selected = self.output_dir_input.text().strip()
        if user_selected:
            if os.path.isabs(user_selected):
                candidates.append(user_selected)
            else:
                candidates.append(os.path.abspath(user_selected))

        if self.input_trees_file_path:
            input_dir = os.path.dirname(self.input_trees_file_path)
            if input_dir:
                candidates.append(input_dir)

        downloads_dir = QStandardPaths.writableLocation(QStandardPaths.StandardLocation.DownloadLocation)
        if downloads_dir:
            candidates.append(downloads_dir)

        documents_dir = QStandardPaths.writableLocation(QStandardPaths.StandardLocation.DocumentsLocation)
        if documents_dir:
            candidates.append(documents_dir)

        fallback_downloads = os.path.join(os.path.expanduser("~"), "Downloads")
        candidates.append(fallback_downloads)

        if os.name == "nt":
            local_app_data = os.environ.get("LOCALAPPDATA")
            if local_app_data:
                candidates.append(os.path.join(local_app_data, APP_NAME, "output"))

        candidates.append(os.path.join(tempfile.gettempdir(), APP_NAME, "output"))
        candidates.append(tempfile.gettempdir())
        candidates.append(os.getcwd())

        for candidate in candidates:
            if _is_writable_directory(candidate):
                return candidate

        return tempfile.gettempdir()

    def browse_output_directory(self):
        """Open a folder dialog to choose the default output directory."""
        initial_dir = self.get_default_dialog_directory()
        entered_dir = self.output_dir_input.text().strip()
        if (not self.last_browse_dir) and entered_dir and os.path.isdir(entered_dir):
            initial_dir = entered_dir
        selected_dir = QFileDialog.getExistingDirectory(self, "Select Output Directory", initial_dir)
        if not selected_dir:
            return
        self._remember_browse_path(selected_dir)
        self.output_dir_input.setText(selected_dir)

    def get_default_dialog_directory(self):
        """Return a readable home-like directory for file/folder pickers."""
        candidates = [
            self.last_browse_dir,
            QStandardPaths.writableLocation(QStandardPaths.StandardLocation.HomeLocation),
            os.path.expanduser("~"),
            os.getcwd(),
        ]
        for candidate in candidates:
            if candidate and os.path.isdir(candidate):
                return candidate
        return os.getcwd()

    def _remember_browse_path(self, path):
        """Remember the last directory used in load/save dialogs."""
        if not path:
            return
        if os.path.isdir(path):
            directory = path
        else:
            directory = os.path.dirname(path)
        if directory and os.path.isdir(directory):
            self.last_browse_dir = directory

    def _resolve_user_output_path(self, filename_text):
        """Resolve a user-entered output filename to a filesystem path."""
        if os.path.isabs(filename_text):
            return filename_text
        if os.path.dirname(filename_text):
            return os.path.abspath(filename_text)
        return os.path.join(self.get_default_output_directory(), filename_text)

    def _make_unique_path(self, file_path):
        """Append a numeric suffix if a path already exists."""
        if not os.path.exists(file_path):
            return file_path

        base, ext = os.path.splitext(file_path)
        counter = 1
        candidate = f"{base}_{counter}{ext}"
        while os.path.exists(candidate):
            counter += 1
            candidate = f"{base}_{counter}{ext}"
        return candidate

    def get_preset_output_path(self, auto_save=False):
        """Return the output path derived from the filename field, or a default path."""
        filename_text = self.output_filename_input.text().strip()
        if not filename_text:
            if not auto_save:
                return None
            filename_text = self.get_default_output_filename()

        file_path = self._resolve_user_output_path(filename_text)
        if auto_save:
            file_path = self._make_unique_path(file_path)
        return file_path

    def get_default_save_dialog_path(self):
        """Build a sensible initial path for the save dialog."""
        preset_path = self.get_preset_output_path(auto_save=False)
        if preset_path:
            return preset_path
        return os.path.join(self.get_default_dialog_directory(), self.get_default_output_filename())

    def load_output_preview_if_small(self):
        """Load output text into memory only when the file is small enough."""
        if self.output_content:
            return True

        if not self.output_temp_file_path or not os.path.exists(self.output_temp_file_path):
            return False

        if self.output_size_bytes > LARGE_OUTPUT_PREVIEW_THRESHOLD:
            self.output_preview_hidden = True
            return False

        try:
            with open(self.output_temp_file_path, "r") as f:
                self.output_content = f.read()
            self.output_preview_hidden = False
            return True
        except Exception as e:
            logging.error(f"Failed to load output preview: {e}")
            return False

    def refresh_output_view(self):
        """Refresh the output panel according to current display mode and output state."""
        if self.output_saved and self.last_saved_output_path:
            self.output_path_label.setText(f"Output location: {self.last_saved_output_path}")
        elif self.has_output_generated():
            preset_path = self.get_preset_output_path(auto_save=False)
            if preset_path:
                self.output_path_label.setText(f"Output location (default): {preset_path}")
            else:
                self.output_path_label.setText("Output location: choose with 'Save Output'")
        else:
            self.output_path_label.setText("Output location: not saved yet")

        if self.output_display_checkbox.isChecked():
            if not self.has_output_generated():
                self.output_text.setPlainText("")
                return

            self.load_output_preview_if_small()

            if self.output_preview_hidden:
                lines = [
                    "Output generated."
                ]
                if self.output_size_bytes:
                    lines[0] = f"Output generated ({self.format_file_size(self.output_size_bytes)})."
                lines.append("Preview hidden for performance (large output).")
                if self.output_saved and self.last_saved_output_path:
                    lines.append(f"Saved to: {self.last_saved_output_path}")
                else:
                    lines.append("Click 'Save Output' to export the Newick file.")
                    preset_path = self.get_preset_output_path(auto_save=False)
                    if preset_path:
                        lines.append(f"Preset save path: {preset_path}")
                self.output_text.setPlainText("\n".join(lines))
                return

            if self.output_content:
                self.output_text.setPlainText(self.output_content)
            else:
                self.output_text.setPlainText(
                    "Output generated, but preview is unavailable. Use 'Save Output' to export."
                )
            return

        if self.output_saved and self.last_saved_output_path:
            msg = f"Output saved to: {self.last_saved_output_path}"
            if self.output_size_bytes:
                msg += f" ({self.format_file_size(self.output_size_bytes)})"
            self.output_text.setPlainText(msg)
            return

        if self.has_output_generated():
            msg = "Output generated"
            if self.output_size_bytes:
                msg += f" ({self.format_file_size(self.output_size_bytes)})"
            preset_path = self.get_preset_output_path(auto_save=False)
            if preset_path:
                msg += f"\nPreset save path: {preset_path}"
            else:
                msg += "\nUse 'Save Output' to choose a location."
            self.output_text.setPlainText(msg)
            return

        self.output_text.setPlainText("No output generated")

    def save_output_to_path(self, file_path, show_message=True):
        """Save output to the specified path, preferring file copy from temp output."""
        if not self.has_output_generated():
            if show_message:
                QMessageBox.warning(self, "Warning", "No output to save.")
            return False

        try:
            self._remember_browse_path(file_path)
            if self.output_temp_file_path and os.path.exists(self.output_temp_file_path):
                shutil.copyfile(self.output_temp_file_path, file_path)
            else:
                with open(file_path, "w") as f:
                    f.write(self.output_content)

            support_copies = []
            support_suffix_labels = {
                ".fsupp": "frequency support",
                ".tsupp": "transfer support",
            }
            if self.output_temp_file_path:
                for suffix, label in support_suffix_labels.items():
                    src = self.output_temp_file_path + suffix
                    dst = file_path + suffix
                    if os.path.exists(src):
                        shutil.copyfile(src, dst)
                        support_copies.append((dst, label))

            self.output_saved = True
            self.last_saved_output_path = file_path
            self.output_filename_input.setText(os.path.basename(file_path))
            self.refresh_output_view()

            if show_message:
                message = f"Output saved to: {file_path}"
                if support_copies:
                    message += "\n\nTree files with support values:"
                    for copied, label in support_copies:
                        message += f"\n- {copied} ({label})"
                QMessageBox.information(self, "File Saved", message)

            logging.info(f"Saved output to: {file_path}")
            return True

        except Exception as e:
            if show_message:
                QMessageBox.critical(self, "Error", f"Failed to save file: {str(e)}")
            logging.error(f"Failed to save output file {file_path}: {e}")
            return False

    def toggle_input_display(self, checked):
        """Toggle between displaying content or filename for input trees."""
        if checked:
            # Show content and make editable
            self.input_text.setPlainText(self.input_trees_content)
            self.input_text.setReadOnly(False)
        else:
            # Show filename and make non-editable
            if self.input_trees_filename:
                self.input_text.setPlainText(f"File loaded: {self.input_trees_filename}")
            else:
                self.input_text.setPlainText("No file loaded")
            self.input_text.setReadOnly(True)

    def toggle_starting_display(self, checked):
        """Toggle between displaying content or filename for starting tree."""
        if checked:
            # Show content and make editable
            self.starting_tree_text.setPlainText(self.starting_tree_content)
            self.starting_tree_text.setReadOnly(False)
        else:
            # Show filename and make non-editable
            if self.starting_tree_filename:
                self.starting_tree_text.setPlainText(f"File loaded: {self.starting_tree_filename}")
            else:
                self.starting_tree_text.setPlainText("No file loaded")
            self.starting_tree_text.setReadOnly(True)

    def toggle_output_display(self, checked):
        """Toggle between displaying content or filename for output."""
        self.refresh_output_view()

    def browse_input_file(self):
        """Open a file dialog to browse for a file containing Newick trees."""
        initial_dir = self.get_default_dialog_directory()
        if (not self.last_browse_dir) and self.input_trees_file_path and os.path.isfile(self.input_trees_file_path):
            initial_dir = os.path.dirname(self.input_trees_file_path)
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Open File",
            initial_dir,
            "All Files (*);;Text Files (*.txt)",
        )
        if file_path:
            try:
                self._remember_browse_path(file_path)
                # Get file size
                file_size = os.path.getsize(file_path)
                is_large_file = file_size > LARGE_FILE_THRESHOLD
                
                with open(file_path, "r") as file:
                    self.input_trees_content = file.read()
                    self.input_trees_filename = os.path.basename(file_path)
                    self.input_trees_file_path = file_path
                    
                    # Update size label
                    self.input_size_label.setText(f"({self.format_file_size(file_size)})")
                    
                    # Auto-hide large files
                    if is_large_file:
                        self.input_display_checkbox.setChecked(False)
                        self.input_text.setPlainText(f"File loaded: {self.input_trees_filename} (Large file - content hidden)")
                        self.input_text.setReadOnly(True)
                    else:
                        if self.input_display_checkbox.isChecked():
                            self.input_text.setPlainText(self.input_trees_content)
                            self.input_text.setReadOnly(False)
                        else:
                            self.input_text.setPlainText(f"File loaded: {self.input_trees_filename}")
                            self.input_text.setReadOnly(True)
                    self.refresh_output_filename_placeholder()
                    logging.info(f"Loaded file: {file_path} ({self.format_file_size(file_size)})")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to load file: {str(e)}")
                logging.error(f"Failed to load file: {str(e)}")

    def browse_starting_file(self):
        """Open a file dialog to browse for a starting tree file."""
        initial_dir = self.get_default_dialog_directory()
        if (not self.last_browse_dir) and self.starting_tree_file_path and os.path.isfile(self.starting_tree_file_path):
            initial_dir = os.path.dirname(self.starting_tree_file_path)
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Open Starting Tree File",
            initial_dir,
            "All Files (*);;Text Files (*.txt)",
        )
        if file_path:
            try:
                self._remember_browse_path(file_path)
                # Get file size
                file_size = os.path.getsize(file_path)
                is_large_file = file_size > LARGE_FILE_THRESHOLD
                
                with open(file_path, "r") as file:
                    self.starting_tree_content = file.read()
                    self.starting_tree_filename = os.path.basename(file_path)
                    self.starting_tree_file_path = file_path
                    
                    # Update size label
                    self.starting_size_label.setText(f"({self.format_file_size(file_size)})")
                    
                    # Auto-hide large files
                    if is_large_file:
                        self.starting_display_checkbox.setChecked(False)
                        self.starting_tree_text.setPlainText(f"File loaded: {self.starting_tree_filename} (Large file - content hidden)")
                        self.starting_tree_text.setReadOnly(True)
                    else:
                        if self.starting_display_checkbox.isChecked():
                            self.starting_tree_text.setPlainText(self.starting_tree_content)
                            self.starting_tree_text.setReadOnly(False)
                        else:
                            self.starting_tree_text.setPlainText(f"File loaded: {self.starting_tree_filename}")
                            self.starting_tree_text.setReadOnly(True)
                    self.refresh_output_filename_placeholder()
                    logging.info(f"Loaded starting tree file: {file_path} ({self.format_file_size(file_size)})")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to load file: {str(e)}")
                logging.error(f"Failed to load file: {str(e)}")


    def save_output_file(self):
        """Save the output to a file."""
        if not self.has_output_generated():
            QMessageBox.warning(self, "Warning", "No output to save.")
            return
            
        output_filename = self.output_filename_input.text().strip()
        if not output_filename:
            # Suggest .nwk (Newick) as default but allow user to choose
            file_path, _ = QFileDialog.getSaveFileName(
                self, 
                "Save Output File", 
                self.get_default_save_dialog_path(),
                "Newick Files (*.nwk);;Tree Files (*.tre);;Text Files (*.txt);;All Files (*)"
            )
            if not file_path:
                return
            self._remember_browse_path(file_path)
        else:
            file_path = self.get_preset_output_path(auto_save=False)
            # Don't force any extension - let user decide

        self.save_output_to_path(file_path, show_message=True)


    def show_progress(self, message):
        """Show progress bar and status message."""
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)  # Indeterminate progress
        self.status_label.setText(message)
        self.status_label.setVisible(True)
        self.generate_button.setEnabled(False)
        self.stop_button.setEnabled(True)

    def hide_progress(self):
        """Hide progress bar and status message."""
        self.progress_bar.setVisible(False)
        self.status_label.setVisible(False)
        self.generate_button.setEnabled(True)
        self.stop_button.setEnabled(False)

    def stop_consensus_generation(self):
        """Stop the consensus generation process."""
        if self.worker_thread and self.worker_thread.isRunning():
            self.worker_thread.stop()
            self.status_label.setText("Stopping process...")
            QTimer.singleShot(3000, self.force_cleanup)  # Force cleanup after 3 seconds

    def force_cleanup(self):
        """Force cleanup if thread doesn't stop gracefully."""
        if self.worker_thread and self.worker_thread.isRunning():
            self.worker_thread.terminate()
            self.worker_thread.wait(1000)  # Wait up to 1 second
        self.cleanup_after_stop()

    def cleanup_after_stop(self):
        """Clean up after stopping the process."""
        self.hide_progress()
        self.cleanup_temp_dir()
        self.worker_thread = None
        self.status_label.setText("Process stopped by user.")
        self.status_label.setVisible(True)
        QTimer.singleShot(3000, lambda: self.status_label.setVisible(False))

    def generate_consensus_tree(self):
        # Use stored content instead of text area content
        input_trees = self.input_trees_content
        method = self.method_combo.currentText()
        method_token = self._current_method_token()
        strategy = self._current_strategy_token()
        starting_tree = self.starting_tree_content

        if not input_trees.strip():
            QMessageBox.warning(self, "Error", "Input trees are required.")
            return

        if method != "Majority Rule" and strategy == "prune" and not starting_tree.strip():
            QMessageBox.warning(self, "Error", "Starting tree is required for the Prune Only strategy.")
            return

        # Reset previous output state and temp files
        self.cleanup_temp_dir()
        self.reset_output_state()

        # Create a persistent temporary directory
        self.temp_dir_path = tempfile.mkdtemp()
        input_file = os.path.join(self.temp_dir_path, "input_trees.txt")
        output_file = os.path.join(self.temp_dir_path, "output_consensus_tree")
        starting_tree_file = os.path.join(self.temp_dir_path, "starting_tree.txt")

        try:
            with open(input_file, "w") as f:
                f.write(input_trees)

            if starting_tree.strip():
                with open(starting_tree_file, "w") as f:
                    f.write(starting_tree)

            args = [
                "--input-trees", input_file,
                "--method", method_token,
                "--output-filename", "output_consensus_tree.nwk",
                "--output-dir", self.temp_dir_path,
            ]

            if starting_tree.strip():
                args.extend(["--starting-tree", starting_tree_file])

            if strategy:
                args.extend(["--strategy", strategy])

            # Show progress
            self.show_progress(f"Generating consensus tree using {method}...")

            # Start the worker thread
            self.worker_thread = ConsensusWorkerThread(args, self.temp_dir_path)
            self.worker_thread.finished.connect(self.on_consensus_finished)
            self.worker_thread.error.connect(self.on_consensus_error)
            self.worker_thread.start()

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to prepare input files: {str(e)}")
            self.hide_progress()

    def on_consensus_finished(self, result):
        """Handle successful completion of consensus calculation."""
        self.hide_progress()
        
        if result == "success":
            try:
                output_file = os.path.join(self.temp_dir_path, "output_consensus_tree.nwk")
                self.output_temp_file_path = output_file
                self.output_size_bytes = os.path.getsize(self.output_temp_file_path)
                self.output_content = ""
                self.output_preview_hidden = self.output_size_bytes > LARGE_OUTPUT_PREVIEW_THRESHOLD
                self.output_saved = False
                self.last_saved_output_path = ""
                
                # Enable save button
                self.save_output_button.setEnabled(True)

                if self.auto_save_output_checkbox.isChecked():
                    auto_path = self.get_preset_output_path(auto_save=True)
                    auto_saved = self.save_output_to_path(auto_path, show_message=True)
                    if not auto_saved:
                        QMessageBox.warning(
                            self,
                            "Auto-save failed",
                            "Output was generated, but automatic saving failed. You can still use 'Save Output'."
                        )

                self.refresh_output_view()
                        
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to read output: {str(e)}")
        
        # Clean up
        self.worker_thread = None

    def on_consensus_error(self, error_message):
        """Handle error in consensus calculation."""
        self.hide_progress()
        QMessageBox.critical(self, "Error", error_message)
        self.worker_thread = None


if __name__ == "__main__":
    logging.info("Starting Consensus Tree Generator...")
    app = QApplication(sys.argv)
    window = ConsensusApp()
    window.show()
    sys.exit(app.exec())
